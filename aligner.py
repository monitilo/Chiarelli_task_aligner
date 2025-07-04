import subprocess
import os
import sys
import pysam
import re
import gzip
import argparse
import resource
import psutil
import time

#%% Functions


def set_memory_limit_gb(max_gb):
    """Limit the memory usage of the script (including child processes)."""
    max_bytes = max_gb * 1024 ** 3
    resource.setrlimit(resource.RLIMIT_AS, (max_bytes, max_bytes))

def parse_args():
    """Collect the imput files"""
    parser = argparse.ArgumentParser(description="Aligner wrapper script using bwa mem")

    parser.add_argument("-r1","--read1", required=True, help="Path to Read 1 FASTQ file (required)")
    parser.add_argument("-r2","--read2", default="", help="Path to Read 2 FASTQ file (optional for paired-end)")
    parser.add_argument("-f","--reference", required=True, help="Path to reference genome FASTA file")
    parser.add_argument("-o","--output_dir", default="results", help="Directory to save output files")
    parser.add_argument("-t","--threads", type=int, default=4, help="Number of threads to use (default: 4)")
    parser.add_argument("--stats", default=None, help="Stats report filename or path")
    parser.add_argument("--keep_intermediates", type=str, choices=["yes", "no"], default="no", help="Keep intermediate SAM/BAM files (yes/no)")

    args = parser.parse_args() 
    # Initialize it and redefine stats in case of None, filename, or path
    if args.stats is None:
        # No stats arg provided, default to output_dir/report.txt
        args.stats = os.path.join(args.output_dir, "report.txt")
    else:
        # stats arg provided, check if it has a directory component
        if os.path.dirname(args.stats) == "":
            # Only filename given, prepend output_dir
            args.stats = os.path.join(args.output_dir, args.stats)
        else:
            # Full or relative path provided, use it
            args.stats = args.stats

    print(f"Stats report will be saved to: {args.stats}")


    # Validate inputs (better be safe than sorry)
    if not os.path.exists(args.read1):
        parser.error(f"Read1 FASTQ file not found: {args.read1}")
    if args.read2 and not os.path.exists(args.read2):
        parser.error(f"Read2 FASTQ file not found: {args.read2}")
    if not os.path.exists(args.reference):
        parser.error(f"Reference genome file not found: {args.reference}")

    # Make sure output directory exists, or create it
    os.makedirs(args.output_dir, exist_ok=True)

    return args


def parse_flagstat(flagstat_path, output_report_path):
    """ Parse the report (stats) file, saving the useful information.
    Work for the output of samtools flagstat """

    stats = {}
    with open(flagstat_path, "r") as f:
        lines = f.readlines()

    total_reads = None
    mapped_reads = None
    mapped_pct = None

    for line in lines:
        # Total reads
        if "in total" in line:
            m = re.search(r"(\d+) \+ \d+ in total", line)
            if m:
                total_reads = int(m.group(1))
                stats["Total reads"] = total_reads

        # Mapped reads
        elif "mapped (" in line and "primary mapped" not in line:
            m = re.search(r"(\d+) \+ \d+ mapped \(([\d\.]+)%", line)
            if m:
                mapped_reads = int(m.group(1))
                mapped_pct = float(m.group(2))
                stats["Mapped reads"] = f"{mapped_reads} ({mapped_pct}%)"

        # Duplicated reads
        elif "duplicates" in line and "primary duplicates" not in line:
            m = re.search(r"(\d+) \+ \d+ duplicates", line)
            if m and total_reads is not None:
                dup_count = int(m.group(1))
                dup_pct = 100 * dup_count / total_reads
                stats["Duplicated reads"] = f"{dup_count} ({dup_pct:.2f}%)"

        # Singletons
        elif "singletons (" in line:
            m = re.search(r"(\d+) \+ \d+ singletons \(([\d\.]+)%", line)
            if m:
                stats["Singletons"] = f"{m.group(1)} ({m.group(2)}%)"

    # Calculate Unmapped reads and percentage
    if total_reads is not None and mapped_reads is not None and mapped_pct is not None:
        unmapped_reads = total_reads - mapped_reads
        unmapped_pct = 100 - mapped_pct
        stats["Unmapped reads"] = f"{unmapped_reads} ({unmapped_pct:.2f}%)"

    # Write report
    with open(output_report_path, "w") as out:
        for key in ["Total reads", "Mapped reads", "Unmapped reads", "Duplicated reads", "Singletons"]:
            if key in stats:
                out.write(f"{key}: {stats[key]}\n")


def average_base_quality(fastq_path):
    """REQ04, -5 request the average base quality (Phred score) """
    open_func = gzip.open if fastq_path.endswith('.gz') else open
    total_quality = 0
    total_bases = 0
    with open_func(fastq_path, 'rt') as f:  # 'rt' means read text mode even for gzip
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().rstrip()
            plus = f.readline()
            qual = f.readline().rstrip()
            qual_scores = [ord(ch) - 33 for ch in qual]  # Phred+33 decoding
            total_quality += sum(qual_scores)
            total_bases += len(qual_scores)
    return total_quality / total_bases if total_bases > 0 else 0

def average_mapping_quality(bam_path):
    """REQ04, -6 request the average mapping quality (MAPQ) """
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    total_mapq = 0
    count = 0
    for read in bamfile.fetch():
        if not read.is_unmapped:
            total_mapq += read.mapping_quality
            count += 1
    bamfile.close()
    return total_mapq / count if count > 0 else 0

#%% Setting the input

set_memory_limit_gb(16)  # Set limit to 16 GB of RAM

args = parse_args() # Read the imputs

READ1_FASTQ_PATH = args.read1
READ2_FASTQ_PATH = args.read2
REFERENCE_GENOME_PATH = args.reference
RESULTS_DIR = args.output_dir
NUM_THREADS = args.threads

IS_PAIRED_END = bool(READ2_FASTQ_PATH)

OUTPUT_BAM_BASENAME = "aligned_reads.bam"
OUTPUT_STATS_BASENAME = "temp_raw_stats"

UNSORTED_SAM = os.path.join(RESULTS_DIR, "temp_aligned.sam")
UNSORTED_BAM = os.path.join(RESULTS_DIR, "temp_unsorted.bam")
SORTED_BAM = os.path.join(RESULTS_DIR, OUTPUT_BAM_BASENAME)
SORTED_BAI = SORTED_BAM + ".bai"
ALIGNMENT_STATS_FILE = os.path.join(RESULTS_DIR, OUTPUT_STATS_BASENAME)

#%% Main function
# --- Helper function for running commands ---
def run_command(command, stdout_file=None, monitor_resources=False):
    """
    Runs a shell command.
    If monitor_resources=True, tracks CPU and memory usage during execution.
    Returns (success: bool, min_cpu, max_cpu, avg_cpu, min_mem, max_mem, avg_mem) if monitored,
    otherwise returns (success: bool).

    """
    if not monitor_resources:
        # Simple run without monitoring
        try:
            if stdout_file:
                with open(stdout_file, 'w') as f_out:
                    process = subprocess.run(command, stdout=f_out, stderr=subprocess.PIPE, text=True, check=True)
            else:
                process = subprocess.run(command, capture_output=True, text=True, check=True)

            if process.stderr:
                print(f"Stderr (if any):\n{process.stderr}")
            return True
        # Error handling
        except subprocess.CalledProcessError as e:
            print(f"Error: Command failed with exit code {e.returncode}")
            print(f"Command: {' '.join(e.cmd)}")
            if e.stdout:
                print(f"Stdout:\n{e.stdout}")
            if e.stderr:
                print(f"Stderr:\n{e.stderr}")
            return False
        except FileNotFoundError:
            print(f"Error: Command '{command[0]}' not found. Make sure it's in your PATH.")
            return False
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            return False

    else:
        # Run with resource monitoring
        max_cpu = 0.0
        max_mem = 0.0
        cpu_samples = []
        mem_samples = []
        try:
            process = subprocess.Popen(command, stdout=subprocess.PIPE if not stdout_file else open(stdout_file, 'w'), stderr=subprocess.PIPE)
            p = psutil.Process(process.pid)

            while True:
                if process.poll() is not None:
                    break

                cpu = p.cpu_percent(interval=1) # Check CPU usage
                mem = p.memory_info().rss / (1024 * 1024)  # Memory in MB

                cpu_samples.append(cpu)
                mem_samples.append(mem)

            # Final communicate to clean buffers
            stdout, stderr = process.communicate()

            # Calculate min, max, avg
            min_cpu = min(cpu_samples) if cpu_samples else 0.0
            max_cpu = max(cpu_samples) if cpu_samples else 0.0
            avg_cpu = sum(cpu_samples) / len(cpu_samples) if cpu_samples else 0.0

            min_mem = min(mem_samples) if mem_samples else 0.0
            max_mem = max(mem_samples) if mem_samples else 0.0
            avg_mem = sum(mem_samples) / len(mem_samples) if mem_samples else 0.0

            print(f"Max CPU usage: {max_cpu:.2f}%, Average CPU usage: {avg_cpu:.2f}%")
            print(f"Max Memory usage: {max_mem:.2f} MB, Average Memory usage: {avg_mem:.2f} MB")

            if stderr:
                print(f"Stderr:\n{stderr.decode()}")
            if process.returncode != 0:
                print(f"Error: Command failed with exit code {process.returncode}")
                return False, min_cpu, max_cpu, avg_cpu, min_mem, max_mem, avg_mem

            return True, min_cpu, max_cpu, avg_cpu, min_mem, max_mem, avg_mem

        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            return False, 0, max_cpu, 0, 0, max_mem, 0
 

print("\nStarting BWA-MEM Alignment Pipeline with predefined paths...")
print(f"Reference: {REFERENCE_GENOME_PATH}")
print(f"Read 1: {READ1_FASTQ_PATH}")
if IS_PAIRED_END:
    print(f"Read 2: {READ2_FASTQ_PATH}")
print(f"Output BAM: {SORTED_BAM}")
print(f"Output Stats: {ALIGNMENT_STATS_FILE}")
print(f"Threads: {NUM_THREADS}")
print(f"Intermediate SAM: {UNSORTED_SAM}")
print(f"Intermediate Unsorted BAM: {UNSORTED_BAM}")



#%%

# --- STEP 1: Aligning reads with bwa mem (saving to SAM) ---------------
print(f"\n--- STEP 1: Aligning reads with bwa mem (output SAM to {UNSORTED_SAM}) ---")

bwa_mem_cmd = ["bwa", "mem", "-t", str(NUM_THREADS), "-M", REFERENCE_GENOME_PATH, READ1_FASTQ_PATH]
if IS_PAIRED_END:
    bwa_mem_cmd.append(READ2_FASTQ_PATH)

# bwa mem with resource usage tracking
success, min_cpu, max_cpu, avg_cpu, min_mem, max_mem, avg_mem = run_command(bwa_mem_cmd, stdout_file=UNSORTED_SAM, monitor_resources=True)

if not success:
    print("Pipeline failed at BWA MEM step.")
    sys.exit(1)

max_cpu_limit = NUM_THREADS * 100 # cpu in percentage
if max_cpu_limit is not None and max_cpu > max_cpu_limit:
    print(f"WARNING: CPU usage exceeded limit: {max_cpu:.2f}% > {max_cpu_limit}%")

max_mem_limit = 16 * 1024  # Memory is limited to 16 GB anyway from line 150 set_memory_limit_gb(16)
if max_mem_limit is not None and max_mem > max_mem_limit:
    print(f"WARNING: Memory usage exceeded limit: {max_mem:.2f} MB > {max_mem_limit} MB")



# --- STEP 2: Converting SAM to Unsorted BAM ------------------
print(f"\n--- STEP 2: Converting SAM to Unsorted BAM  with samtools (output to {UNSORTED_BAM}) ---")
samtools_view_cmd = ["samtools", "view", "-b", UNSORTED_SAM] # -b means output BAM

if not run_command(samtools_view_cmd, stdout_file=UNSORTED_BAM):
    print("Pipeline failed at SAM to BAM conversion step.")
    sys.exit(1)
print(f"Unsorted BAM file created: {UNSORTED_BAM}")

# --- STEP 3: Sorting BAM file by coordinate ------------------
print(f"\n--- STEP 3: Sorting BAM file with samtools (output to {SORTED_BAM}) ---")
samtools_sort_cmd = ["samtools", "sort", UNSORTED_BAM, "-o", SORTED_BAM]

if not run_command(samtools_sort_cmd):
    print("Pipeline failed at BAM sorting step.")
    sys.exit(1)
print(f"Sorted BAM file created: {SORTED_BAM}")


# --- STEP 4: Indexing the Sorted BAM file ------------------
print(f"\n--- STEP 4: Indexing the Sorted BAM file with samtools (output to {SORTED_BAI}) ---")
samtools_index_cmd = ["samtools", "index", SORTED_BAM]

if not run_command(samtools_index_cmd):
    print("Pipeline failed at BAM indexing step.")
    sys.exit(1)
print(f"BAM index file created: {SORTED_BAI}")


# --- STEP 5: Collecting Alignment Statistics ------------------
print(f"\n--- STEP 5: Collecting Alignment Statistics with samtools (output to {ALIGNMENT_STATS_FILE}) ---")
samtools_flagstat_cmd = ["samtools", "flagstat", SORTED_BAM]

if not run_command(samtools_flagstat_cmd, stdout_file=ALIGNMENT_STATS_FILE):
    print("Pipeline failed at alignment statistics collection step.")
    sys.exit(1)
print(f"Alignment statistics saved to: {ALIGNMENT_STATS_FILE}")

print("\nPipeline finished successfully!")

# --- STEP 5: Prepare the human-readable report ------------------
report_txt = args.stats
parse_flagstat(ALIGNMENT_STATS_FILE, report_txt)

# Calculate average base quality (Phred) for read1 and possibly read2
avg_base_qual_1 = average_base_quality(READ1_FASTQ_PATH)
if READ2_FASTQ_PATH:
    avg_base_qual_2 = average_base_quality(READ2_FASTQ_PATH)
    avg_base_qual = (avg_base_qual_1 + avg_base_qual_2) / 2
else:
    avg_base_qual = avg_base_qual_1

# Calculate average MAPQ from BAM
avg_mapq = average_mapping_quality(SORTED_BAM)

# Read the existing alignment statistics
with open(report_txt, "r") as f:
    flagstat_report = f.read()

# Append the additional stats to the report
with open(report_txt, "a") as f:
    f.write("\n")
    f.write(f"Average base quality (Phred): {avg_base_qual:.2f}\n")
    f.write(f"Average mapping quality (MAPQ): {avg_mapq:.2f}\n")
    f.write(f"Max CPU usage (%): {max_cpu:.2f}\n")
    f.write(f"CPU usage (%) - min: {min_cpu:.2f}, max: {max_cpu:.2f}, avg: {avg_cpu:.2f}\n")
    f.write(f"Max Memory usage (MB): {max_mem:.2f}\n")
    f.write(f"Memory usage (MB) - min: {min_mem:.2f}, max: {max_mem:.2f}, avg: {avg_mem:.2f}\n")

print(f"Average base quality (Phred): {avg_base_qual:.2f}")
print(f"Average mapping quality (MAPQ): {avg_mapq:.2f}")



#%%
# --- Optional: Clean up intermediate files ---
if args.keep_intermediates == "no":
    # print("\nCleaning up intermediate files...")
    try:
        if os.path.exists(UNSORTED_SAM):
            os.remove(UNSORTED_SAM)
            print(f"Removed: {UNSORTED_SAM}")
        if os.path.exists(UNSORTED_BAM):
            os.remove(UNSORTED_BAM)
            print(f"Removed: {UNSORTED_BAM}")
        if os.path.exists(ALIGNMENT_STATS_FILE):
            os.remove(ALIGNMENT_STATS_FILE)
            print(f"Removed: {ALIGNMENT_STATS_FILE}")
    except OSError as e:
        print(f"Error during cleanup: {e}")

