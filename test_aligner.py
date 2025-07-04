import os
import subprocess
import pytest
from pathlib import Path
import re
import time
import glob

aligning_script = "aligner.py"

reference_hg38  = "test_cases/reference_genome/hg38/hg38.fa.gz"
reference_chr21 = "test_cases/reference_genome/chr21/chr21.fa.gz"
reference_chr22 = "test_cases/reference_genome/chr22/chr22.fa"

case_1_R1 = "test_cases/test_case_1_R1.fastq.gz"
case_1_R2 = "test_cases/test_case_1_R2.fastq.gz"
case_2_R1 = "test_cases/test_case_2_R1.fastq.gz"
case_2_R2 = "test_cases/test_case_2_R2.fastq.gz"
case_3_R1 = "test_cases/test_case_3_R1.fastq.gz"
case_4_R1 = "test_cases/test_case_4_R1.fastq.gz"
case_4_R2 = "test_cases/test_case_4_R2.fastq.gz"

threads = 4

def merge_parts_if_missing(base_filename):
    """ This function is needed because I am limited to upload data to github. And i want the docker to be built as requested,
     being able to run end-to-end using hg38 full reference genome"""
    if os.path.exists(base_filename):
        print(f"{base_filename} already exists, skipping merge.")
        return
    
    pattern = base_filename + "_part_*"
    parts = sorted(glob.glob(pattern))
    if not parts:
        print(f"No parts found matching {pattern}. Cannot merge.")
        return
    
    print(f"Merging parts into {base_filename}...")
    with open(base_filename, 'wb') as wfd:
        for part in parts:
            with open(part, 'rb') as fd:
                while True:
                    chunk = fd.read(1024*1024)
                    if not chunk:
                        break
                    wfd.write(chunk)
    print("Merge completed.")

@pytest.fixture(scope="session", autouse=True)
def merge_reference():
    merge_parts_if_missing("test_cases/reference_genome/hg38/hg38.fa.gz.bwt")

@pytest.fixture
def setup_results_dir(tmp_path, request):
    test_name = getattr(request.function, "test_name", "default_test")
    dir_path = tmp_path / f"test_{test_name}"
    dir_path.mkdir()
    return dir_path

def named_test(name):
    def decorator(func):
        func.test_name = name
        return func
    return decorator

@named_test("req01")
def test_req01_paired_end_alignment(setup_results_dir):
    """
    REQ01: The aligner must support paired-end reads and produce a stats file with more than 0 total reads.
    Acceptance Criteria:
    - Aligner runs successfully (return code 0)
    - Bam file is created
    - Stats file is created
    - Stats file contains "Total reads" with value > 0
    """
    
    results_dir = str(setup_results_dir)
    test_name = "req01"
    report_name = "{}_report.txt".format(test_name)

    test_req01_paired_end_alignment.metadata = {
        "reads_R1": os.path.basename(case_1_R1),
        "reads_R2": os.path.basename(case_1_R2),
        "reference": os.path.basename(reference_hg38),
        "threads": threads
    }

    # Command to run aligner wrapper script
    cmd = [
        "python", aligning_script,
        "-r1", case_1_R1,
        "-r2", case_1_R2,
        "-f", reference_hg38,
        "-o", results_dir,
        "-t", str(threads),
        "--stats", report_name
    ]

    # Run the aligner wrapper
    result = subprocess.run(cmd, capture_output=True, text=True)
    # Check that it ran
    assert result.returncode == 0, f"Aligner failed:\n{result.stderr}"

    # Check the Bam files was created
    bam_files = [f for f in os.listdir(results_dir) if f.endswith(".bam")]
    assert bam_files, "No BAM file was generated for single-end input"

    # Check the stats file exists
    stats_file = os.path.join(results_dir, report_name)
    assert os.path.exists(stats_file), "Stats file missing"

    # Read stats and verify minimal criteria
    with open(stats_file) as f:
        stats = f.read()

    # Extract number of total reads need to be bigger than 0
    match = re.search(r"Total reads:\s*(\d+)", stats)
    assert match, "'Total reads' not found in stats"
    total_reads = int(match.group(1))
    assert total_reads > 0, f"Expected total reads > 0, got {total_reads}"



@named_test("req02")
def test_req02_single_end_alignment(setup_results_dir):
    """
    REQ02: The aligner should align reads in single-ended mode and produce a stats file with more than 0 total reads.
    Acceptance Criteria:
    - Aligner runs successfully (return code 0)
    - Bam file is created
    - Stats file is created
    - Stats file contains "Total reads" with value > 0
    """
    results_dir = str(setup_results_dir)
    test_name = "req02"
    report_name = "{}_report.txt".format(test_name)

    test_req02_single_end_alignment.metadata = {
        "reads_R1": os.path.basename(case_3_R1),
        "reads_R2": "--",
        "reference": os.path.basename(reference_hg38),
        "threads": threads
    }

    # Build command
    cmd = [
        "python", aligning_script,
        "-r1", case_3_R1,
        "-f", reference_hg38,
        "-o", results_dir,
        "-t", str(threads),
        "--stats", report_name
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"Aligner failed:\n{result.stderr}"

    bam_files = [f for f in os.listdir(results_dir) if f.endswith(".bam")]
    assert bam_files, "No BAM file was generated for single-end input"

    # Check the stats file exists
    stats_file = os.path.join(results_dir, report_name)
    assert os.path.exists(stats_file), "Stats file missing"

    # Read stats and verify minimal criteria
    with open(stats_file) as f:
        stats = f.read()

    # Extract number of total reads need to be bigger than 0
    match = re.search(r"Total reads:\s*(\d+)", stats)
    assert match, "'Total reads' not found in stats"
    total_reads = int(match.group(1))
    assert total_reads > 0, f"Expected total reads > 0, got {total_reads}"

@named_test("req03")
def test_req03_different_references(setup_results_dir):
    """
    REQ03: The aligner should work with different reference genomes
    Acceptance Criteria:
    - Aligner runs successfully for all different genomes (return code 0)
    - Bam files are created
    - Stats files are created
    - Stats files contains "Total reads" with value > 0 in all cases
    """
    results_dir = str(setup_results_dir)
    test_name = "req03"
    report_name = "{}_report.txt".format(test_name)
    
    
    references = [reference_chr21, reference_chr22, reference_hg38]
    names_references = f"{os.path.basename(reference_chr21)} + {os.path.basename(reference_chr22)} + {os.path.basename(reference_hg38)}"

    test_req03_different_references.metadata = {
        "reads_R1": os.path.basename(case_2_R1),
        "reads_R2": os.path.basename(case_2_R2),
        "reference": names_references,
        "threads": threads
    }

    print("\n ")
    for ref in references:
        print("RUN req03 with", ref)
        ref_name = Path(ref).stem
        output_dir = setup_results_dir / f"output_req03_{ref_name}"
        output_dir.mkdir(parents=True, exist_ok=True)

        report_name = "req3_report" + ref_name + ".txt"

        cmd = [
            "python", aligning_script,
            "-r1", case_2_R1,
            "-r2", case_2_R2,
            "-f", ref,
            "-o", output_dir,
            "-t", str(threads),
            "--stats", report_name
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Aligner failed with reference {ref}:\n{result.stderr}"

        # Confirm BAM exists
        bam_files = [f for f in os.listdir(output_dir) if f.endswith(".bam")]
        assert bam_files, f"No BAM file generated for reference {ref}"

        # Check the stats file exists
        stats_file = os.path.join(output_dir, report_name)
        assert os.path.exists(stats_file), f"No Stats file generated for reference {ref}"

        # Read stats and verify minimal criteria
        with open(stats_file) as f:
            stats = f.read()

        # Extract number of total reads need to be bigger than 0
        match = re.search(r"Total reads:\s*(\d+)", stats)
        assert match, "'Total reads' not found in stats"
        total_reads = int(match.group(1))
        assert total_reads > 0, f"Expected total reads > 0, got {total_reads}, for reference {ref}"


@named_test("req04")
def test_req04_alignment_stats_present(setup_results_dir):
    """REQ04: The aligner wrapper should correctly report all required alignment statistics
    Acceptance Criteria:
    - Aligner runs successfully for all different genomes (return code 0)
    - Bam files are created
    - Stats files are created
    - Stats files contains the information parsed in expected patterns
    """

    results_dir = str(setup_results_dir)
    test_name = "req04"
    report_name = "{}_report.txt".format(test_name)


    test_req04_alignment_stats_present.metadata = {
        "reads_R1": os.path.basename(case_2_R1),
        "reads_R2": os.path.basename(case_2_R2),
        "reference": os.path.basename(reference_hg38),
        "threads": threads
    }

    cmd = [
        "python", aligning_script,
        "-r1", case_2_R1,
        "-r2", case_2_R2,
        "-f", reference_hg38,
        "-o", results_dir,
        "-t", str(threads),
        "--stats", report_name
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"Aligner failed: {result.stderr}"

    stats_file = os.path.join(results_dir, report_name)
    assert os.path.exists(stats_file), "Stats file missing"


    with open(stats_file, "r") as f:
        stats_content = f.read()

    # Define expected patterns to find in the report  file
    expected_patterns = [
        r"Total reads\s*:\s*\d+",
        r"Mapped reads\s*:\s*\d+ \(\d+(\.\d+)?%\)",
        r"Unmapped reads\s*:\s*\d+ \(\d+(\.\d+)?%\)",
        r"Duplicated reads\s*:\s*\d+ \(\d+(\.\d+)?%\)",
        r"Singletons\s*:\s*\d+ \(\d+(\.\d+)?%\)",
        r"Average base quality \(Phred\):\s*\d+(\.\d+)?",
        r"Average mapping quality \(MAPQ\):\s*\d+(\.\d+)?",
    ]

    for pattern in expected_patterns:
        assert re.search(pattern, stats_content), f"Missing or wrong field: {pattern}"


@named_test("req05")
def test_req05_resource_limits(setup_results_dir):
    """REQ05: The aligner wrapper should not exceed the pre-defined computational resources (CPU & memory) allocated to it, for 1M reads using hg38
    Acceptance Criteria:
    - Aligner runs successfully for all different genomes (return code 0)
    - Bam files are created
    - The process should report less than 400% CPU usage and 16Mb memory
    """
    # Define your resource limits for this test
    MAX_CPU_PERCENT = 400
    MAX_MEMORY_MB = 16 * 1024

    results_dir = str(setup_results_dir)
    test_name = "req05"
    report_name = "{}_report.txt".format(test_name)

    test_req05_resource_limits.metadata = {
        "reads_R1": os.path.basename(case_4_R1),
        "reads_R2": os.path.basename(case_4_R2),
        "reference": os.path.basename(reference_hg38),
        "threads": threads
    }

    cmd = [
        "python3", aligning_script,
        "-r1", case_4_R1,
        "-r2", case_4_R2,
        "-f", reference_hg38,
        "-o", results_dir,
        "-t", str(threads),
        "--stats", report_name
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"Aligner failed: {result.stderr}"

    stats_file = os.path.join(results_dir, report_name)
    assert os.path.exists(stats_file), "Stats file missing"

    with open(stats_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith("Max CPU usage (%):"):
                    max_cpu = float(line.split(":")[1].strip())
                    print(max_cpu)
                elif line.startswith("Max Memory usage (MB):"):
                    max_mem = float(line.split(":")[1].strip())

    assert max_cpu <= MAX_CPU_PERCENT, f"CPU usage exceeded limit: {max_cpu} > {MAX_CPU_PERCENT}"
    assert max_mem <= MAX_MEMORY_MB, f"Memory usage exceeded limit: {max_mem} > {MAX_MEMORY_MB}"




@named_test("req06")
def test_req06_runtime_small_input(setup_results_dir):
    """REQ06: The aligner wrapper should run fast when provided with small input data
    Acceptance Criteria:
    - Aligner runs successfully (return code 0)
    - Bam files are created
    - The process should run in less than 120 seconds for small imputs (~10000 mapped reads)
    """

    results_dir = str(setup_results_dir)
    test_name = "req06"
    report_name = "{}_report.txt".format(test_name)

    test_req06_runtime_small_input.metadata = {
        "reads_R1": os.path.basename(case_1_R1),
        "reads_R2": os.path.basename(case_1_R2),
        "reference": os.path.basename(reference_hg38),
        "threads": threads
    }

    # Build the command with all needed flags
    cmd = [
        "python", aligning_script,
        "-r1", case_1_R1,
        "-r2", case_1_R2,
        "-f", reference_hg38,
        "-t", str(threads),
        "-o", results_dir,
        "--stats", report_name
    ]
    
    start_time = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    end_time = time.time()
    
    runtime = end_time - start_time
    print(f"Runtime for small input: {runtime:.2f} seconds")
    
    assert result.returncode == 0, "Aligner failed to run successfully"
    assert runtime < 120, f"Runtime exceeded 2 minutes: {runtime:.2f} seconds"

@named_test("extra_cores")
def test_extra_different_cores(setup_results_dir):
    """
    Extra test: The aligner should work with different setting for threads (cores)
    Acceptance Criteria:
    - Aligner runs successfully using 2, 3 and 4 threads. (can be extended to more if needed)
    - Bam files are created
    - generate a report for all the cases
    """
    results_dir = str(setup_results_dir)
    test_name = "extra_cores"
    report_name = "{}_report.txt".format(test_name)

    threads_list = [2, 3, 4]
    threads_list_names = " + ".join(str(x) for x in threads_list)

    test_extra_different_cores.metadata = {
        "reads_R1": os.path.basename(case_1_R1),
        "reads_R2": os.path.basename(case_1_R2),
        "reference": os.path.basename(reference_hg38),
        "threads": threads_list_names
    }

    print("\n ")
    for thread in threads_list:
        print("RUN Using threads = ", thread)
        ref_name = str(thread)
        output_dir = setup_results_dir / f"output_extra_core_{ref_name}"
        output_dir.mkdir(parents=True, exist_ok=True)

        report_name = "Extra_core_report" + ref_name + ".txt"

        cmd = [
            "python", aligning_script,
            "-r1", case_2_R1,
            "-r2", case_2_R2,
            "-f", reference_hg38,
            "-o", output_dir,
            "-t", str(thread),
            "--stats", report_name
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Aligner failed with threads = {thread}:\n{result.stderr}"

        # Confirm BAM exists
        bam_files = [f for f in os.listdir(output_dir) if f.endswith(".bam")]
        assert bam_files, f"No BAM file generated for threads = {thread}"

        # Confirm report exists
        report_path = output_dir / report_name
        assert report_path.exists(), f"No report generated for threads = {thread}"

