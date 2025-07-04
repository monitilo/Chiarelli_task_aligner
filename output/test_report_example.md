# Date 2025-07-04 08:47:02

# Test results:

## Test REQ01
Acceptance criteria: REQ01: The aligner must support paired-end reads and produce a stats file with more than 0 total reads.
Reads used: test_case_1_R1.fastq.gz + test_case_1_R2.fastq.gz
Reference genome: hg38.fa.gz
Threads: 4
Result: PASS

## Test REQ02
Acceptance criteria: REQ02: The aligner should align reads in single-ended mode and produce a stats file with more than 0 total reads.
Reads used: test_case_3_R1.fastq.gz + --
Reference genome: hg38.fa.gz
Threads: 4
Result: PASS

## Test REQ03
Acceptance criteria: REQ03: The aligner should work with different reference genomes
Reads used: test_case_2_R1.fastq.gz + test_case_2_R2.fastq.gz
Reference genome: chr21.fa.gz + chr22.fa + hg38.fa.gz
Threads: 4
Result: PASS

## Test REQ04
Acceptance criteria: REQ04: The aligner wrapper should correctly report all required alignment statistics
Reads used: test_case_2_R1.fastq.gz + test_case_2_R2.fastq.gz
Reference genome: hg38.fa.gz
Threads: 4
Result: PASS

## Test REQ05
Acceptance criteria: REQ05: The aligner wrapper should not exceed the pre-defined computational resources (CPU & memory) allocated to it, for 1M reads using hg38
Reads used: test_case_4_R1.fastq.gz + test_case_4_R2.fastq.gz
Reference genome: hg38.fa.gz
Threads: 4
Result: FAIL

## Test REQ06
Acceptance criteria: REQ06: The aligner wrapper should run fast when provided with small input data
Reads used: test_case_1_R1.fastq.gz + test_case_1_R2.fastq.gz
Reference genome: hg38.fa.gz
Threads: 4
Result: PASS

## Test EXTRA_CORES
Acceptance criteria: Extra test: The aligner should work with different setting for threads (cores)
Reads used: test_case_1_R1.fastq.gz + test_case_1_R2.fastq.gz
Reference genome: hg38.fa.gz
Threads: 2 + 3 + 4
Result: PASS
