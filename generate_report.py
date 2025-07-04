import csv
import argparse
from datetime import datetime
import os

def generate_markdown(input_csv, output_md):
    with open(input_csv, newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = list(reader)

    date_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    report_lines = []
    report_lines.append(f"# Date {date_str}\n")
    report_lines.append("# Test results:\n")

    for row in rows:
        req = row["requirement"]
        criteria = row["acceptance_criteria"]
        reads = row.get("reads_R1", "")
        if row.get("reads_R2"):
            reads += f" + {row['reads_R2']}"
        reference = row.get("reference", "")
        threads = row.get("threads", "")
        result = row["result"]

        report_lines.append(f"## Test {req}")
        report_lines.append(f"Acceptance criteria: {criteria}")
        report_lines.append(f"Reads used: {reads}")
        report_lines.append(f"Reference genome: {reference}")
        report_lines.append(f"Threads: {threads}")
        report_lines.append(f"Result: {result}\n")

    with open(output_md, "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))

    print(f"Markdown report generated: {output_md}")
    # # Delete the input CSV file after report generation
    # try:
    #     os.remove(input_csv)
    #     print(f"Deleted input CSV file: {input_csv}")
    # except Exception as e:
    #     print(f"Could not delete {input_csv}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Generate Markdown test report from CSV results")
    parser.add_argument("-i", "--input", required=True, help="Input CSV file (e.g., test_results.csv)")
    parser.add_argument("-o", "--output", default="test_report.md", help="Output Markdown file (e.g., test_report.md)")
    args = parser.parse_args()

    generate_markdown(args.input, args.output)

if __name__ == "__main__":
    main()
