import csv
import subprocess
import os
from pathlib import Path


test_results = []

def pytest_runtest_makereport(item, call):
    # This hook is called when a test phase finishes
    # I want only the 'call' phase (test execution)
    if call.when == "call":
        outcome = "PASS" if call.excinfo is None else "FAIL"

        # Extract custom metadata from test function attributes
        test_name = getattr(item.function, "test_name", item.name)

        # get the docstring as acceptance criteria (if you want)
        acceptance_criteria = (item.function.__doc__ or "").strip().split('\n')[0]

        metadata = getattr(item.function, "metadata", {})

        test_results.append({
            "requirement": test_name.upper(),
            "acceptance_criteria": acceptance_criteria,
            "test_case": item.name,
            "result": outcome,
            "reads_R1": metadata.get("reads_R1", ""),
            "reads_R2": metadata.get("reads_R2", ""),
            "reference": metadata.get("reference", ""),
            "threads": metadata.get("threads", "")
        })

def pytest_sessionfinish(session, exitstatus):
    
    cwd = Path(os.getcwd())
    test_result_dir = cwd / "output"
    test_result_dir.mkdir(parents=True, exist_ok=True)

    # Write CSV file
    csv_file = os.path.join(test_result_dir, "test_results.csv")

    with open(csv_file, "w", newline='', encoding="utf-8") as f:
        fieldnames = ["requirement", "acceptance_criteria", "test_case", "result",
                    "reads_R1", "reads_R2", "reference", "threads"]
        writer = csv.DictWriter(f, fieldnames)
        writer.writeheader()
        for row in test_results:
            writer.writerow(row)
    print(" CSV file written: test_results.csv")

    # Auto-generate Markdown report calling generate_report.py
    report_file = os.path.join(test_result_dir, "test_report.md")
    try:
        subprocess.run(["python", "generate_report.py", "-i", csv_file, "-o", report_file], check=True)
        print(" Markdown report generated: test_report.md")
    except Exception as e:
        print(" Could not generate Markdown report:", e)

