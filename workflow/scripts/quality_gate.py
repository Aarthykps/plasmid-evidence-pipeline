#!/usr/bin/env python3

import argparse
import json
import sys


parser = argparse.ArgumentParser()
parser.add_argument("--json", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--minimum-reads", type=int, required=True)
parser.add_argument("--minimum-q30-rate", type=float, required=True)
args = parser.parse_args()

with open(args.json) as handle:
    report = json.load(handle)

after = report["summary"]["after_filtering"]

total_reads = int(after["total_reads"])
q30_rate = float(after["q30_rate"])

if "read1_before_filtering" in report["summary"]:
    before_r1 = report["summary"]["read1_before_filtering"]["total_reads"]
    before_r2 = report["summary"]["read2_before_filtering"]["total_reads"]
    after_r1 = report["summary"]["read1_after_filtering"]["total_reads"]
    after_r2 = report["summary"]["read2_after_filtering"]["total_reads"]
    paired_counts_match = before_r1 == before_r2 and after_r1 == after_r2
else:
    paired_counts_match = True

checks = {
    "minimum_reads": total_reads >= args.minimum_reads,
    "minimum_q30_rate": q30_rate >= args.minimum_q30_rate,
    "paired_counts_match": paired_counts_match,
}

status = "PASS" if all(checks.values()) else "FAIL"

with open(args.output, "w") as handle:
    handle.write(f"status\t{status}\n")
    handle.write(f"total_reads\t{total_reads}\n")
    handle.write(f"q30_rate\t{q30_rate:.6f}\n")
    handle.write(f"minimum_reads\t{args.minimum_reads}\n")
    handle.write(f"minimum_q30_rate\t{args.minimum_q30_rate:.6f}\n")
    handle.write(f"paired_counts_match\t{paired_counts_match}\n")

if status == "FAIL":
    sys.exit("Quality gate failed")
