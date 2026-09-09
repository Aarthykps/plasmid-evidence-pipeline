#!/usr/bin/env python3

import argparse
import csv
import re

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--sample", required=True)
args = parser.parse_args()


def write_record(header, sequence, writer):
    contig_id = header.split()[0]

    flag_match = re.search(r"\bflag=([^\s]+)", header)
    multi_match = re.search(r"\bmulti=([^\s]+)", header)

    length = len(sequence)
    a_count = sequence.count("A")
    c_count = sequence.count("C")
    g_count = sequence.count("G")
    t_count = sequence.count("T")
    n_count = sequence.count("N")

    gc_percent = (
        100.0 * (g_count + c_count) / length
        if length else 0.0
    )

    writer.writerow([
        args.sample,
        contig_id,
        length,
        f"{gc_percent:.4f}",
        a_count,
        c_count,
        g_count,
        t_count,
        n_count,
        flag_match.group(1) if flag_match else "",
        multi_match.group(1) if multi_match else "",
    ])


with open(args.input) as fasta, open(args.output, "w", newline="") as output:
    writer = csv.writer(output, delimiter="\t")

    writer.writerow([
        "sample",
        "contig_id",
        "length",
        "gc_percent",
        "a_count",
        "c_count",
        "g_count",
        "t_count",
        "n_count",
        "megahit_flag",
        "megahit_multi",
    ])

    header = None
    sequence_parts = []

    for line in fasta:
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            if header is not None:
                write_record(header, "".join(sequence_parts), writer)

            header = line[1:]
            sequence_parts = []
        else:
            sequence_parts.append(line.upper())

    if header is not None:
        write_record(header, "".join(sequence_parts), writer)
