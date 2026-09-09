#!/usr/bin/env python3

import argparse


def read_lengths(path):
    lengths = []
    current = 0

    with open(path) as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
            else:
                current += len(line)

    if current:
        lengths.append(current)

    return lengths


def n50(lengths):
    total = sum(lengths)
    running = 0

    for length in sorted(lengths, reverse=True):
        running += length
        if running >= total / 2:
            return length

    return 0


parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

lengths = read_lengths(args.input)

if not lengths:
    raise SystemExit("No contigs found")

thresholds = [200, 500, 1000, 1500, 2000, 5000, 10000]

with open(args.output, "w") as out:
    out.write("metric\tvalue\n")
    out.write(f"total_contigs\t{len(lengths)}\n")
    out.write(f"total_bases\t{sum(lengths)}\n")
    out.write(f"minimum_length\t{min(lengths)}\n")
    out.write(f"maximum_length\t{max(lengths)}\n")
    out.write(f"average_length\t{sum(lengths) / len(lengths):.2f}\n")
    out.write(f"N50\t{n50(lengths)}\n")

    for threshold in thresholds:
        selected = [x for x in lengths if x >= threshold]
        out.write(
            f"contigs_ge_{threshold}bp\t{len(selected)}\n"
        )
        out.write(
            f"bases_ge_{threshold}bp\t{sum(selected)}\n"
        )
