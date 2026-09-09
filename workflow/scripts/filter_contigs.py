#!/usr/bin/env python3

import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--min-length", type=int, required=True)
args = parser.parse_args()


def write_record(header, sequence, output):
    if len(sequence) >= args.min_length:
        output.write(header + "\n")

        for start in range(0, len(sequence), 80):
            output.write(sequence[start:start + 80] + "\n")


with open(args.input) as fasta, open(args.output, "w") as output:
    header = None
    sequence_parts = []

    for line in fasta:
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            if header is not None:
                write_record(
                    header,
                    "".join(sequence_parts),
                    output
                )

            header = line
            sequence_parts = []
        else:
            sequence_parts.append(line)

    if header is not None:
        write_record(
            header,
            "".join(sequence_parts),
            output
        )
