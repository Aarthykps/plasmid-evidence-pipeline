import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--proteins", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    candidates = set()

    with open(args.candidates) as f:
        next(f)
        for line in f:
            if not line.strip():
                continue
            candidates.add(line.rstrip("\n").split("\t")[1])

    keep = False

    with open(args.proteins) as fin, open(args.output, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                protein_id = line.split()[0][1:]
                contig_id = protein_id.rsplit("_", 1)[0]
                keep = contig_id in candidates

            if keep:
                fout.write(line)


if __name__ == "__main__":
    main()
