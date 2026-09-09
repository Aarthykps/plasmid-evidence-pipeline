#!/usr/bin/env python3

import argparse
import csv
import re
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a transparent plasmid evidence matrix from contig, coverage, Prodigal and DIAMOND results."
    )

    parser.add_argument("--features", required=True)
    parser.add_argument("--coverage", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--diamond", required=True)
    parser.add_argument("--annotations", required=True)

    parser.add_argument("--evidence-output", required=True)
    parser.add_argument("--candidate-output", required=True)
    parser.add_argument("--candidate-only-output", required=True)

    parser.add_argument("--min-identity", type=float, default=40.0)
    parser.add_argument("--min-query-coverage", type=float, default=50.0)
    parser.add_argument("--min-breadth", type=float, default=90.0)

    return parser.parse_args()


def contig_from_protein(protein_id):
    return re.sub(r"_[0-9]+$", "", protein_id)


def add_unique(mapping, value):
    if value not in mapping:
        mapping.append(value)


def load_features(path):
    data = {}

    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            data[row["contig_id"]] = row

    return data


def load_coverage(path):
    data = {}

    with open(path) as f:
        for line in f:
            if not line.strip():
                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9:
                continue

            contig = fields[0]

            data[contig] = {
                "breadth": float(fields[5]),
                "depth": float(fields[6]),
            }

    return data


def load_protein_counts(path):
    counts = defaultdict(int)

    with open(path) as f:
        for line in f:
            if not line.strip():
                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 1:
                continue

            contig = fields[0]
            counts[contig] += 1

    return counts


def load_annotations(path):
    annotations = {}

    with open(path, newline="") as f:
        reader = csv.reader(f, delimiter="\t")

        for row in reader:
            if len(row) < 2:
                continue

            accession = row[0]
            annotation = " ".join(row[1:]).strip()

            annotations[accession] = annotation

    return annotations


def classify_function(annotation):
    """
    Return a biologically interpretable evidence category.

    Important:
    - RepB/RepC-terminal plasmid partition proteins are maintenance evidence,
      not automatically replication evidence.
    - Relaxases are mobility evidence.
    - Generic 'replication' proteins are excluded unless plasmid-associated.
    """

    text = annotation.lower()

    # ---------------------------------------------------------
    # REPLICATION
    # ---------------------------------------------------------

    replication_patterns = [
        r"plasmid replication initiator",
        r"plasmid replication protein",
        r"plasmid replicon",
        r"plasmid.*replication.*protein",
        r"repabc",
        r"\btrfa\b",
        r"replication initiator trfa",
        r"replication initiation.*plasmid",
    ]

    for pattern in replication_patterns:
        if re.search(pattern, text):
            return "replication"

    # ---------------------------------------------------------
    # MOBILITY / RELAXASE
    # ---------------------------------------------------------

    mobility_patterns = [
        r"conjugative.*relaxase",
        r"relaxase/mobilization nuclease",
        r"mobilization nuclease",
        r"\bmobq\b.*relaxase",
        r"\bmobf\b.*relaxase",
        r"\bmobh\b.*relaxase",
        r"trai/moba.*relaxase",
        r"ti-type conjugative transfer relaxase",
        r"vird2",
    ]

    for pattern in mobility_patterns:
        if re.search(pattern, text):
            return "mobility"

    # ---------------------------------------------------------
    # MAINTENANCE / PARTITION
    # ---------------------------------------------------------

    maintenance_patterns = [
        r"plasmid partitioning protein",
        r"plasmid partition protein",
        r"plasmid.*partition.*repb",
        r"plasmid.*partition.*repc",
        r"plasmid stability",
        r"plasmid maintenance",
        r"parA.*plasmid",
        r"parB.*plasmid",
        r"toxin[- ]antitoxin",
        r"toxin.*antitoxin",
    ]

    for pattern in maintenance_patterns:
        if re.search(pattern, text):
            return "maintenance"

    # ---------------------------------------------------------
    # CONJUGATION / TRANSFER
    # ---------------------------------------------------------

    conjugation_patterns = [
        r"conjugative transfer protein",
        r"conjugation protein",
        r"type iv secretion",
        r"\btra[bcdefghijklmnpqrstuvwxyz]\b.*conjug",
        r"\btrb[bcdefghijklmnopqrs]\b",
    ]

    for pattern in conjugation_patterns:
        if re.search(pattern, text):
            return "conjugation"

    # ---------------------------------------------------------
    # MOBILE GENETIC ELEMENTS
    # ---------------------------------------------------------

    mge_patterns = [
        r"transposase",
        r"insertion sequence",
        r"\bintegrase\b",
        r"site-specific recombinase",
        r"mobile genetic element",
    ]

    for pattern in mge_patterns:
        if re.search(pattern, text):
            return "mge"

    return None


def main():

    args = parse_args()

    features = load_features(args.features)
    coverage = load_coverage(args.coverage)
    protein_counts = load_protein_counts(args.gff)
    annotations = load_annotations(args.annotations)

    evidence = defaultdict(lambda: {
        "replication": [],
        "mobility": [],
        "maintenance": [],
        "conjugation": [],
        "mge": [],
    })

    # ---------------------------------------------------------
    # Parse DIAMOND best hits
    # ---------------------------------------------------------

    with open(args.diamond) as f:

        for line in f:

            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 8:
                continue

            qseqid = fields[0]
            sseqid = fields[1]

            try:
                identity = float(fields[2])
                alignment_length = int(fields[3])
                qlen = int(fields[4])
                evalue = fields[6]
                bitscore = float(fields[7])
            except ValueError:
                continue

            qcov = (
                alignment_length / qlen * 100
                if qlen > 0
                else 0
            )

            if identity < args.min_identity:
                continue

            if qcov < args.min_query_coverage:
                continue

            annotation = annotations.get(sseqid)

            if annotation is None:
                continue

            category = classify_function(annotation)

            if category is None:
                continue

            contig = contig_from_protein(qseqid)

            evidence_string = (
                f"{qseqid}|{sseqid}|"
                f"{identity:.1f}%id|"
                f"{qcov:.1f}%qcov|"
                f"E={evalue}|"
                f"bitscore={bitscore:.0f}|"
                f"{annotation}"
            )

            add_unique(
                evidence[contig][category],
                evidence_string
            )

    # ---------------------------------------------------------
    # Write evidence matrix
    # ---------------------------------------------------------

    evidence_fields = [
        "sample",
        "contig_id",
        "length",
        "gc_percent",
        "coverage_depth",
        "coverage_breadth",
        "protein_count",
        "replication_evidence",
        "mobility_evidence",
        "maintenance_evidence",
        "conjugation_evidence",
        "mge_evidence",
        "evidence_status",
    ]

    candidate_fields = [
        "sample",
        "contig_id",
        "length",
        "gc_percent",
        "coverage_depth",
        "coverage_breadth",
        "protein_count",
        "replication_evidence",
        "mobility_evidence",
        "maintenance_evidence",
        "conjugation_evidence",
        "mge_evidence",
        "evidence_status",
    ]

    with open(args.evidence_output, "w", newline="") as out:

        writer = csv.DictWriter(
            out,
            fieldnames=evidence_fields,
            delimiter="\t"
        )

        writer.writeheader()

        for contig, feature in features.items():

            cov = coverage.get(
                contig,
                {"breadth": 0.0, "depth": 0.0}
            )

            rep = evidence[contig]["replication"]
            mob = evidence[contig]["mobility"]
            maint = evidence[contig]["maintenance"]
            conj = evidence[contig]["conjugation"]
            mge = evidence[contig]["mge"]

            breadth = cov["breadth"]
            depth = cov["depth"]

            good_read_support = breadth >= args.min_breadth

            # -------------------------------------------------
            # Transparent biological classification
            # -------------------------------------------------

            if good_read_support and rep:
                status = "Supported plasmid candidate"

            elif good_read_support and (
                mob or conj or (rep and maint)
            ):
                status = "Possible/unresolved plasmid"

            elif good_read_support and mge and not (
                rep or maint or mob or conj
            ):
                status = "Mobile-element-associated"

            else:
                status = "Unresolved / no strong plasmid evidence"

            row = {
                "sample": feature.get("sample", ""),
                "contig_id": contig,
                "length": feature.get("length", ""),
                "gc_percent": feature.get("gc_percent", ""),
                "coverage_depth": f"{depth:.4f}",
                "coverage_breadth": f"{breadth:.4f}",
                "protein_count": protein_counts.get(contig, 0),

                "replication_evidence":
                    " || ".join(rep),

                "mobility_evidence":
                    " || ".join(mob),

                "maintenance_evidence":
                    " || ".join(maint),

                "conjugation_evidence":
                    " || ".join(conj),

                "mge_evidence":
                    " || ".join(mge),

                "evidence_status":
                    status,
            }

            writer.writerow(row)

    # ---------------------------------------------------------
    # Write final candidate table
    # ---------------------------------------------------------

    candidate_statuses = {
        "Supported plasmid candidate",
        "Possible/unresolved plasmid",
    }

    with open(args.candidate_output, "w", newline="") as out, \
         open(args.candidate_only_output, "w", newline="") as out_only:

        writer = csv.DictWriter(
            out,
            fieldnames=candidate_fields,
            delimiter="\t"
        )

        writer_only = csv.DictWriter(
            out_only,
            fieldnames=candidate_fields,
            delimiter="\t"
        )

        writer.writeheader()
        writer_only.writeheader()

        with open(args.evidence_output, newline="") as evidence_file:

            reader = csv.DictReader(
                evidence_file,
                delimiter="\t"
            )

            for row in reader:

                if row["evidence_status"] not in candidate_statuses:
                    continue

                writer.writerow(row)
                writer_only.writerow(row)


if __name__ == "__main__":
    main()
