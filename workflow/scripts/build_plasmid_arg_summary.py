import argparse
import csv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--amrfinder", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    candidates = {}

    with open(args.candidates, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            candidates[row["contig_id"]] = row

    gff = {}

    with open(args.gff) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue

            contig = fields[0]
            start = fields[3]
            end = fields[4]
            strand = fields[6]
            attributes = fields[8]

            protein_id = None

            for item in attributes.split(";"):
                if item.startswith("ID="):
                    internal_id = item.split("=", 1)[1]
                    protein_id = internal_id

            if protein_id:
                gff[protein_id] = {
                    "contig_id": contig,
                    "gene_start": start,
                    "gene_end": end,
                    "strand": strand,
                }

    amr_hits = []

    with open(args.amrfinder, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            protein_id = row["Protein id"]

            protein_gff = None

            # Match AMRFinder protein IDs to Prodigal CDS records.
            # Prodigal protein IDs are contig_geneNumber, while GFF
            # uses internal IDs such as 55065_6.
            parts = protein_id.rsplit("_", 1)

            if len(parts) == 2:
                contig_id, gene_number = parts
                target = f"{contig_id}_{gene_number}"

                # Convert the Prodigal gene number to the corresponding
                # CDS by counting CDS features on that contig.
                cds_count = 0
                for key, value in gff.items():
                    if value["contig_id"] == contig_id:
                        cds_count += 1
                        if cds_count == int(gene_number):
                            protein_gff = value
                            break

            if contig_id not in candidates:
                continue

            candidate = candidates[contig_id]

            amr_hits.append({
                "sample": candidate["sample"],
                "contig_id": contig_id,
                "length": candidate["length"],
                "gc_percent": candidate["gc_percent"],
                "coverage_depth": candidate["coverage_depth"],
                "coverage_breadth": candidate["coverage_breadth"],
                "protein_count": candidate["protein_count"],
                "replication_evidence": candidate["replication_evidence"],
                "mobility_evidence": candidate["mobility_evidence"],
                "maintenance_evidence": candidate["maintenance_evidence"],
                "conjugation_evidence": candidate["conjugation_evidence"],
                "mge_evidence": candidate["mge_evidence"],
                "evidence_status": candidate["evidence_status"],
                "protein_id": protein_id,
                "arg": row["Element symbol"],
                "arg_name": row["Element name"],
                "arg_class": row["Class"],
                "amrfinder_method": row["Method"],
                "arg_reference_coverage": row["% Coverage of reference"],
                "arg_identity": row["% Identity to reference"],
                "reference_accession": row["Closest reference accession"],
                "gene_start": protein_gff["gene_start"] if protein_gff else "",
                "gene_end": protein_gff["gene_end"] if protein_gff else "",
                "strand": protein_gff["strand"] if protein_gff else "",
            })

    fields = [
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
        "protein_id",
        "arg",
        "arg_name",
        "arg_class",
        "amrfinder_method",
        "arg_reference_coverage",
        "arg_identity",
        "reference_accession",
        "gene_start",
        "gene_end",
        "strand",
    ]

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(amr_hits)


if __name__ == "__main__":
    main()
