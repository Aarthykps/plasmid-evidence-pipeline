# Metagenomic Plasmid-Evidence Pipeline

A reproducible Snakemake workflow for identifying and characterizing candidate plasmid-containing contigs from metagenomic sequencing data using integrated sequence, read-mapping, gene-function, and mobile-element evidence.

## Overview

Plasmids are important vehicles for horizontal gene transfer and can contribute to the dissemination of antimicrobial resistance (AMR) genes. However, identifying plasmids directly from metagenomic assemblies is challenging because sequence characteristics alone do not establish plasmid origin.

This project uses an evidence-integration strategy rather than relying on a single plasmid classifier.

The pipeline integrates:

- Read quality control and trimming
- Metagenomic assembly with MEGAHIT
- Contig sequence-feature analysis
- Read mapping and coverage estimation
- Metagenomic gene prediction with Prodigal
- DIAMOND similarity searches against RefSeq plasmid-associated proteins
- Functional evidence for replication, mobilization, conjugation, maintenance, and mobile genetic elements
- Candidate plasmid extraction
- Antimicrobial resistance gene detection using AMRFinderPlus
- Integrated evidence summaries

> No single feature is treated as definitive evidence of plasmid origin. Candidate classification is based on converging biological and sequencing evidence.

## Workflow

```text
Raw paired-end reads
        |
        v
   FastQC / MultiQC
        |
        v
      fastp
        |
        v
FastQC / MultiQC (trimmed)
        |
        v
    Quality gate
        |
        v
      MEGAHIT
        |
        v
Contig filtering (>=1 kb)
        |
        +-----------------------------+
        |                             |
        v                             v
 Contig features                 Read mapping
 length / GC%                    Bowtie2 + BAM
        |                             |
        +-------------+---------------+
                      |
                      v
                Prodigal -meta
                      |
                      v
               Predicted proteins
                      |
                      v
          DIAMOND against RefSeq
            plasmid proteins
                      |
                      v
              Best-hit filtering
                      |
                      v
              Evidence integration
                      |
          +-----------+-----------+
          |                       |
          v                       v
 Candidate plasmids          Unresolved / MGE
          |
          v
 Candidate protein extraction
          |
          v
      AMRFinderPlus
          |
          v
 ARG-containing plasmid candidates
          |
          v
    Final integrated summary
```

## Scientific rationale

### 1. Read quality control

FastQC and MultiQC are used to inspect raw and trimmed sequencing quality. fastp performs paired-end read preprocessing and quality filtering.

A quality gate is applied before assembly to reduce the risk of poor-quality input propagating into downstream analyses.

### 2. Metagenomic assembly

MEGAHIT is used to assemble the processed reads into contigs.

The main downstream dataset is filtered to contigs of at least 1 kb. This reduces extremely short sequences that provide limited biological context for plasmid-associated gene interpretation.

### 3. Contig sequence features

For each filtered contig, the pipeline records contig length, GC percentage, base composition, and MEGAHIT-derived contig information.

These features provide descriptive sequence context but are not treated as independent proof of plasmid origin.

### 4. Read-mapping support

Processed reads are mapped back to the assembly using Bowtie2.

For each contig, the pipeline calculates mean read depth, coverage breadth, and covered bases.

Read support is used primarily to establish that an assembled contig is supported by the sequencing data. High coverage alone does not imply plasmid origin.

### 5. Gene prediction

Prodigal is run in metagenomic mode to predict coding sequences and translated proteins.

The predicted proteins provide the functional units required for downstream plasmid-associated protein annotation.

### 6. Plasmid-associated protein similarity

Predicted proteins are searched against a RefSeq plasmid protein database using DIAMOND BLASTP.

Qualifying hits are filtered using minimum identity, query coverage, subject coverage, and E-value criteria. The highest-bitscore qualifying hit is retained for each query protein.

Similarity to a protein present in a plasmid reference database is treated as supporting evidence rather than automatic proof that the source contig is a plasmid.

### 7. Evidence integration

The pipeline integrates functional evidence associated with replication, mobilization, conjugation, maintenance, and mobile genetic elements.

Candidate classification uses transparent biological categories rather than an arbitrary numerical score.

### 8. AMR analysis

AMRFinderPlus is applied only to proteins extracted from plasmid candidates.

This ordering is deliberate: plasmid evidence is evaluated first, followed by candidate selection and then AMR gene detection.

An antimicrobial resistance gene is therefore not used by itself to classify a contig as a plasmid.

## Dataset and final results

The workflow was evaluated on one paired-end metagenomic sample.

After assembly and filtering:

- 59,481 contigs >=1 kb were retained
- 201,318 proteins were predicted by Prodigal
- 79,124 proteins had qualifying best plasmid-protein hits after filtering
- 107 contigs were retained as supported or possible/unresolved plasmid candidates
- 1 contig was classified as a supported plasmid candidate
- 1 AMR gene was detected among candidate-plasmid proteins

Final evidence distribution:

| Evidence status | Contigs |
|---|---:|
| Supported plasmid candidate | 1 |
| Possible/unresolved plasmid | 106 |
| Mobile-element-associated | 1,176 |
| Unresolved / no strong plasmid evidence | 58,198 |
| **Total** | **59,481** |

## Strongest plasmid candidate

### Contig: k141_401600

The strongest candidate identified by the integrated evidence framework was k141_401600.

| Feature | Result |
|---|---:|
| Contig length | 115,638 bp |
| GC content | 55.6270% |
| Mean read depth | 22.9537x |
| Coverage breadth | 99.9914% |
| Predicted proteins | 106 |
| Evidence status | Supported plasmid candidate |

The contig contains multiple functionally distinct plasmid-associated proteins, providing convergent evidence rather than relying on a single annotation.

### Functional evidence on k141_401600

| Protein | Functional annotation | Identity | Query coverage | E-value |
|---|---|---:|---:|---:|
| k141_401600_86 | TrfA-like plasmid replication initiator | 57.8% | 86.2% | 1.79e-106 |
| k141_401600_26 | MobH family relaxase | 50.0% | 111.6% | 1.86e-256 |
| k141_401600_23 | TrbC conjugation-associated protein | 59.4% | 76.0% | 3.37e-74 |
| k141_401600_63 | TraC type IV secretion system protein | 83.3% | 99.3% | 0 |
| k141_401600_70 | TraA family conjugative transfer protein | 73.7% | 99.0% | 1.39e-42 |
| k141_401600_2 | Tn3 family transposase | 100.0% | 99.9% | 0 |

The combination of a plasmid replication-associated protein, a relaxase, and several conjugation-associated proteins provides stronger plasmid-associated functional evidence than any individual similarity hit alone. The Tn3-family transposase additionally indicates mobile genetic element activity on the same contig.

The contig also contains a blaOXA antimicrobial resistance gene identified independently by AMRFinderPlus. Its position downstream of the Tn3-family transposase is consistent with a mobile-element-associated genomic context, although proximity alone does not establish direct mobilization of the resistance gene.

## Antimicrobial resistance evidence

AMRFinderPlus was applied to proteins extracted from the 107 plasmid candidates. One antimicrobial resistance hit was identified.

| Feature | Result |
|---|---|
| Contig | k141_401600 |
| Protein | k141_401600_6 |
| AMR gene | blaOXA |
| Gene class | Class D beta-lactamase |
| Detection method | HMM |
| Reference coverage | 100.00% |
| Identity | 74.62% |
| Closest reference | WP_032488481.1 |
| HMM accession | NF012161.0 |

The blaOXA gene is located at coordinates 6594-7376 on k141_401600. The predicted protein is 261 amino acids long and covers 100% of the AMRFinderPlus reference sequence.

## Interpretation

k141_401600 was classified as a supported plasmid candidate because multiple independent evidence types converge on the same contig. Read mapping shows 99.99% coverage breadth with a mean depth of 22.95x, demonstrating strong sequencing support for the assembled sequence.

More importantly, the contig contains a TrfA-like replication initiator together with a MobH relaxase and several conjugation-associated proteins. These annotations provide functional evidence consistent with plasmid replication and horizontal transfer capability.

The presence of a Tn3-family transposase and blaOXA further indicates that the contig contains mobile-element and antimicrobial-resistance-associated features. However, the analysis does not establish complete plasmid circularity or experimental transfer, so the sequence is reported as a supported plasmid candidate rather than a confirmed plasmid.

## Conclusion

This project demonstrates an evidence-integration approach for plasmid identification from metagenomic assemblies. Rather than treating a single classifier, coverage value, or AMR annotation as definitive, the workflow combines assembly quality, read support, protein-level similarity, plasmid-associated functions, and targeted AMR annotation to identify candidates with interpretable biological evidence.
## Software and methods

| Analysis stage | Tool / script | Role in the pipeline |
|---|---|---|
| Workflow management | Snakemake | Defines dependencies, executes analysis rules, and supports reproducible reruns |
| Raw-read quality control | FastQC | Evaluates quality characteristics of the raw paired-end reads |
| QC report aggregation | MultiQC | Combines FastQC and other QC outputs into summary reports |
| Read preprocessing | fastp | Performs paired-end read preprocessing, adapter detection, and quality filtering |
| Post-trimming quality control | FastQC + MultiQC | Verifies read quality after preprocessing |
| Quality gate | `quality_gate.py` | Checks whether the processed reads meet the configured minimum read-count and Q30 criteria before assembly |
| Metagenomic assembly | MEGAHIT | Assembles processed reads into metagenomic contigs |
| Contig statistics | `contig_stats.py` | Summarizes assembly contig characteristics |
| Contig filtering | `filter_contigs.py` | Retains contigs meeting the configured minimum length threshold |
| Contig feature analysis | `contig_features.py` | Calculates sequence-level features including length and GC content |
| Read mapping | Bowtie2 | Maps processed reads back to the assembled contigs |
| BAM processing | SAMtools | Sorts, checks, indexes, and calculates coverage statistics from mapped reads |
| Gene prediction | Prodigal | Predicts coding sequences and translated proteins in metagenomic mode |
| Protein similarity search | DIAMOND BLASTP | Searches predicted proteins against the RefSeq plasmid protein database |
| Best-hit selection | `diamond_best_hits` rule | Applies identity and coverage thresholds and retains the highest-bitscore qualifying hit for each protein |
| Plasmid evidence integration | `plasmid_evidence.py` | Integrates contig, coverage, and plasmid-associated functional evidence into candidate classifications |
| Candidate protein extraction | `extract_candidate_proteins.py` | Extracts predicted proteins belonging to selected plasmid candidates |
| AMR detection | AMRFinderPlus | Detects antimicrobial resistance genes among proteins from candidate plasmids |
| Final AMR summary | `build_plasmid_arg_summary.py` | Produces the final integrated plasmid and AMR summary |
