configfile: "config/config.yaml"

import csv

with open(config["samples"], newline="") as handle:
    samples = {
        row["sample"]: row
        for row in csv.DictReader(handle, delimiter="\t")
    }

SAMPLES = list(samples.keys())


rule plasmid_arg_summary:
    input:
        candidates="results/final/{sample}.plasmid_candidates_only.tsv",
        amrfinder="results/arg/{sample}.amrfinder.tsv",
        gff="results/genes/{sample}.gff"
    output:
        summary="results/final/{sample}.plasmid_arg_summary.tsv"
    log:
        "logs/final/{sample}.plasmid_arg_summary.log"
    shell:
        """
        mkdir -p results/final logs/final

        python workflow/scripts/build_plasmid_arg_summary.py \
            --candidates {input.candidates} \
            --amrfinder {input.amrfinder} \
            --gff {input.gff} \
            --output {output.summary} \
            > {log} 2>&1

        test -s {output.summary}
        """

rule all:
    input:
        expand(
            "results/qc/raw/{sample}_R1_fastqc.html",
            sample=SAMPLES
        ),
        expand(
            "results/qc/raw/{sample}_R1_fastqc.zip",
            sample=SAMPLES
        ),
        expand(
            "results/qc/raw/{sample}_R2_fastqc.html",
            sample=SAMPLES
        ),
        expand(
            "results/qc/raw/{sample}_R2_fastqc.zip",
            sample=SAMPLES
        ),
        "results/qc/multiqc/raw/multiqc_report.html",
        expand(
            "results/trimmed/{sample}_R1.trimmed.fastq.gz",
            sample=SAMPLES
        ),
        expand(
            "results/trimmed/{sample}_R2.trimmed.fastq.gz",
            sample=SAMPLES
        ),
        expand(
            "results/qc/fastp/{sample}.json",
            sample=SAMPLES
        ),
        expand(
            "results/qc/fastp/{sample}.html",
            sample=SAMPLES
        ),
        expand(
            "results/qc/quality_gate/{sample}.quality_gate.txt",
            sample=SAMPLES
        ),
        expand(
            "results/assembly/{sample}/final.contigs.fa",
            sample=SAMPLES
        ),
        expand(
            "results/features/{sample}.contig_features.tsv",
            sample=SAMPLES
        ),
        expand(
            "results/coverage/{sample}.coverage.tsv",
            sample=SAMPLES
        ),
        expand(
              "results/genes/{sample}.fna",
              sample=SAMPLES
        ),
        expand(
              "results/genes/{sample}.faa",
              sample=SAMPLES
        ),
        expand(
              "results/genes/{sample}.gff",
              sample=SAMPLES
        ),
         expand(
            "results/annotation/{sample}.diamond.tsv",
            sample=SAMPLES
        ),
        expand(
            "results/annotation/{sample}.diamond.best.tsv",
            sample=SAMPLES
        ),
        expand(
            "results/evidence/{sample}.plasmid_evidence.tsv",
            sample=SAMPLES
        ),
        expand(
            "results/final/{sample}.plasmid_candidates.tsv",
            sample=SAMPLES
        ),
        expand(
            "results/final/{sample}.plasmid_candidates_only.tsv",
            sample=SAMPLES
        ),
        expand(
            "results/arg/{sample}.amrfinder.tsv",
            sample=SAMPLES
        ),

        expand(
            "results/final/{sample}.plasmid_arg_summary.tsv",
            sample=SAMPLES
        )
rule fastqc_raw:
    input:
        r1=lambda wildcards: samples[wildcards.sample]["fq1"],
        r2=lambda wildcards: samples[wildcards.sample]["fq2"]
    output:
        html1="results/qc/raw/{sample}_R1_fastqc.html",
        zip1="results/qc/raw/{sample}_R1_fastqc.zip",
        html2="results/qc/raw/{sample}_R2_fastqc.html",
        zip2="results/qc/raw/{sample}_R2_fastqc.zip"
    threads: config["threads"]["fastqc"]
    log:
        "logs/fastqc_raw/{sample}.log"
    conda:
        "workflow/envs/qc.yaml"
    shell:
        """
        mkdir -p results/qc/raw logs/fastqc_raw

        fastqc \
            --threads {threads} \
            --outdir results/qc/raw \
            {input.r1} {input.r2} \
            > {log} 2>&1
        """


rule multiqc_raw:
    input:
        expand(
            "results/qc/raw/{sample}_R1_fastqc.zip",
            sample=SAMPLES
        ),
        expand(
            "results/qc/raw/{sample}_R2_fastqc.zip",
            sample=SAMPLES
        )
    output:
        "results/qc/multiqc/raw/multiqc_report.html"
    threads: config["threads"]["multiqc"]
    log:
        "logs/multiqc_raw.log"
    conda:
        "workflow/envs/qc.yaml"
    shell:
        """
        mkdir -p results/qc/multiqc/raw logs

        multiqc \
            results/qc/raw \
            --outdir results/qc/multiqc/raw \
            --filename multiqc_report.html \
            --force \
            > {log} 2>&1
        """


rule fastp:
    input:
        r1=lambda wildcards: samples[wildcards.sample]["fq1"],
        r2=lambda wildcards: samples[wildcards.sample]["fq2"]
    output:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz",
        json="results/qc/fastp/{sample}.json",
        html="results/qc/fastp/{sample}.html"
    threads: config["threads"]["fastp"]
    params:
        quality=config["fastp"]["qualified_quality_phred"],
        min_length=config["fastp"]["length_required"]
    log:
        "logs/fastp/{sample}.log"
    conda:
        "workflow/envs/fastp.yaml"
    shell:
        """
        mkdir -p results/trimmed results/qc/fastp logs/fastp

        fastp \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1} \
            --out2 {output.r2} \
            --thread {threads} \
            --detect_adapter_for_pe \
            --qualified_quality_phred {params.quality} \
            --length_required {params.min_length} \
            --json {output.json} \
            --html {output.html} \
            > {log} 2>&1

        test -s {output.r1}
        test -s {output.r2}
        test -s {output.json}
        test -s {output.html}
        """


rule fastqc_trimmed:
    input:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz"
    output:
        html1="results/qc/trimmed/{sample}_R1.trimmed_fastqc.html",
        zip1="results/qc/trimmed/{sample}_R1.trimmed_fastqc.zip",
        html2="results/qc/trimmed/{sample}_R2.trimmed_fastqc.html",
        zip2="results/qc/trimmed/{sample}_R2.trimmed_fastqc.zip"
    threads: config["threads"]["fastqc"]
    log:
        "logs/fastqc_trimmed/{sample}.log"
    conda:
        "workflow/envs/qc.yaml"
    shell:
        """
        mkdir -p results/qc/trimmed logs/fastqc_trimmed

        fastqc \
            --threads {threads} \
            --outdir results/qc/trimmed \
            {input.r1} {input.r2} \
            > {log} 2>&1
        """


rule multiqc_trimmed:
    input:
        expand(
            "results/qc/trimmed/{sample}_R1.trimmed_fastqc.zip",
            sample=SAMPLES
        ),
        expand(
            "results/qc/trimmed/{sample}_R2.trimmed_fastqc.zip",
            sample=SAMPLES
        )
    output:
        "results/qc/multiqc/trimmed/multiqc_report.html"
    threads: config["threads"]["multiqc"]
    log:
        "logs/multiqc_trimmed.log"
    conda:
        "workflow/envs/qc.yaml"
    shell:
        """
        mkdir -p results/qc/multiqc/trimmed logs

        multiqc \
            results/qc/trimmed \
            --outdir results/qc/multiqc/trimmed \
            --filename multiqc_report.html \
            --force \
            > {log} 2>&1
        """

rule quality_gate:
    input:
        json="results/qc/fastp/{sample}.json"
    output:
        "results/qc/quality_gate/{sample}.quality_gate.txt"
    params:
        minimum_reads=config["quality_gate"]["minimum_reads"],
        minimum_q30_rate=config["quality_gate"]["minimum_q30_rate"]
    log:
        "logs/quality_gate/{sample}.log"
    shell:
        """
        mkdir -p results/qc/quality_gate logs/quality_gate

        python workflow/scripts/quality_gate.py \
            --json {input.json} \
            --output {output} \
            --minimum-reads {params.minimum_reads} \
            --minimum-q30-rate {params.minimum_q30_rate} \
            > {log} 2>&1

        test -s {output}
        """

rule megahit:
    input:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz",
        quality_gate="results/qc/quality_gate/{sample}.quality_gate.txt"
    output:
        contigs="results/assembly/{sample}/final.contigs.fa"
    threads: config["threads"]["megahit"]
    log:
        "logs/megahit/{sample}.log"
    conda:
        "workflow/envs/megahit.yaml"
    shell:
        """
        mkdir -p logs/megahit
        rm -rf results/assembly/{wildcards.sample}

        megahit \
            -1 {input.r1} \
            -2 {input.r2} \
            --out-dir {output.contigs.dirname} \
            --out-prefix final \
            --num-cpu-threads {threads} \
            > {log} 2>&1

        test -s {output.contigs}
        """

rule contig_stats:
    input:
        contigs="results/assembly/{sample}/final.contigs.fa"
    output:
        stats="results/assembly/{sample}/contig_stats.tsv"
    log:
        "logs/contig_stats/{sample}.log"
    shell:
        """
        mkdir -p logs/contig_stats

        python workflow/scripts/contig_stats.py \
            --input {input.contigs} \
            --output {output.stats} \
            > {log} 2>&1

        test -s {output.stats}
        """

rule filter_contigs:
    input:
        contigs="results/assembly/{sample}/final.contigs.fa"
    output:
        filtered500="results/assembly/{sample}/filtered_contigs_500bp.fa",
        filtered1000="results/assembly/{sample}/filtered_contigs_1000bp.fa"
    log:
        "logs/filter_contigs/{sample}.log"
    shell:
        """
        mkdir -p logs/filter_contigs

        python workflow/scripts/filter_contigs.py \
            --input {input.contigs} \
            --output {output.filtered500} \
            --min-length 500 \
            > {log} 2>&1

        python workflow/scripts/filter_contigs.py \
            --input {input.contigs} \
            --output {output.filtered1000} \
            --min-length 1000 \
            >> {log} 2>&1

        test -s {output.filtered500}
        test -s {output.filtered1000}
        """

rule contig_features:
    input:
        # Stage 2 uses the >=1000 bp filtered contigs from Stage 1
        contigs="results/assembly/{sample}/filtered_contigs_1000bp.fa"

    output:
        # One feature table containing one row per contig
        features="results/features/{sample}.contig_features.tsv"

    log:
        # Keep a separate log so the analysis is reproducible/debuggable
        "logs/contig_features/{sample}.log"

    shell:
        """
        # Create directories needed for Stage 2
        mkdir -p results/features logs/contig_features

        # Run the version-controlled feature extraction script
        python workflow/scripts/contig_features.py \
            --input {input.contigs} \
            --output {output.features} \
            --sample {wildcards.sample} \
            > {log} 2>&1

        # Make sure the expected output was actually created
        test -s {output.features}
        """
rule bowtie2_index:
    input:
        contigs="results/assembly/{sample}/filtered_contigs_1000bp.fa"
    output:
        idx=multiext(
            "results/coverage/index/{sample}.filtered_contigs",
            ".1.bt2", ".2.bt2", ".3.bt2",
            ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"
        )
    log:
        "logs/coverage/index/{sample}.log"
    conda:
        "workflow/envs/mapping.yaml"
    shell:
        """
        mkdir -p results/coverage/index logs/coverage/index

        bowtie2-build \
            {input.contigs} \
            results/coverage/index/{wildcards.sample}.filtered_contigs \
            > {log} 2>&1

        for f in {output.idx}; do
            test -s "$f"
        done
        """


rule map_reads:
    input:
        index=rules.bowtie2_index.output.idx,
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz"
    output:
        bam="results/coverage/{sample}.sorted.bam"
    threads: config["threads"]["megahit"]
    log:
        "logs/coverage/map/{sample}.log"
    conda:
        "workflow/envs/mapping.yaml"
    shell:
        """
        mkdir -p results/coverage logs/coverage/map

        bowtie2 \
            -x results/coverage/index/{wildcards.sample}.filtered_contigs \
            -1 {input.r1} \
            -2 {input.r2} \
            --threads {threads} \
            2> {log} \
        | samtools sort \
            -@ {threads} \
            -o {output.bam} -

        samtools quickcheck -v {output.bam}
        test -s {output.bam}
        """


rule index_bam:
    input:
        bam="results/coverage/{sample}.sorted.bam"
    output:
        bai="results/coverage/{sample}.sorted.bam.bai"
    log:
        "logs/coverage/index_bam/{sample}.log"
    conda:
        "workflow/envs/mapping.yaml"
    shell:
        """
        mkdir -p logs/coverage/index_bam

        samtools index \
            {input.bam} \
            {output.bai} \
            > {log} 2>&1

        test -s {output.bai}
        """


rule contig_coverage:
    input:
        bam="results/coverage/{sample}.sorted.bam",
        bai="results/coverage/{sample}.sorted.bam.bai"
    output:
        coverage="results/coverage/{sample}.coverage.tsv"
    log:
        "logs/coverage/coverage/{sample}.log"
    conda:
        "workflow/envs/mapping.yaml"
    shell:
        """
        mkdir -p logs/coverage/coverage

        samtools coverage \
            {input.bam} \
            > {output.coverage} \
            2> {log}

        test -s {output.coverage}
        """
rule prodigal:
    input:
        fasta="results/assembly/{sample}/filtered_contigs_1000bp.fa"
    output:
        cds="results/genes/{sample}.fna",
        proteins="results/genes/{sample}.faa",
        gff="results/genes/{sample}.gff"
    log:
        "logs/prodigal/{sample}.log"
    conda:
        "workflow/envs/prodigal.yaml"
    shell:
        """
        mkdir -p results/genes logs/prodigal

        prodigal \
            -i {input.fasta} \
            -d {output.cds} \
            -a {output.proteins} \
            -o {output.gff} \
            -f gff \
            -p meta \
            > {log} 2>&1
        """


rule diamond_plasmid_search:
    input:
        proteins="results/genes/{sample}.faa",
        db="resources/diamond/refseq_plasmid_proteins.dmnd"
    output:
        hits="results/annotation/{sample}.diamond.tsv"
    log:
        "logs/diamond/{sample}.log"
    conda:
        "workflow/envs/diamond.yaml"
    threads: 8
    shell:
        """
        mkdir -p results/annotation logs/diamond

        diamond blastp \
            --query {input.proteins} \
            --db {input.db} \
            --out {output.hits} \
            --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \
            --evalue 1e-5 \
            --query-cover 50 \
            --subject-cover 50 \
            --max-target-seqs 10 \
            --threads {threads} \
            > {log} 2>&1
        """
rule diamond_best_hits:
    input:
        hits="results/annotation/{sample}.diamond.tsv"
    output:
        best="results/annotation/{sample}.diamond.best.tsv"
    log:
        "logs/diamond/{sample}.best.log"
    shell:
        """
        mkdir -p results/annotation logs/diamond

        awk 'BEGIN{{FS="\\t"; OFS="\\t"}}
        {{
            qcov=$4/$5*100;
            scov=$4/$6*100;
            if ($3>=40 && qcov>=50 && scov>=50)
                print $0, qcov, scov
        }}' {input.hits} |
        sort -k1,1 -k8,8nr |
        awk '!seen[$1]++' > {output.best}

        test -s {output.best}
        """


rule plasmid_evidence:
    input:
        features="results/features/{sample}.contig_features.tsv",
        coverage="results/coverage/{sample}.coverage.tsv",
        gff="results/genes/{sample}.gff",
        diamond="results/annotation/{sample}.diamond.best.tsv",
        annotations="results/annotation/{sample}.hit_annotations.tsv"
    output:
        evidence="results/evidence/{sample}.plasmid_evidence.tsv",
        candidates="results/final/{sample}.plasmid_candidates.tsv",
        candidates_only="results/final/{sample}.plasmid_candidates_only.tsv"
    log:
        "logs/evidence/{sample}.log"
    shell:
        """
        mkdir -p results/evidence results/final logs/evidence

        python workflow/scripts/plasmid_evidence.py \
            --features {input.features} \
            --coverage {input.coverage} \
            --gff {input.gff} \
            --diamond {input.diamond} \
            --annotations {input.annotations} \
            --evidence-output {output.evidence} \
            --candidate-output {output.candidates} \
            --candidate-only-output {output.candidates_only} \
            --min-identity 40 \
            --min-query-coverage 50 \
            --min-breadth 90 \
            > {log} 2>&1

        test -s {output.evidence}
        test -s {output.candidates}
        test -s {output.candidates_only}
        """
rule extract_candidate_proteins:
    input:
        candidates="results/final/{sample}.plasmid_candidates_only.tsv",
        proteins="results/genes/{sample}.faa"
    output:
        proteins="results/arg/{sample}.plasmid_candidate_proteins.faa"
    log:
        "logs/arg/{sample}.extract_proteins.log"
    conda:
        "workflow/envs/amrfinder.yaml"
    shell:
        """
        mkdir -p results/arg logs/arg

        python workflow/scripts/extract_candidate_proteins.py \
            --candidates {input.candidates} \
            --proteins {input.proteins} \
            --output {output.proteins} \
            > {log} 2>&1

        test -s {output.proteins}
        """
rule amrfinder:
    input:
        proteins="results/arg/{sample}.plasmid_candidate_proteins.faa"
    output:
        hits="results/arg/{sample}.amrfinder.tsv"
    log:
        "logs/arg/{sample}.amrfinder.log"
    conda:
        "workflow/envs/amrfinder.yaml"
    threads: 4
    shell:
        """
        mkdir -p results/arg logs/arg

        amrfinder -u >> {log} 2>&1

        amrfinder \
            -p {input.proteins} \
            -o {output.hits} \
            >> {log} 2>&1

        test -s {output.hits}
        """
