import os
import glob

# Get the list of transcriptomes dynamically
transcriptomes = [d for d in os.listdir("data/transcriptomes") if os.path.isdir(os.path.join("data/transcriptomes", d))]

rule all:
    input:
        "output/filtered_insertions.csv",
        "output/filtered_insertions.fasta",
        expand("blast_db/{transcriptome}.nhr", transcriptome=transcriptomes),
        expand("blast_db/{transcriptome}.nin", transcriptome=transcriptomes),
        expand("blast_db/{transcriptome}.nsq", transcriptome=transcriptomes),
        expand("output/blast_results/{transcriptome}.tbl", transcriptome=transcriptomes),
        "output/merged_blast_results.csv"  # Add this line

# Original insertion processing rule
rule process_insertions:
    input:
        "data/insertions.json"
    output:
        csv="output/filtered_insertions.csv",
        fasta="output/filtered_insertions.fasta"
    params:
        script="scripts/process_insertions.py"
    shell:
        "python {params.script} --input {input} --output_csv {output.csv} --output_fasta {output.fasta}"

# Concatenation of .fa.gz files for each transcriptome
rule concat_fasta:
    input:
        lambda wildcards: glob.glob("data/transcriptomes/" + wildcards.transcriptome + "/*.fa.gz")
    output:
        concat="data/transcriptomes/{transcriptome}/concatenated.fasta.gz"
    shell:
        "cat {input} > {output.concat}"

rule build_blast_db:
    input:
        "data/transcriptomes/{transcriptome}/concatenated.fasta.gz"
    output:
        nhr="blast_db/{transcriptome}.nhr",
        nin="blast_db/{transcriptome}.nin",
        nsq="blast_db/{transcriptome}.nsq"
    params:
        db_prefix="blast_db/{transcriptome}"
    log:
        "logs/blast_db/{transcriptome}.log"
    shell:
        """
        gunzip -c {input} > {output.nhr}.temp.fasta
        makeblastdb -in {output.nhr}.temp.fasta -dbtype nucl -out {params.db_prefix} 2> {log}
        rm {output.nhr}.temp.fasta
        """

rule run_blast:
    input:
        fasta="output/filtered_insertions.fasta",
        db_nhr="blast_db/{transcriptome}.nhr",
        db_nin="blast_db/{transcriptome}.nin",
        db_nsq="blast_db/{transcriptome}.nsq"
    output:
        tbl="output/blast_results/{transcriptome}.tbl"
    params:
        db_prefix="blast_db/{transcriptome}"
    shell:
        """
        blastn -evalue 10 -word_size 7  -query {input.fasta} -db {params.db_prefix} -outfmt 7 -num_threads 4 > {output.tbl}.temp
        grep -v '^#' {output.tbl}.temp > {output.tbl}
        rm {output.tbl}.temp
        """

rule merge_blast_results:
    input:
        csv="output/filtered_insertions.csv",
        tbl=expand("output/blast_results/{transcriptome}.tbl", transcriptome=transcriptomes)
    output:
        merged_csv="output/merged_blast_results.csv"
    params:
        script="scripts/merge_blast_results.py"
    shell:
        "python {params.script} --input_csv {input.csv} --input_tbl_dir output/blast_results --output_csv {output.merged_csv}"
