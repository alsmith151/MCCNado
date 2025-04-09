from mcc.helpers import define_time_requested


rule fastqc_raw:
    input:
        fq1="seqnado_output/fastqs/{sample}_1.fastq.gz",
        fq2="seqnado_output/fastqs/{sample}_2.fastq.gz",
    output:
        html1="seqnado_output/qc/fastqc_raw/{sample}_1_fastqc.html",
        html2="seqnado_output/qc/fastqc_raw/{sample}_2_fastqc.html",
        zip1="seqnado_output/qc/fastqc_raw/{sample}_1_fastqc.zip",
        zip2="seqnado_output/qc/fastqc_raw/{sample}_2_fastqc.zip",
    params:
        extra="--quiet",
        output_dir="seqnado_output/qc/fastqc_raw/",
        temp_prefix="seqnado_output/qc/fastqc_raw/{sample}",
    threads: 1
    resources:
        mem="1.5GB",
         runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        "seqnado_output/logs/fastqc_raw/{sample}.log",
    shell:
        """
        fastqc -o {params.output_dir} {input.fq1} {input.fq2} > {log} 2>&1
        """

use rule fastqc_raw as fastqc_trimmed with:
    input:
        fq1="seqnado_output/trimmed/{sample}_1.fastq.gz",
        fq2="seqnado_output/trimmed/{sample}_2.fastq.gz",
    output:
        html1="seqnado_output/qc/fastqc_trimmed/{sample}_1_fastqc.html",
        html2="seqnado_output/qc/fastqc_trimmed/{sample}_2_fastqc.html",
        zip1="seqnado_output/qc/fastqc_trimmed/{sample}_1_fastqc.zip",
        zip2="seqnado_output/qc/fastqc_trimmed/{sample}_2_fastqc.zip",
    params:
        extra="--quiet",
        output_dir="seqnado_output/qc/fastqc_trimmed/",
        temp_prefix="seqnado_output/qc/fastqc_trimmed/{sample}",
    log:
        "seqnado_output/logs/fastqc_trimmed/{sample}.log",
    
