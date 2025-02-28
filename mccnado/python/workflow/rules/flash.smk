rule flash:
    input:
        fq1="seqnado_output/trimmed/{sample}_1.fastq.gz",
        fq2="seqnado_output/trimmed/{sample}_2.fastq.gz",
    output:
        flashed="seqnado_output/flashed/{sample}/{sample}.extendedFrags.fastq.gz",
        pe1="seqnado_output/flashed/{sample}/{sample}.notCombined_1.fastq.gz",
        pe2="seqnado_output/flashed/{sample}/{sample}.notCombined_2.fastq.gz",
        hist=temp(
            "seqnado_output/flashed/{sample}/{sample}.hist"
        ),
        histogram=temp(
            "seqnado_output/flashed/{sample}/{sample}.histogram"
        ),
    params:
        outdir="seqnado_output/flashed/{sample}/{sample}",
    threads: 4
    resources:
        mem_mb=1000,
    container: None #alternative == "docker://pegi3s/flash"
    log:
        "seqnado_output/logs/flash/{sample}.log",
    shell:
        """
        flash {input.fq1} {input.fq2} -o {params.outdir} -t {threads} -z --compress-prog-args pigz > {log} 2>&1
        """