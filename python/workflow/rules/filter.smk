from mcc.helpers import check_options

rule filter_bam:
    input:
        bam="seqnado_output/aligned/raw/{sample}.bam",
        bai="seqnado_output/aligned/raw/{sample}.bam.bai",
    output:
        bam="seqnado_output/aligned/filtered/{sample}.bam",
        bai="seqnado_output/aligned/filtered/{sample}.bam.bai",
    threads: config["samtools"]["threads"]
    resources:
        mem="500MB",
    log:
        "seqnado_output/logs/filter/{sample}.log",
    params:
        options=check_options(config["samtools"]["filter_options"]),
    shell:
        """
        samtools view -@ {threads} -h -b {input.bam} {params.options} > {output.bam} &&
        samtools index {output.bam} &&
        echo 'Filtered reads' > {log} 2>&1
        """