from mcc.helpers import check_options, define_time_requested, define_memory_requested

rule minimap2_to_viewpoints:
    input:
        fq="seqnado_output/flashed/{sample}/{sample}.extendedFrags.fastq.gz",
        viewpoints="seqnado_output/viewpoints.fa",
    output:
        bam=temp("seqnado_output/aligned/{sample}/{sample}.viewpoints.bam"),
    threads: 4
    resources:
        mem="4GB",
    container: None
    log:
        "seqnado_output/logs/aligned/{sample}.log",
    shell:
        """
        minimap2 -x sr -a -k 8 -w 1 --cs=long {input.viewpoints} {input.fq} 2> {log} |
        samtools view -h -b -o {output.bam} 2>> {log} &&
        samtools sort -o {output.bam}.sorted {output.bam} 2>> {log} &&
        mv {output.bam}.sorted {output.bam} &&
        samtools index {output.bam}
        """ 

rule align_to_genome:
    input:
        fq1="seqnado_output/split_reads/{sample}/{sample}.viewpoints.fastq.gz",
    params:
        index=config["genome"]["indices"],
        options=check_options(config["bowtie2"]["options"]),
    output:
        bam=temp("seqnado_output/aligned/raw/{sample}.bam"),
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    threads: config["bowtie2"]["threads"]
    log:
        "seqnado_output/logs/align/{sample}.log",
    shell:
        """bowtie2 -p {threads} -x {params.index} -U {input.fq1} {params.options} 2> {log} |
            samtools view -bS - > {output.bam} &&
            samtools sort -@ {threads} -o {output.bam}_sorted {output.bam} &&
            mv {output.bam}_sorted {output.bam} &&
            samtools index {output.bam}
        """

rule realign_unmapped:
    input:
        bam="seqnado_output/aligned/raw/{sample}.bam",
    output:
        bam=temp("seqnado_output/aligned/realigned/{sample}.bam"),
        bai=temp("seqnado_output/aligned/realigned/{sample}.bam.bai"),
    threads: config["samtools"]["threads"]
    resources:
        mem="4GB",
    log:
        "seqnado_output/logs/realign/{sample}.log",
    params:
        index=config["genome"]["indices"],
        options=check_options(config["bowtie2"]["options"]),

    shell:
        """
        samtools view -b -f 4 {input.bam} | bowtie2 -p {threads} -x {params.index} -b - --very-sensitive-local 2>> {log} |
        samtools view -bS - > {output.bam} &&
        samtools sort -@ {threads} -o {output.bam}_sorted {output.bam} &&
        mv {output.bam}_sorted {output.bam} &&
        samtools index {output.bam}
        """
    
rule combine_aligned:
    input:
        bam1="seqnado_output/aligned/raw/{sample}.bam",
        bam2="seqnado_output/aligned/realigned/{sample}.bam",
    output:
        bam=temp("seqnado_output/aligned/combined/{sample}.bam"),
    threads: config["samtools"]["threads"]
    resources:
        mem="2GB",
    log:
        "seqnado_output/logs/combine/{sample}.log",
    shell:
        """
        samtools merge -@ {threads} {output.bam} {input.bam1} {input.bam2} &&
        samtools view -F 4 -b {output.bam} > {output.bam}.tmp &&
        mv {output.bam}.tmp {output.bam} &&
        samtools sort -n -@ {threads} -o {output.bam}_sorted {output.bam} &&
        mv {output.bam}_sorted {output.bam}
        """