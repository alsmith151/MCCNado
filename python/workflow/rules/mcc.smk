rule split_viewpoint_reads:
    input:
        bam="seqnado_output/aligned/{sample}/{sample}.viewpoints.bam",
    output:
        fq="seqnado_output/split_reads/{sample}/{sample}.viewpoints.fastq.gz",
    threads: 1
    resources:
        mem="1GB",
    log:
        "seqnado_output/logs/split_reads/{sample}.log",
    container: None
    run:
        # import logging
        # FORMAT = '%(levelname)s %(name)s %(asctime)-15s %(filename)s:%(lineno)d %(message)s'
        # logging.basicConfig(format=FORMAT)
        # logging.getLogger().setLevel(logging.INFO)
        from mcc import mcc        
        
        mcc.split_viewpoint_reads(input.bam, output.fq)

checkpoint split_genomic_reads:
    input:
        bam="seqnado_output/aligned/combined/{sample}.bam",
    output:
        bams=directory("seqnado_output/aligned/split_genomic_reads/{sample}"),
    params:
        output_dir="seqnado_output/aligned/split_genomic_reads/{sample}/",
    threads: 1
    resources:
        mem="1GB",
    log:
        "seqnado_output/logs/split_genomic_reads/{sample}.log",
    container: None
    run:
        from mcc import mcc
        import pathlib

        outdir = pathlib.Path(params.output_dir)
        outdir.mkdir(exist_ok=True, parents=True)
        mcc.split_genomic_reads(input.bam, params.output_dir)


rule index_bam:
    input:
        bam="seqnado_output/aligned/split_genomic_reads/{sample}/{viewpoint}.bam",
    output:
        bai="seqnado_output/aligned/split_genomic_reads/{sample}/{viewpoint}.bam.bai",
    threads: 1
    resources:
        mem="1GB",
    log:
        "seqnado_output/logs/index_bam/{sample}/{viewpoint}.log",
    shell:
        """
        samtools sort -@ {threads} -o {input.bam}.sorted {input.bam} &&
        mv {input.bam}.sorted {input.bam} &&
        samtools index {input.bam} {output.bai}
        """


def identify_extracted_bam_files(wildcards):
    import pathlib

    checkpoint_output = checkpoints.split_genomic_reads.get(**wildcards)
    outdir = pathlib.Path(checkpoint_output.output.bams)
    viewpoints = glob_wildcards(str(outdir / "{viewpoint}.bam")).viewpoint
    return expand(str(outdir / "{viewpoint}.bam"), viewpoint=viewpoints)

def define_pileup_files(wildcards):
    import pathlib

    checkpoint_output = checkpoints.split_genomic_reads.get(**wildcards)
    outdir = pathlib.Path(checkpoint_output.output.bams)
    viewpoints = glob_wildcards(str(outdir / "{viewpoint}.bam")).viewpoint
    return expand("seqnado_output/bigwigs/deeptools/unscaled/{sample}/{viewpoint}.bigWig", viewpoint=viewpoints, sample=wildcards.sample)  

def redefine_viewpoints(samples):
    viewpoint_set = set()
    
    for ii, sample in enumerate(samples):
        checkpoint_output = checkpoints.split_genomic_reads.get(sample=sample)
        outdir = pathlib.Path(checkpoint_output.output.bams)
        viewpoints = glob_wildcards(str(outdir / "{viewpoint}.bam")).viewpoint
        
        if ii == 0:
            viewpoint_set = set(viewpoints)
        else:
            viewpoint_set = viewpoint_set.intersection(viewpoints)
    return list(viewpoint_set)


rule pileup_sentinel:
    input:
        bams=define_pileup_files,
    output:
        sentinel="seqnado_output/pileup/sentinel/{sample}.sentinel",
    params:
        output_dir="seqnado_output/pileup/{sample}/",
    threads: 1
    resources:
        mem="1GB",
    log:
        "seqnado_output/logs/pileup/{sample}.log",
    container: None
    shell:
        """
        touch {output.sentinel}
        """


ruleorder:
    split_genomic_reads > index_bam