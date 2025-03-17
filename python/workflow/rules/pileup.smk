from typing import List, Dict

rule pileup:
    input:
        bam="seqnado_output/aligned/split_genomic_reads/{sample}/{viewpoint}.bam",
        bai="seqnado_output/aligned/split_genomic_reads/{sample}/{viewpoint}.bam.bai",
        exclusion_regions="seqnado_output/exclusion_regions.bed",
    output:
        pileup="seqnado_output/bigwigs/deeptools/unscaled/{sample}/{viewpoint}.bigWig",
    log:
        "seqnado_output/logs/pileup/{sample}/{viewpoint}.log",
    threads: 8
    resources:
        mem="4GB",
    params:
        options="",
    shell:
        """
        bamCoverage -p {threads} -b {input.bam} -o {output.pileup} {params.options} --blackListFileName {input.exclusion_regions} 2> {log}
        """


def define_pileup_files_with_shared_viewpoints(samples) -> Dict[str, Dict[str, str]]:
    import pathlib

    pileups = dict()
    for sample in samples:
        checkpoint_output = checkpoints.split_genomic_reads.get(sample=sample)
        outdir = pathlib.Path(checkpoint_output.output.bams)
        viewpoints = glob_wildcards(str(outdir / "{viewpoint}.bam")).viewpoint
    
        pileups[sample] = {
            'viewpoints': set(viewpoints),
            'files': {viewpoint: f"seqnado_output/bigwigs/deeptools/unscaled/{sample}/{viewpoint}.bigWig" for viewpoint in viewpoints}
        }
    
    return pileups


def define_inputs_for_pileup(wildcards, samples) -> List[str]:
    pileups = define_pileup_files_with_shared_viewpoints(samples)
    files = []
    for sample in samples:
        if wildcards.viewpoint in pileups[sample]['viewpoints']:
            files.append(pileups[sample]['files'][wildcards.viewpoint])
    return files

rule combine_pileups:
    input:
        bws=lambda wc: define_inputs_for_pileup(wc, samples),
        sentinels=expand("seqnado_output/pileup/sentinel/{samples}.sentinel", samples=samples),
    output:
        bw="seqnado_output/bigwigs/deeptools/combined/combine_{viewpoint}.bigWig",
    log:
        "seqnado_output/logs/combine_pileups_{viewpoint}.log",
    threads: 8
    resources:
        mem="3GB",
    shell:
        """
        bigwigAverage -b {input.bws} -o {output.bw} -p {threads} 2> {log}
        """