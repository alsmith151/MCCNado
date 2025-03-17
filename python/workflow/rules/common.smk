rule viewpoints_to_fasta:
    input:
        bed=config["viewpoints"],
        genome=config["genome"]["fasta"],
    output:
        fasta="seqnado_output/viewpoints.fa",
    log:
        "seqnado_output/logs/bed_to_fasta/viewpoints.log",
    
    shell:
        """
        bedtools getfasta -fi {input.genome} -bed {input.bed} -fo {output.fasta} -name 2> {log} &&
        cat {output.fasta} | sed -E 's/:+/-/g' > {output.fasta}.tmp &&
        mv {output.fasta}.tmp {output.fasta}
        """


rule fasta_index:
    input:
        fasta="seqnado_output/viewpoints.fa",
    output:
        index="seqnado_output/viewpoints.fa.fai",
    log:
        "seqnado_output/logs/bed_to_fasta/index.log",
    shell:
        """
        samtools faidx {input.fasta} -o {output.index}
        """

rule exclusion_regions:
    input:
        bed=config['viewpoints'],
    output:
        bed="seqnado_output/exclusion_regions.bed"
    log:
        "seqnado_output/logs/exclusion_regions.log"
    params:
        genome=config['genome']['chromosome_sizes'],
        exclusion_zone=config.get("exclusion_zone", 500)
    shell:
        """
        bedtools slop -i {input.bed} -g {params.genome}  -b {params.exclusion_zone} > {output.bed}
        """
