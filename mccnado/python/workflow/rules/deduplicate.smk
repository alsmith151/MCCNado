
rule unzip:
    input:
        "seqnado_output/fastqs/{sample}_1.fastq.gz",
        "seqnado_output/fastqs/{sample}_2.fastq.gz",
    output:
        fq1=temp("seqnado_output/fastqs/{sample}_1.fastq"),
        fq2=temp("seqnado_output/fastqs/{sample}_2.fastq"),
    threads: 1
    resources:
        mem="1GB",
    shell:
        """
        gunzip -c {input[0]} > {output.fq1}
        gunzip -c {input[1]} > {output.fq2}
        """



rule deduplicate_fastq_raw:
    input:
        fq1="seqnado_output/fastqs/{sample}_1.fastq",
        fq2="seqnado_output/fastqs/{sample}_2.fastq",
    output:
        deduped1=temp("seqnado_output/deduped/{sample}/{sample}_1.fastq.gz"),
        deduped2=temp("seqnado_output/deduped/{sample}/{sample}_2.fastq.gz"),
    threads: 1
    resources:
        mem="1GB",
    log:
        "seqnado_output/logs/deduplication/{sample}.log",
    run:
        from mcc import mcc
        import json

        stats = mcc.deduplicate_fastq(str(input.fq1), output.deduped1, str(input.fq2), output.deduped2)

        with open(log[0], "w") as f:
            json.dump(stats, f)



# rule deduplicate_fastq_flashed:
#     input:
#         flashed="seqnado_output/flashed/{sample}/{sample}.extendedFrags.fastq.gz",
#     output:
#         deduped="seqnado_output/deduped/{sample}/{sample}.extendedFrags.deduped.fastq.gz",
#     threads: 1
#     resources:
#         mem="1GB",
#     log:
#         "seqnado_output/logs/deduplication/{sample}.log",
#     run:
#         from mcc import mcc
#         import json

#         stats = mcc.deduplicate_fastq(str(input.flashed), output.deduped)

#         with open(log[0], "w") as f:
#             json.dump(stats, f)
