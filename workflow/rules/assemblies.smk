rule seqkit_stats:
    input:
        fasta=get_all_assemblies(),
    output:
        stats="results/summary/seqkit_stats.tsv",
    log:
        "logs/summary/seqkit.log",
    params:
        extra="--tabular",
    conda:
        "../envs/seqkit.yaml"
    shell:
        "seqkit stats {input} {params.extra} > {output} 2> {log}"


rule compute_lowest_genome_size:
    input:
        stats="results/summary/seqkit_stats.tsv",
    output:
        stats=temp("results/summary/lowest_genome_size.txt"),
    log:
        "logs/summary/lowest_genome_size.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        "(cat {input} | tail -n +2 | cut -f 5 | sort -n | head -n 1 > {output}) 2> {log}"
