rule seqkit_stats:
    input:
        fasta=get_all_assemblies(),
    output:
        stats="results/summary/seqkit_stats.tsv",
    log:
        "logs/summary/seqkit.log",
    params:
        command="stats",
        extra="--tabular",
    wrapper:
        "v3.3.6/bio/seqkit"


rule compute_lowest_genome_size:
    input:
        stats="results/summary/seqkit_stats.tsv",
    output:
        stats="results/summary/lowest_genome_size.txt",
    log:
        "logs/summary/lowest_genome_size.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        "(cat {input} | tail -n +2 | cut -f 5 | sort | head -n 1 > {output}) 2> {log}"
