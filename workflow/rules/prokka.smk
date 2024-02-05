
rule prokka_run:
    input:
        infer_assembly_fasta,
    output:
        multiext("results/prokka/{sample}/{sample}", ".faa", ".ffn", ".fna", ".fsa", ".gbk", ".gff"),
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        # --mincontiglen=200, # Minimum contig size [NCBI needs 200] (default '1')
        # --centre="CUSP", # Sequencing centre ID. (default '')
        # --compliant       Force Genbank/ENA/DDJB compliance:
        extra="--mincontiglen=200 --centre=CUSP --compliant --noanno",
    conda:
        "../envs/prokka.yaml"
    threads: min(config["threads"]["prokka"], config["max_threads"])
    log:
        "logs/prokka/{sample}.log",
    shell:
        "prokka --force --prefix {wildcards.sample} --cpus {threads} --outdir {params.outdir} {input} > {log} 2>&1"


rule roary_run:
    input:
        GFFs=expand("results/prokka/{sample}/{sample}.gff", sample=get_sample_names()),
    output:
        aln="results/roary/_1705504897/core_gene_alignment.aln",
    params:
        outdir=lambda wildcards, output: os.path.dirname(os.path.dirname(output[0])),
    conda:
        "../envs/roary.yaml"
    threads: min(config["threads"]["roary"], config["max_threads"])
    log:
        "logs/roary.log",
    shell:
        "roary -r -e --mafft -p {threads} {input} -f {output} > {log} 2>&1"


rule snpdists_compute:
    input:
        aln="results/roary/_1705504897/core_gene_alignment.aln",
    output:
        "results/roary/_1705504897/snps_distance_matrix.tsv",
    conda:
        "../envs/snpdists.yaml"
    threads: min(config["threads"]["snp_dists"], config["max_threads"])
    log:
        "logs/snp_dists.log",
    shell:
        '(snp-dists -b -j {threads} {input} | sed -e "s/_Unicycler_scaffolds_annot//g" > {output})  2> {log}'


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


rule iqtree_phylogeny:
    input:
        aln="results/roary/_1705504897/core_gene_alignment.aln",
    output:
        tree="results/roary/_1705504897/core_gene_alignment.aln.treefile",
    log:
        "logs/iqtree.log",
    threads: min(config["threads"]["iqtree"], config["max_threads"])
    conda:
        "../envs/iqtree.yaml"
    shell:
        "iqtree -s {input.aln} -T {threads} > {log} 2>&1"


# phylogeny
# sed -e "s/_Unicycler_scaffolds_annot//g" roary/*/core_gene_alignment.aln.treefile > cga_IQtree.newick
# coverage of genome calculation
# core_genome_size=$(tail -n 1 snps_distance.log | rev | cut -d ' ' -f 1 | rev)
# coverage=$(echo "scale=3;  $core_genome_size / $lowest_genome_size" | bc)
# echo "Core genome coverage of the computed genomes is 0${coverage}" > core_genome_coverage.txt
# echo "Core genome size is ${core_genome_size} bp." >> core_genome_coverage.txt
# echo "Genome with smallest size is ${lowest_genome_size} bp long. This includes noncoding regions and paralogous genes." >> core_genome_coverage.txt
