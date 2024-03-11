rule prokka_run:
    input:
        infer_assembly_fasta,
    output:
        multiext("results/prokka/{sample}/{sample}", ".faa", ".ffn", ".fna", ".fsa", ".gbk", ".gff"),
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        mincontiglen=config["prokka"]["mincontiglen"],
        centre=config["prokka"]["centre"],
    conda:
        "../envs/prokka.yaml"
    threads: min(config["threads"]["prokka"], config["max_threads"])
    log:
        "logs/prokka/{sample}.log",
    shell:
        "prokka --force --prefix {wildcards.sample} --noanno --compliant --mincontiglen {params.mincontiglen}"
        " --centre {params.centre} --cpus {threads} --outdir {params.outdir} {input} > {log} 2>&1"


rule panaroo_download_mash_db:
    output:
        protected(os.path.join(config["panaroo_qc"]["mash_db"], "refseq.genomes.k21s1000.msh")),
    params:
        url="https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh",
    conda:
        "../envs/panaroo.yaml"
    log:
        os.path.join(config["panaroo_qc"]["mash_db"], "logs", "mash_db.log"),
    shell:
        "wget -O {output} {params.url} > {log} 2>&1"


rule panaroo_qc:
    input:
        GFFs=expand("results/prokka/{sample}/{sample}.gff", sample=get_sample_names()),
        mash_db=os.path.join(config["panaroo_qc"]["mash_db"], "refseq.genomes.k21s1000.msh"),
    output:
        "results/panaroo_qc/mash_contamination_barplot.html",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
    conda:
        "../envs/panaroo.yaml"
    threads: min(config["threads"]["panaroo_QC"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_panaroo_QC,
    log:
        "logs/panaroo_qc.log",
    shell:
        "(rm -rf {params.outdir} && mkdir -p {params.outdir} && panaroo-qc -t {threads} --graph_type all -i {input.GFFs} --ref_db {input.mash_db} -o {params.outdir}) > {log} 2>&1"


rule panaroo_run:
    input:
        GFFs=expand("results/prokka/{sample}/{sample}.gff", sample=get_sample_names()),
    output:
        aln="results/panaroo/output/core_gene_alignment.aln",
        aln_filt="results/panaroo/output/core_gene_alignment_filtered.aln",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        aligner=config["panaroo"]["aligner"],
        mode=config["panaroo"]["mode"],
        core_threshold=config["panaroo"]["core_threshold"],
    conda:
        "../envs/panaroo.yaml"
    threads: min(config["threads"]["panaroo"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_panaroo_run,
    log:
        "logs/panaroo.log",
    shell:
        "(rm -rf {params.outdir} && panaroo -i {input.GFFs} -o {params.outdir} -a core --aligner {params.aligner}"
        " --core_threshold {params.core_threshold} -t {threads} --clean-mode {params.mode}) > {log} 2>&1"


rule snpdists_compute:
    input:
        aln="results/panaroo/output/core_gene_alignment_filtered.aln",
    output:
        tsv="results/panaroo/snps_distance/snps_distance_matrix.tsv",
        log="results/panaroo/snps_distance/snps_distance.log",
    conda:
        "../envs/snpdists.yaml"
    threads: min(config["threads"]["snp_dists"], config["max_threads"])
    log:
        "logs/snp_dists.log",
    shell:
        '(snp-dists -b -j {threads} {input} | sed -e "s/_Unicycler_scaffolds_annot//g" > {output.tsv})  2> {output.log}'
        ' && echo "snpdists log redirected into {output.log}" > {log}'


rule iqtree_phylogeny:
    input:
        aln="results/panaroo/output/core_gene_alignment_filtered.aln",
    output:
        tree="results/panaroo/output/core_gene_alignment_filtered.aln.treefile",
    params:
        bootstrap=get_iqtree_bootstrap_params(),
    log:
        "logs/iqtree.log",
    threads: min(config["threads"]["iqtree"], config["max_threads"])
    conda:
        "../envs/iqtree.yaml"
    shell:
        "iqtree2 -s {input.aln} -T {threads} {params.bootstrap} > {log} 2>&1"


rule find_core_genome_size:
    input:
        log="results/panaroo/snps_distance/snps_distance.log",
    output:
        temp("results/summary/core_genome_size.txt"),
    log:
        "logs/find_core_genome_size.log",
    conda:
        "../envs/coreutils.yaml"
    localrule: True
    shell:
        "(tail -n 1 {input.log} | rev | cut -d ' ' -f 1 | rev) > {output} 2> {log}"


rule find_coverage_and_summarize:
    input:
        core="results/summary/core_genome_size.txt",
        lowest="results/summary/lowest_genome_size.txt",
    output:
        "results/summary/summary.tsv",
    conda:
        "../envs/python.yaml"
    localrule: True
    log:
        "logs/find_coverage_and_summarize.log",
    script:
        "../scripts/find_coverage_and_summarize.py"


rule visualize_tree:
    input:
        tree="results/panaroo/output/core_gene_alignment_filtered.aln.treefile",
    output:
        tree_img="results/panaroo/output/outbreak_phylogeny_rectangular.jpg",
    localrule: True
    conda:
        "../envs/newick_plot.yaml"
    log:
        "logs/visualize_tree.log",
    script:
        "../scripts/plot_newick_tree.R"
