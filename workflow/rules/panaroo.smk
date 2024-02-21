
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


rule panaroo_download_mash_db:
    output:
        os.path.join(config["panaroo"]["mash_db"], "refseq.genomes.k21s1000.msh"),
    params:
        url="https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh",
    conda:
        "../envs/panaroo.yaml"
    log:
        os.path.join(config["panaroo"]["mash_db"], "logs", "mash_db.log"),
    shell:
        "wget -O {output} {params.url} > {log} 2>&1"


rule panaroo_QC:
    input:
        GFFs=expand("results/prokka/{sample}/{sample}.gff", sample=get_sample_names()),
        mash_db=os.path.join(config["panaroo"]["mash_db"], "refseq.genomes.k21s1000.msh"),
    output:
        "results/panaroo/QC/contam.graph",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
    conda:
        "../envs/panaroo.yaml"
    threads: min(config["threads"]["panaroo"], config["max_threads"])
    log:
        "logs/panaroo_qc.log",
    shell:
        "(rm -rf {params.outdir} && panaroo-qc -t {threads} --graph_type all -i {input.GFFs} --ref_db {input.mash_db} -o {params.outdir}) > {log} 2>&1"


rule panaroo_run:
    input:
        GFFs=expand("results/prokka/{sample}/{sample}.gff", sample=get_sample_names()),
        panaroo_qc="results/panaroo_QC/contam.graph",
    output:
        aln="results/panaroo/output/core_gene_alignment.aln",
        aln_filt="results/panaroo/output/core_gene_alignment_filtered.aln",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        mode="sensitive",
        core_threshold=0.98,
    conda:
        "../envs/panaroo.yaml"
    threads: min(config["threads"]["panaroo"], config["max_threads"])
    log:
        "logs/panaroo.log",
    shell:
        "(rm -rf {params.outdir} && panaroo -i {input.GFFs} -o {params.outdir} -a core --aligner mafft"
        " --core_threshold {params.core_threshold} -t {threads} --clean-mode {params.mode}) > {log} 2>&1"


# rule panaroo_postfilter:
#     panaroo-filter-pa -i ./gene_presence_absence.csv -o ./ --type pseudo,length


rule snpdists_compute:
    input:
        aln="results/panaroo/output/core_gene_alignment_filtered.aln",
    output:
        "results/panaroo/snps_distance/snps_distance_matrix.tsv",
    conda:
        "../envs/snpdists.yaml"
    threads: min(config["threads"]["snp_dists"], config["max_threads"])
    log:
        "logs/snp_dists.log",
    shell:
        '(snp-dists -b -j {threads} {input} | sed -e "s/_Unicycler_scaffolds_annot//g" > {output})  2> {log}'


rule iqtree_phylogeny:
    input:
        aln="results/panaroo/output/core_gene_alignment_filtered.aln",
    output:
        tree="results/panaroo/phylogeny/core_gene_alignment_filtered.aln.treefile",
    log:
        "logs/iqtree.log",
    threads: min(config["threads"]["iqtree"], config["max_threads"])
    conda:
        "../envs/iqtree.yaml"
    shell:
        "iqtree -s {input.aln} -nt {threads} > {log} 2>&1"


# phylogeny
# sed -e "s/_Unicycler_scaffolds_annot//g" roary/*/core_gene_alignment.aln.treefile > cga_IQtree.newick
# coverage of genome calculation
# core_genome_size=$(tail -n 1 snps_distance.log | rev | cut -d ' ' -f 1 | rev)
# coverage=$(echo "scale=3;  $core_genome_size / $lowest_genome_size" | bc)
# echo "Core genome coverage of the computed genomes is 0${coverage}" > core_genome_coverage.txt
# echo "Core genome size is ${core_genome_size} bp." >> core_genome_coverage.txt
# echo "Genome with smallest size is ${lowest_genome_size} bp long. This includes noncoding regions and paralogous genes." >> core_genome_coverage.txt
# # needs conda activate R_new environment
# # simple tree drawing
# Rscript --vanilla newick_tree_plott.R cga_IQtree.newick outbreak_phylogeny_rectangular.jpg
# echo 'conda activate R_new'
# echo 'Rscript --vanilla newick_tree_plott.R cga_IQtree.newick outbreak_phylogeny_rectangular.jpg'
