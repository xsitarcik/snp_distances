from snakemake.utils import validate


configfile: "config/config.yaml"


validate(config, "../schemas/config.schema.yaml", set_default=False)


pepfile: config.get("pepfile", "config/pep/config.yaml")


validate(pep.sample_table, "../schemas/samples.schema.yaml")


### Layer for adapting other workflows  ###############################################################################


### Data input handling independent of wildcards ######################################################################


def get_sample_names() -> list[str]:
    return list(pep.sample_table["sample_name"].values)


def get_fasta_for_sample_from_pep(sample: str) -> str:
    return pep.sample_table.loc[sample][["fasta"]][0]


def get_all_assemblies() -> list[str]:
    return list(pep.sample_table["fasta"].values)


### Global rule-set stuff #############################################################################################


def infer_assembly_fasta(wildcards) -> str:
    return get_fasta_for_sample_from_pep(wildcards.sample)


def get_outputs():
    return {
        "stats": "results/summary/lowest_genome_size.txt",
        "tree": "results/panaroo/phylogeny/core_gene_alignment_filtered.aln.treefile",
        "snps": "results/panaroo/snps_distance/snps_distance_matrix.tsv",
        "summary": "results/summary/summary.txt",
    }


### Contract for other workflows ######################################################################################


### Parameter parsing from config #####################################################################################


### Resource handling #################################################################################################


def get_mem_mb_for_XY(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["XY_mem_mb"] * attempt)
