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
    tree = "results/panaroo/output/core_gene_alignment_filtered.aln.treefile"
    newick = "results/panaroo/output/outbreak_phylogeny_rectangular.jpg"
    return {
        "tree": newick if config["panaroo"]["newick_tree"] else tree,
        "snps": "results/panaroo/snps_distance/snps_distance_matrix.tsv",
        "summary": "results/summary/summary.tsv",
    }


### Contract for other workflows ######################################################################################


### Parameter parsing from config #####################################################################################


def get_iqtree_bootstrap_param():
    value = config["iqtree"].get("bootstrap", 0)
    if value is None or value == 0:
        return ""
    else:
        return f"-bb {value}"


### Resource handling #################################################################################################


def get_mem_mb_for_XY(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["XY_mem_mb"] * attempt)
