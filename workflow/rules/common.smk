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


if config["validate_taxa_and_mlst"]:
    if "GTDBtk_taxa" not in pep.sample_table or "mlst" not in pep.sample_table:
        raise ValueError(
            "validate_taxa_and_mlst is set to True, but GTDBtk_taxa or mlst are not present in the sample table"
        )
    taxonomies = set(list(pep.sample_table["GTDBtk_taxa"].values))
    mlst_types = set(list(pep.sample_table["mlst"].values))
    if len(taxonomies) > 1:
        raise ValueError(f"Multiple taxonomies found: {taxonomies}")
    if len(mlst_types) > 1:
        raise ValueError(f"Multiple MLST types found: {mlst_types}")


### Global rule-set stuff #############################################################################################


def infer_assembly_fasta(wildcards) -> str:
    return get_fasta_for_sample_from_pep(wildcards.sample)


def get_outputs():
    outputs = {
        "snps": "results/snp_dists/snp_distance_matrix.tsv",
        "summary": "results/summary/summary.tsv",
    }
    tree = "results/panaroo/core_gene_alignment_filtered.aln.treefile"
    newick = "results/panaroo/outbreak_phylogeny_rectangular.jpg"
    outputs["tree"] = newick if config["panaroo"]["newick_tree"] else tree
    if config["panaroo_qc"]["do"]:
        outputs["panaroo_qc"] = "results/panaroo_qc/mash_contamination_barplot.html"
    return outputs


### Contract for other workflows ######################################################################################


### Parameter parsing from config #####################################################################################


def get_iqtree_bootstrap_params():
    value = config["iqtree"]["bootstrap_replicates"]
    if value is None or value == 0:
        return ""

    method = config["iqtree"]["method"]
    if method == "UFBoot":
        boot_arg = f"--ufboot {value}"
    elif method == "standard":
        boot_arg = f"--boot {value}"
    else:
        raise ValueError(f"Unknown bootstrap method: {method}")

    tests_arg = ""
    if sbt_value := config["iqtree"]["single_branch_tests"]:
        if sbt_value == "SH-aLRT":
            tests_arg = f"--alrt {value}"
        elif sbt_value == "abayes":
            tests_arg = "--abayes"
        elif sbt_value == "lbp":
            tests_arg = f"--lbp {value}"
        elif sbt_value != "null":
            raise ValueError(f"Unknown single branch test: {sbt_value}")

    return f"{boot_arg} {tests_arg}"


### Resource handling #################################################################################################


def get_mem_mb_for_panaroo_QC(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["panaroo_qc__mem__mb"] * attempt)


def get_mem_mb_for_panaroo_run(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["panaroo__mem_mb"] * attempt)
