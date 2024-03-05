import os
import sys


def find_coverage(core_size_file: str, lowest_size_file: str, output_file: str):
    with open(core_size_file) as f:
        core_size = int(f.read().strip())
    with open(lowest_size_file) as f:
        lowest_size = int(f.read().strip())

    if not os.path.exists(os.path.dirname(output_file)):
        os.makedirs(os.path.dirname(output_file))

    coverage = core_size / lowest_size

    with open(output_file, "w") as f:
        f.write("Core genome size:\t{}\n".format(core_size))
        f.write("Lowest genome size:\t{}\n".format(lowest_size))
        f.write("Core genome coverage:\t{:.3f}\n".format(coverage))


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    find_coverage(
        snakemake.input.core,
        snakemake.input.lowest,
        snakemake.output[0],
    )
