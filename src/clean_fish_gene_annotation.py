"""
Filter a gene annotation file for the TE Density algorithm
"""

__author__ = "Scott Teresi"

import pandas as pd
import argparse
import os
import logging
import coloredlogs


def import_genes(genes_input_path):
    """Import genes file.

    Args:
        TODO
    """

    col_names = [
        "Chromosome",
        "Software",
        "Feature",
        "Start",
        "Stop",
        "Score",
        "Strand",
        "Frame",
        "FullName",
    ]

    col_to_use = [
        "Chromosome",
        "Software",
        "Feature",
        "Start",
        "Stop",
        "Strand",
        "FullName",
    ]

    gene_anno = pd.read_csv(
        genes_input_path,
        sep="\t+",
        header=None,
        engine="python",
        names=col_names,
        usecols=col_to_use,
        dtype={"Stop": "float64", "Start": "float64", "Strand": str, "FullName": str},
        comment="#",
    )

    # TODO check column names

    # rows in annotation
    gene_anno = gene_anno[gene_anno.Feature == "gene"]  # drop non-gene rows

    # clean the names and set as the index (get row wrt name c.f. idx)
    gene_anno["Gene_Name"] = gene_anno["FullName"].str.extract(r"ID=?(.*)")
    gene_anno.set_index("Gene_Name", inplace=True)

    # Drop extraneous columns
    gene_anno = gene_anno.drop(columns=["FullName", "Software"])

    # Edit gene length
    gene_anno["Length"] = gene_anno.Stop - gene_anno.Start + 1

    gene_anno.sort_values(by=["Chromosome", "Start"], inplace=True)
    return gene_anno


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Reformat gene annotation file")
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    output_default = os.path.join(dir_main, "../results")
    parser.add_argument(
        "gene_input_file", type=str, help="Parent path of gene annotation file"
    )

    parser.add_argument("genome_name", type=str, help="Name for the genome")

    parser.add_argument(
        "output_dir",
        type=str,
        help="Parent directory to output results",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.gene_input_file = os.path.abspath(args.gene_input_file)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Execute
    cleaned_genes = import_genes(args.gene_input_file)

    # Save the file
    file_name = os.path.join(
        args.output_dir, ("Cleaned_" + args.genome_name + "_Genes.tsv")
    )

    logger.info("Writing cleaned gene file to: %s" % file_name)
    cleaned_genes.to_csv(file_name, sep="\t", header=True, index=True)
