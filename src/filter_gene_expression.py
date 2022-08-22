"""
Filter the gene expression tables so that I can analyze them easier in
conjunction with the TE data
"""

__author__ = "Scott Teresi"

import pandas as pd
import argparse
import os
import logging
import coloredlogs


def clean_tpm(input_tpm_file, genome_name):
    """
    Take a TPM file and turn it into a pandas dataframe
    """
    data = pd.read_csv(input_tpm_file, sep="\t", header="infer")
    data.rename(columns={"ID": "Gene_Name"}, inplace=True)
    chromosomes_i_want = [
        "Chr01B",
        "Chr01A",
        "Chr02B",
        "Chr02A",
        "Chr03B",
        "Chr03A",
        "Chr04B",
        "Chr04A",
        "Chr05B",
        "Chr05A",
        "Chr06B",
        "Chr06A",
        "Chr07B",
        "Chr07A",
        "Chr08B",
        "Chr08A",
        "Chr09B",
        "Chr09A",
        "Chr10B",
        "Chr10A",
        "Chr11B",
        "Chr11A",
        "Chr12B",
        "Chr12A",
        "Chr13B",
        "Chr13A",
        "Chr14B",
        "Chr14A",
        "Chr15B",
        "Chr15A",
        "Chr16B",
        "Chr16A",
        "Chr17B",
        "Chr17A",
        "Chr18B",
        "Chr18A",
        "Chr19B",
        "Chr19A",
        "Chr20B",
        "Chr20A",
        "Chr21B",
        "Chr21A",
        "Chr22B",
        "Chr22A",
        "Chr23B",
        "Chr23A",
        "Chr24B",
        "Chr24A",
        "Chr25B",
        "Chr25A",
    ]
    data["Gene_Name"] = data["Gene_Name"].str.rstrip(".m1")

    # MAGIC, get the pseudomolecule name by splitting on the first period in the
    # gene name
    data["Chromosome"] = data["Gene_Name"].str.split(".").str[0]

    # Add the genome name as string
    data["Genome"] = genome_name
    # Filter on my chromosome whitelist
    data = data.loc[data["Chromosome"].isin(chromosomes_i_want)]
    data.sort_values(by=["Chromosome"], inplace=True)

    return data


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Read gene expression files as pandas")
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser.add_argument("LC_TPM", type=str, help="LC gene expression file")
    parser.add_argument("PR_TPM", type=str, help="PR gene expression file")
    parser.add_argument("SS_TPM", type=str, help="SS gene expression file")

    parser.add_argument(
        "output_dir",
        type=str,
        help="Parent directory to output results",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.LC_TPM = os.path.abspath(args.LC_TPM)
    args.PR_TPM = os.path.abspath(args.PR_TPM)
    args.SS_TPM = os.path.abspath(args.SS_TPM)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Execute
    LC_TPM = clean_tpm(args.LC_TPM, "LC")
    PR_TPM = clean_tpm(args.PR_TPM, "PR")
    SS_TPM = clean_tpm(args.SS_TPM, "SS")

    LC_TPM.to_csv(
        os.path.join(args.output_dir, "LC_TPM_Filtered.tsv"),
        sep="\t",
        header=True,
        index=False,
    )
    PR_TPM.to_csv(
        os.path.join(args.output_dir, "PR_TPM_Filtered.tsv"),
        sep="\t",
        header=True,
        index=False,
    )
    SS_TPM.to_csv(
        os.path.join(args.output_dir, "SS_TPM_Filtered.tsv"),
        sep="\t",
        header=True,
        index=False,
    )
