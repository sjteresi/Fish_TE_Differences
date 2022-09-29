"""
Analyze the TE density values of 'regular' syntelogs in various genomes
"""

__author__ = "Scott Teresi"

import pandas as pd
import argparse
import os
import logging
import coloredlogs
import numpy as np
import re
from configparser import ConfigParser

from transposon.import_filtered_genes import import_filtered_genes
from transposon.density_data import DensityData
from transposon.gene_data import GeneData

from src.fish_utils import (
    get_gene_data_as_list,
    add_hdf5_indices_to_gene_data_from_list_hdf5,
    add_te_vals_to_gene_info_pandas_from_list_hdf5,
    add_te_vals_to_syntelogs,
    parse_analysis_config,
)

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def read_syntelog_file(filepath, genome_name):
    """
    Read the syntelog file provided by collaborators as a pandas.DataFrame

    Args:
        filepath (str): String representing the filepath on disk.

        genome_name (str): String repreenting the shorthand genome name. This
        is used to subset the file with the appropriate columns so that we get
        the right A and B subgenomes. Acceptable values provided by
        args.genome_name

    Returns:
            syntelogs (pandas.DataFrame)
                Index:
                    RangeIndex
                Columns:
                    Name: genome_name + '_A': dtype, object
                    Name: genome_name + '_B', dtype: object
    """
    # MAGIC column names inherent to the file
    syntelogs = pd.read_csv(
        filepath,
        header=None,
        sep="\t",
        names=[
            "Drop",
            "LC_B",
            "LC_A",
            "Drop2",
            "PR_B",
            "PR_A",
            "Drop3",
            "SS_A",
            "SS_B",
        ],
        usecols=["LC_B", "LC_A", "PR_B", "PR_A", "SS_A", "SS_B"],
    )
    cols_i_want = [genome_name + "_A", genome_name + "_B"]
    syntelogs = syntelogs[cols_i_want]

    # Remove the '.m1' string from the end
    # NOTE, from now on use the replace syntax because the .strip can have
    # weird unintended functionality
    syntelogs = syntelogs.apply(lambda x: x.str.replace("\.m1", "", regex=True))

    # Remove the genome name prefix, this is inherent to the input data
    # MAGIC string split
    syntelogs = syntelogs.apply(lambda x: x.str.split("_").str[1])

    return syntelogs


def graph_barplot_density_differences(
    values,
    te_type,
    window_val,
    direction,
    number_of_zeros,
    genome_name_str,
    genome_1,
    genome_2,
    output_dir,
    logger,
    display=False,
    align="left",
):
    """
    Plot a histogram of TE density differences between syntelog pairs

    Args:
        values (list): A list of values representing the TE density differences
        between syntelog pairs

        te_type (str): String representing the TE type being plotted

        window_val (int): Integer representing the current window of which the
        data is being plotted

        direction (str): string representing whether or not the graphs are
        coming from upstream or downstream TE density data

        number_of_zeros ():

        logger (logging.Logger): Object to log information to

        display (boolean): Defaults to False, if True shows the plot upon
        generation with the plt.show() command
    """

    # MAGIC, the bins to group density values for the histogram AND the values
    # for the xticks on the xaxis

    # NOTE could use linspace but I want the numbers to show up nice for the
    # figure, so I had to hard-code
    # tick_bins = np.linspace(-1.0, 1.0, num=21, endpoint=True).tolist()
    tick_bins = [
        -1.0,
        -0.9,
        -0.8,
        -0.7,
        -0.6,
        -0.5,
        -0.4,
        -0.3,
        -0.2,
        -0.1,
        0,
        0.1,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.8,
        0.9,
        1.0,
    ]

    plt.figure(figsize=(8, 6))
    n, bins, patches = plt.hist(
        values, bins=tick_bins, facecolor="blue", ec="black", alpha=0.5
    )
    plt.rcParams["xtick.labelsize"] = 7  # MAGIC set size of axis ticks
    plt.ylabel("Number of Genes")
    plt.xlabel("Difference in TE Density Values")

    plt.title(
        genome_1.replace("_", " Subgenome ")
        + " - "
        + genome_2.replace("_", " Subgenome ")
    )  # MAGIC genome name order here, Be
    # careful that it corresponds to your subtraction scheme

    N = mpatches.Patch(
        color=None,
        alpha=0,
        label="Total Plotted Genes: %s \nTE type: %s \nWindow: %s \nDirection: %s \nNo. 0 Differences: %s"
        % (len(values), te_type, window_val, direction, str(number_of_zeros)),
    )
    plt.xticks(tick_bins)
    plt.legend(handles=[N])

    # Create a black vertical line at 0 for visual purposes
    plt.axvline(x=0, color="black")

    path = os.path.join(
        output_dir,
        "Barplot",
        (
            genome_name_str
            + "_"
            + te_type
            + "_"
            + str(window_val)
            + "_"
            + direction
            + "_DensityDifferences.svg"
        ),
    )
    logger.info("Saving graph to: %s" % path)
    os.makedirs(os.path.join(args.output_dir, "Barplot"), exist_ok=True)
    plt.savefig(path, bbox_inches="tight")
    if display:
        plt.show()
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "syntelog_file",
        type=str,
        help="""syntelog file provided
                        by collaborators""",
    )

    parser.add_argument(
        "genome_name", type=str, help="genome name", choices=["LC", "SS", "PR"]
    )
    parser.add_argument(
        "gene_data",
        type=str,
        help="""gene data file that is the input to TE Density""",
    )

    parser.add_argument(
        "HDF5_dir",
        type=str,
        help=""".h5 files corresponding to the genome""",
    )

    parser.add_argument(
        "config_file",
        type=str,
        help="parent path of config file",
    )

    parser.add_argument(
        "output_dir",
        type=str,
        help="Parent directory to output results",
    )
    logger = logging.getLogger(__name__)
    args = parser.parse_args()
    args.syntelog_file = os.path.abspath(args.syntelog_file)
    args.gene_data = os.path.abspath(args.gene_data)
    args.config_file = os.path.abspath(args.config_file)
    args.HDF5_dir = os.path.abspath(args.HDF5_dir)
    args.output_dir = os.path.abspath(args.output_dir)
    windows, directions, orders = parse_analysis_config(args.config_file)

    # -----------------------------------------------------------------------
    # Read and initialize various data

    # Read the universal syntelog table
    syntelogs = read_syntelog_file(args.syntelog_file, args.genome_name)

    # Read cleaned genes for the given genome as pandas
    cleaned_genes = import_filtered_genes(args.gene_data, logger)

    # Get list of GeneData for each genome to enable initialization of
    # DensityData
    genedata_list = get_gene_data_as_list(cleaned_genes)

    # Initialize DensityData for each genome
    processed_dd_data = DensityData.from_list_gene_data_and_hdf5_dir(
        genedata_list, args.HDF5_dir, args.genome_name + "_(.*?).h5", logger
    )

    # Reset index to make it easier to add the HDF5 indices to a pandas frame
    cleaned_genes.reset_index(inplace=True)

    # Add HDF5 indices to a pandas dataframe to enable easier HDF5 TE value
    # access later
    # This is used for the special and regular gene set independently
    gene_frame_with_indices = add_hdf5_indices_to_gene_data_from_list_hdf5(
        cleaned_genes, processed_dd_data
    )

    # -----------------------------------------------------------------------
    # Begin analysis of regular genes

    # Do the loops to create graphs for the general set
    for window in windows:
        for direction in directions:
            for order in orders:

                # NOTE this is only for the orders. I could do it for the
                # superfamilies but I don't think that is wanted or needed
                cleaned_with_te_vals = add_te_vals_to_gene_info_pandas_from_list_hdf5(
                    gene_frame_with_indices,
                    processed_dd_data,
                    "Order",
                    order,
                    direction,
                    window,
                )

                syntelogs_w_te_vals = add_te_vals_to_syntelogs(
                    syntelogs, cleaned_with_te_vals, args.genome_name
                )

                # MAGIC string formatting
                # NOTE A genome MINUS B genome. Negative values means that the
                # B copy had more TE
                syntelogs_w_te_vals["Difference"] = (
                    syntelogs_w_te_vals[
                        order + "_" + str(window) + "_" + direction + "_A"
                    ]
                    - syntelogs_w_te_vals[
                        order + "_" + str(window) + "_" + direction + "_B"
                    ]
                )

                total_length = len(syntelogs_w_te_vals)

                # NB subset to have only rows with a nonzero difference

                syntelogs_w_te_vals = syntelogs_w_te_vals.loc[
                    syntelogs_w_te_vals["Difference"] != 0
                ]
                number_of_zeros = total_length - len(syntelogs_w_te_vals)

                graph_barplot_density_differences(
                    syntelogs_w_te_vals["Difference"].to_list(),
                    order,
                    window,
                    direction,
                    number_of_zeros,
                    args.genome_name,
                    args.genome_name + "_A",
                    args.genome_name + "_B",
                    args.output_dir,
                    logger,
                    display=False,
                )
