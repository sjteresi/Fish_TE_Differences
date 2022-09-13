"""
Analyze the TE density values of syntelogs in various genomes
"""

__author__ = "Scott Teresi"

import pandas as pd
from collections import namedtuple
import argparse
import os
import logging
import coloredlogs
import numpy as np
import re

from transposon.import_filtered_genes import import_filtered_genes
from transposon.density_data import DensityData
from transposon.gene_data import GeneData

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import scipy
from statannotations.Annotator import Annotator
import statsmodels
from collections import OrderedDict


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
            + "_DensityDifferences.png"
        ),
    )
    logger.info("Saving graph to: %s" % path)
    os.makedirs(os.path.join(args.output_dir, "Barplot"), exist_ok=True)
    plt.savefig(path, bbox_inches="tight")
    if display:
        plt.show()
    plt.close()


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


def get_gene_data_as_list(cleaned_genes):
    """
    Take a cleaned genes annotation file from TE Density
    (import_filtered_genes) and break it into a list of GeneData objects by
    chromosome ID. This is used to initialize all of the DensityData objects in
    a list.

    Args:
        cleaned_genes (pandas.DataFrame)
            Index:
                Name: Gene_Name, strings of gene names
            Columns:
                Name: Chromosome, object
                Name: Feature, object
                Name: Start, float64
                Name: Stop, float64
                Name: Strand, object
                Name: Length, float64

    Returns:
        genedata_list (list of GeneData)
    """
    # MAGIC group by column Chromosome
    gene_dataframe_list = [
        dataframe for k, dataframe in cleaned_genes.groupby("Chromosome")
    ]

    # MAGIC initialize GeneData iteratively using the magic unique chromosome
    genedata_list = [
        GeneData(dataframe, dataframe["Chromosome"].unique()[0])
        for dataframe in gene_dataframe_list
    ]
    return genedata_list


def add_hdf5_indices_to_gene_data_from_list_hdf5(cleaned_genes, list_processed_dd_data):
    """
    NOTE this could be refactored to DensityData... Others may find it useful
    """
    to_concat = []
    for chrom, dataframe in cleaned_genes.groupby(["Chromosome"]):
        for processed_dd_datum in list_processed_dd_data:
            if processed_dd_datum.unique_chromosome_id == chrom:
                x = processed_dd_datum.add_hdf5_indices_to_gene_data(dataframe)
                to_concat.append(x)
    gene_frame_with_indices = pd.concat(to_concat)
    return gene_frame_with_indices


def add_te_vals_to_gene_info_pandas_from_list_hdf5(
    gene_frame_with_indices,
    list_processed_dd_data,
    te_group,
    te_name,
    direction,
    window,
    gene_name_col="Gene_Name",
    chrom_col="Chromosome",
    index_col="Index_Val",
):
    """
    NOTE this could be refactored to DensityData...
    This is a really bad 'wrapper' for a function.

    """
    to_concat = []
    for chrom, dataframe in gene_frame_with_indices.groupby([chrom_col]):
        for processed_dd_datum in list_processed_dd_data:
            if processed_dd_datum.unique_chromosome_id == chrom:
                x = processed_dd_datum.add_te_vals_to_gene_info_pandas(
                    dataframe,
                    te_group,
                    te_name,
                    direction,
                    window,
                    gene_name_col,
                    chrom_col,
                    index_col,
                )
                to_concat.append(x)
    gene_frame_w_ind_te_vals = pd.concat(to_concat)

    return gene_frame_w_ind_te_vals


def read_special_gene_set(filepath, genome_name):
    """
    Read a csv file received from Kevin Bird that is a special set of syntelogs
    that are biased in expression towards one subgenome or the other.

    Args:
        filepath (str): String representing the filepath to the 'raw' data,
        provided by collaborators

        genome_name (str): String repreenting the shorthand genome name. This
        is used to subset the file with the appropriate columns so that we get
        the right A and B subgenomes. Acceptable values provided by
        args.genome_name

    Returns:
        special_syntelogs (pandas.Data.Frame)
            Index:
                RangeIndex
            Columns:
                Name: Sub_A_Gene, dtype: object
                Name: Sub_B_Gene, dtype: object
                Name: Bias, dtype: object
                Name: Tissue, dtype: object
                Name: Species, dtype: object

        tissue_list (list): List of strings. The strings are tissues and are
        are where gene a or B was biased
    """
    special_syntelogs = pd.read_csv(filepath, header="infer", sep=",")
    special_syntelogs = special_syntelogs.iloc[:, 1:]  # MAGIC iloc command to
    # get rid of the first column which is just row numbers

    # Get a list of tissues from the dataset, this will be used to subset out
    # genes later. It is important that we do it at the beginning so we don't
    # miss any tissues later after we subset.
    tissue_list = special_syntelogs["Tissue"].unique().tolist()

    # Replace values with the shortened gene name
    special_syntelogs["Species"].replace(
        {
            "Spinibarbus sinensis": "SS",
            "Luciobarbus capito": "LC",
            "Procypris rabaudi": "PR",
        },
        inplace=True,
    )
    special_syntelogs = special_syntelogs.loc[
        special_syntelogs["Species"] == genome_name
    ]
    special_syntelogs.drop(["L2FC"], axis=1, inplace=True)

    # Remove the '.m1' string from the end, same as the other syntelog read
    # func
    special_syntelogs[["SubA", "SubB"]] = special_syntelogs[["SubA", "SubB"]].apply(
        lambda x: x.str.replace("\.m1", "", regex=True)
    )

    # Remove the genome name prefix, same as the other syntelog read func
    # MAGIC column names for subgenome specificity and the string split to get
    # the portion after the '_'
    special_syntelogs[["SubA", "SubB"]] = special_syntelogs[["SubA", "SubB"]].apply(
        lambda x: x.str.split("_").str[1]
    )

    # Sort for legibility, MAGIC
    special_syntelogs = special_syntelogs.sort_values(by=["SubA", "SubB"])

    # Rename columns to be easier to understand and work with input args
    # MAGIC
    special_syntelogs.loc[
        special_syntelogs["Bias"] == "SubA Biased", "Bias"
    ] = "Sub_A_Biased"

    special_syntelogs.loc[
        special_syntelogs["Bias"] == "SubB Biased", "Bias"
    ] = "Sub_B_Biased"

    # Rename columns to be easier to understand and work with input args
    special_syntelogs.rename(
        columns={"SubA": "Sub_A_Gene", "SubB": "Sub_B_Gene"}, inplace=True
    )

    return special_syntelogs, tissue_list


def gen_boxplot(
    cleaned_with_te_vals,
    bias_str,
    other,
    order,
    window,
    direction,
    te_column_str,
    n_tissues,
    long_genome_name,
    p_val,
    number_genes,
    output_dir,
):
    """
    Create a boxlot of the nonzero TE values of the A and B subgenomes and test
    for differences

    Args:
        cleaned_with_te_vals

        bias_str

        other

        order

        window

        direction

        te_column_str

        n_tissues

        long_genome_name

        p_val

        number_genes,

        output_dir

    Returns: None, generates a boxplot and saves to disk
    """

    te_col_other = str(te_column_str + "_" + other)
    te_col_main = str(te_column_str + "_" + bias_str)

    # Subset out the nonzero values
    cleaned_with_te_vals = cleaned_with_te_vals[
        (cleaned_with_te_vals[te_col_other] != 0)
        & (cleaned_with_te_vals[te_col_main] != 0)
    ]

    # Return out of this if we detect that there is empty dataframe
    # This would happen if we don't have any nonzero TE values for the genes
    if cleaned_with_te_vals.empty:
        print(f"{te_column_str} yields an empty dataframe: {cleaned_with_te_vals}")
        return None

    # TODO rename variable
    melted_genes_w_te_data = pd.melt(
        cleaned_with_te_vals,
        id_vars=[bias_str + "_Gene", other + "_Gene"],
        value_vars=[te_col_main, te_col_other],
    )
    melted_genes_w_te_data.sort_values(
        by=["variable"], inplace=True
    )  # NOTE sort the values so
    # that A always comes before B so it is nice when we graph

    # Set colors and order for consistency across loops of bias (A or B)
    # Collaborators want A to always be RED
    if "_B" in bias_str:
        bias_color = "tab:blue"
        other_color = "tab:red"
        order_of_things = [te_col_other, te_col_main]
    if "_A" in bias_str:
        bias_color = "tab:red"
        other_color = "tab:blue"
        order_of_things = [te_col_main, te_col_other]

    plt.figure(figsize=(6, 8))  # MAGIC figure size
    ax = sns.boxplot(
        data=melted_genes_w_te_data,
        x="variable",
        y="value",
        order=order_of_things,
        palette={te_col_main: bias_color, te_col_other: other_color},
    )

    # Make the annotator for significance
    annotator = Annotator(
        ax,
        [(te_col_main, te_col_other)],
        data=melted_genes_w_te_data,
        x="variable",
        y="value",
        order=order_of_things,
        palette={te_col_main: bias_color, te_col_other: other_color},
    )

    # NOTE use this if you want the actual p-val displayed
    # Transform each p-value to "p=" in scientific notation
    # formatted_pvalue = f"p={p_val:.2e}"
    # annotator.custom_annotations([formatted_pvalue])

    # NOTE this is to show ns, *, **, etc for significance. Does not show
    # actual p-value
    annotator.configure(text_format="star")  # NB if you want actual val use
    # 'simple'
    annotator.set_pvalues([p_val])
    annotator.annotate()
    plt.ylabel("TE Density Values")
    plt.ylim((0, 1.19))  # NB technically we cant get values over 1 but making
    # the graph taller helps with the significance line not looking bad
    locs, labels = plt.xticks()

    # Do some fancy regex string formatting to make the x-axis labels nice
    plt.xticks(
        ticks=locs,
        labels=[re.sub("(.*?_Sub_)", "Subgenome ", x.get_text()) for x in labels],
    )
    plt.yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    plt.xlabel(None)

    N = mpatches.Patch(
        color=None,
        alpha=0,
        label="""TE type: %s \nWindow: %s \nDirection: %s \nNo. Min Tissues: %s \nNo. Gene Pairs: %s"""
        % (
            order.replace("_", " "),  # underscore sub to be pretty
            str(window),
            direction,
            str(n_tissues),
            str(number_genes),
        ),
    )

    # Is this wanted?
    plt.legend(handles=[N], loc="upper right")
    plt.rc("legend", fontsize="small")

    plt.title(
        "TE Density profiles of "
        + long_genome_name
        + " orthologs that are biased towards subgenome "
        + args.bias_str.replace("Sub_", "")
    )
    plt.tight_layout()
    save_file_path = os.path.join(
        args.output_dir,
        "Boxplot",
        args.bias_str,
        (
            args.genome_name
            + "_Ntissues"
            + str(args.num_tissues)
            + "_"
            + order
            + "_"
            + str(window)
            + "_"
            + direction
            + "_"
            + args.bias_str
            + "_DensityBoxPlot.png"
        ),
    )
    os.makedirs(os.path.join(args.output_dir, "Boxplot", args.bias_str), exist_ok=True)
    # plt.show()
    plt.savefig(save_file_path, bbox_inches="tight")
    plt.close()
    return None


def fdr_correct_p_vals(p_val_dict):
    """
    FDR correct the pvalues for the box plot comparisons.

    Args:
        p_val_dict (collections.OrderedDict): Dictionary with keys as TE
        window, direction i.e 'LTR_1000_Upstream. Values are p-val of ttest_ind
        between A and B subgenomes

    Returns:
        boolean_array (numpy.ndarray): True False array that is
        len(p_val_dict). True values are indices where the p_val is less than
        or equal to alpha

        p_val_dict_corrected (collections.OrderedDict): Dictionary with keys as TE
        window, direction i.e 'LTR_1000_Upstream. Values are p-val of ttest_ind
        between A and B subgenomes. Value are FDR corrected.
    """
    names = list(p_val_dict.keys())
    values = list(p_val_dict.values())
    boolean_array, corrected_p_val_array = statsmodels.stats.multitest.fdrcorrection(
        pvals=list(values), alpha=0.05
    )

    p_val_dict_corrected = OrderedDict()  # initialize empty dict
    for name, val in zip(names, corrected_p_val_array):
        p_val_dict_corrected[name] = val

    return boolean_array, p_val_dict_corrected


def add_te_vals_to_syntelogs(syntelogs, te_vals, genome_name):
    """
    Add TE information to the syntelog table, unique because of the format of
    the syntelog table

    Args:
        syntelogs

        te_vals

        genome_name (
    """
    syntelogs_w_te_vals = syntelogs.merge(
        te_vals,
        how="inner",
        left_on=genome_name + "_A",
        right_on="Gene_Name",
        suffixes=["_A", "_B"],
    ).drop(columns=["Gene_Name"])

    # Do again to get the B column
    syntelogs_w_te_vals = syntelogs_w_te_vals.merge(
        te_vals,
        how="inner",
        left_on=genome_name + "_B",
        right_on="Gene_Name",
        suffixes=["_A", "_B"],
    ).drop(columns=["Gene_Name"])
    return syntelogs_w_te_vals


def gen_scatter(
    cleaned_with_te_vals,
    bias_str,
    other,
    order,
    window,
    direction,
    te_column_str,
    num_tissues,
    long_genome_name,
    number_unique_subgenome_tissue_specific_genes,
    output_dir,
):
    # NOTE this code is probably mildly broken because I never moved it out of
    # the exploratory analysis phase because the goalpoasts of my analysis were
    # shifted
    raise NotImplementedError
    plt.figure(figsize=(8, 8))
    sns.scatterplot(
        data=cleaned_with_te_vals,
        y=str(te_column_str + "_" + other),
        x=str(te_column_str + "_" + args.bias_str),
        # line_kws={"color": "black"},
    )
    # NOTE be careful here
    # NOTE this order must correspond to the values given above
    # String formatting to make the x and y axis read nicely.
    plt.xlabel(
        "TE Density Values (Subgenome " + args.bias_str.replace("Sub_", "") + ")"
    )
    plt.ylabel("TE Density Values (Subgenome " + other.replace("Sub_", "") + ")")
    N = mpatches.Patch(
        label="""Total Plotted Genes: %s \nTE type: %s \nWindow: %s \nDirection %s \nNo. Minimum Tissues: %s"""
        % (
            number_unique_subgenome_tissue_specific_genes,
            order,
            window,
            direction,
            str(args.num_tissues),
        )
    )
    plt.legend(handles=[N], loc="upper left")
    plt.ylim(0.0, 1.0)
    plt.xlim(0.0, 1.0)
    plt.title(
        "TE Density profiles of "
        + long_genome_name
        + " orthologs that are biased towards subgenome "
        + args.bias_str.replace("Sub_", "")
    )

    # Remove this if we don't want the 1:1 diagonal line
    diagonal = [0.0, 1.0]
    plt.plot(diagonal, diagonal, color="gray")
    # MAGIC 'Scatter' filepath
    save_file_path = os.path.join(
        args.output_dir,
        "Scatter",
        (
            genome_name
            + "_Ntissues_"
            + str(num_tissues)
            + "_"
            + order
            + "_"
            + str(window)
            + "_"
            + direction
            + "_"
            + bias_str
            + "_DensityScatterPlot.png"
        ),
    )
    plt.savefig(save_file_path)
    # plt.show()
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
        "special_syntelog_file",
        type=str,
        help="""special syntelog file provided
                by collaborators, this is a subset
                and determined by an expression cutoff
                that was not made by me""",
    )

    parser.add_argument(
        "num_tissues",
        type=int,
        help="""Number of tissues that
                you want to subset the HEB syntelog (special syntelog
                file) with. This is the minimum number of tissues a gene must
                be biased in in order for it to be included in
                downstream analyses.""",
    )
    parser.add_argument(
        "bias_str",
        type=str,
        help="""The string
                representing the subgenome genome you want to analyze""",
        choices=["Sub_A", "Sub_B"],
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
        "output_dir",
        type=str,
        help="Parent directory to output results",
    )
    logger = logging.getLogger(__name__)
    args = parser.parse_args()
    args.syntelog_file = os.path.abspath(args.syntelog_file)
    args.special_syntelog_file = os.path.abspath(args.special_syntelog_file)
    args.gene_data = os.path.abspath(args.gene_data)
    args.HDF5_dir = os.path.abspath(args.HDF5_dir)

    # NOTE this is kind of ugly, but I need the B subgenome if I supplied A as
    # an input arg, and vice versa, so:
    if args.bias_str == "Sub_A":
        other = "Sub_B"
    elif args.bias_str == "Sub_B":
        other = "Sub_A"
    else:
        raise ValueError("Need correct string")

    # -----------------------------------------------------------------------
    # Read and initialize various data

    # Read the universal syntelog table
    syntelogs = read_syntelog_file(args.syntelog_file, args.genome_name)

    # Read the special syntelogs table
    special_syntelogs, tissue_list = read_special_gene_set(
        args.special_syntelog_file, args.genome_name
    )

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

    # MAGIC, just get the TE names, and windows for one density data,
    # they correspond for all other density datas too
    # This is used for both the special gene set and the regular gene set.
    all_windows = processed_dd_data[0].window_list  # MAGIC
    directions = ["Upstream", "Downstream"]
    orders = processed_dd_data[0].order_list
    # supers = processed_dd_data[0].super_list  # NB unused, only going to
    # focus on the TE orders

    # MAGIC declare the full genome name so I have something nice to label my
    # figure axes with
    if args.genome_name == "LC":
        long_genome_name = "Luciocarpus capito"
    if args.genome_name == "PR":
        long_genome_name = "Procypris_rabaudi"
    if args.genome_name == "SS":
        long_genome_name = "Spinibarbus_sinensis"

    # -----------------------------------------------------------------------
    # Begin analysis of special genes

    # Transform the special table, group by the genes and make sure they appear in
    # at LEAST the tissue num, so a num_tissues of 3 would mean that Gene X must
    # be observed in at least 3 unique tissues (there is a total of 6)
    special_genes = special_syntelogs[
        special_syntelogs.groupby(args.bias_str + "_Gene")["Tissue"].transform(
            "nunique"
        )
        >= args.num_tissues
    ]

    # Make sure the bias is consistent across all tissues, this, coupled with
    # the command above, ensures that Gene X is observed in AT LEAST
    # args.tissu_num tissues and is biased towards Subgenome A or B (depending
    # on what we are currently testing for)
    special_genes = special_genes.loc[
        special_genes["Bias"] == args.bias_str + "_Biased"
    ]

    # Subset the gene frame with indices with my special set
    shortened_gene_frame_with_indices = gene_frame_with_indices.drop(
        columns=["Start", "Strand", "Stop", "Feature", "Length"]
    )

    # Add the TE HDF5 indices for the main genome to the special gene set
    # Do one subgenome
    # MAGIC gene name column from the GeneData
    special_genes = special_genes.merge(
        shortened_gene_frame_with_indices,
        how="inner",
        left_on=args.bias_str + "_Gene",
        right_on="Gene_Name",
    )

    # Do the other subgenome
    # MAGIC gene name column from the GeneData
    special_genes = special_genes.merge(
        shortened_gene_frame_with_indices,
        how="inner",
        left_on=other + "_Gene",
        right_on="Gene_Name",
        suffixes=("_" + args.bias_str, "_" + other),
    )

    # Drop some columns to remove redundant information that was collected
    # during the above merges
    # MAGIC gene name column from the GeneData
    special_genes.drop(
        columns=["Gene_Name_" + args.bias_str, "Gene_Name_" + other], inplace=True
    )

    # -----------------------------------------------------------------------
    # Begin graphing analysis of special genes
    # TODO should this be put in a function?

    # NOTE, start to loop over windows, upstream/downstream, and TE orders and
    # use the TE Density HDF5 indices to grab the appropriate TE values.
    # Further in the loop we then graph
    windows_i_want = [1000, 2500, 5000, 10000]  # MAGIC these are the windows we are
    # using
    orders_i_want = ["Total_TE_Density", "LTR", "TIR"]
    p_val_dict = OrderedDict()
    for window in windows_i_want:
        for direction in directions:
            for order in orders:
                if order not in orders_i_want:
                    # This is just to avoid generating graphs for the TE types
                    # we don't really care about for this project.
                    continue

                # NOTE this is the string ID of the column that was added to the above pandas
                # data frame, it represents the TE type, window, and direction
                # of values
                te_column_str = "_".join([order, str(window), direction])

                # Add TE vals in now
                cleaned_with_te_vals = add_te_vals_to_gene_info_pandas_from_list_hdf5(
                    special_genes,
                    processed_dd_data,
                    "Order",
                    order,
                    direction,
                    window,
                    args.bias_str + "_Gene",
                    chrom_col="Chromosome_" + args.bias_str,
                    index_col="Index_Val_" + args.bias_str,
                )

                # Rename the TE val column for the previously added
                # args.bias_str subgenome (A or B). This way we can keep track
                # of whether the values are the A or B values
                cleaned_with_te_vals.rename(
                    columns={te_column_str: te_column_str + "_" + args.bias_str},
                    inplace=True,
                )

                # TODO, is there a better way to do this so that I am not repeating code
                # here? This is kind of ugly and bad.
                # TODO I also don't like that I am using this shitty wrapper I
                # wrote.
                # Add in the TE values for other, and chain with a merge to get the total
                # dataframe
                cleaned_with_te_vals = add_te_vals_to_gene_info_pandas_from_list_hdf5(
                    special_genes,
                    processed_dd_data,
                    "Order",
                    order,
                    direction,
                    window,
                    other + "_Gene",
                    chrom_col="Chromosome_" + other,
                    index_col="Index_Val_" + other,
                ).merge(cleaned_with_te_vals)

                # NOTE repeated somewhat, this time for the 'other' subgenome
                # Rename the TE val column for the previously added
                # args.bias_str subgenome (A or B). This way we can keep track
                # of whether the values are the A or B values
                cleaned_with_te_vals.rename(
                    columns={te_column_str: te_column_str + "_" + other},
                    inplace=True,
                )

                ######################################################
                t_val, p_val = scipy.stats.ttest_ind(
                    cleaned_with_te_vals[te_column_str + "_" + args.bias_str],
                    cleaned_with_te_vals[te_column_str + "_" + other],
                )
                p_val_dict[te_column_str] = p_val
                ######################################################

    # FDR correct the pvalues
    # The dictionary preserves order when using the Ordered Dict
    boolean_sig_array, p_val_dict_corrected = fdr_correct_p_vals(p_val_dict)

    # NOTE, RESTART the loop to actually make the graphs this time
    for window in windows_i_want:
        for direction in directions:
            for order in orders:
                if order not in orders_i_want:
                    # This is just to avoid generating graphs for the TE types
                    # we don't really care about for this project.
                    continue

                # NOTE this is the string ID of the column that was added to the above pandas
                # data frame, it represents the TE type, window, and direction
                # of values
                te_column_str = "_".join([order, str(window), direction])

                # Add TE vals in now
                cleaned_with_te_vals = add_te_vals_to_gene_info_pandas_from_list_hdf5(
                    special_genes,
                    processed_dd_data,
                    "Order",
                    order,
                    direction,
                    window,
                    args.bias_str + "_Gene",
                    chrom_col="Chromosome_" + args.bias_str,
                    index_col="Index_Val_" + args.bias_str,
                )

                # Rename the TE val column for the previously added
                # args.bias_str subgenome (A or B). This way we can keep track
                # of whether the values are the A or B values
                cleaned_with_te_vals.rename(
                    columns={te_column_str: te_column_str + "_" + args.bias_str},
                    inplace=True,
                )

                # TODO, is there a better way to do this so that I am not repeating code
                # here? This is kind of ugly and bad.
                # TODO I also don't like that I am using this shitty wrapper I
                # wrote.
                # Add in the TE values for other, and chain with a merge to get the total
                # dataframe
                cleaned_with_te_vals = add_te_vals_to_gene_info_pandas_from_list_hdf5(
                    special_genes,
                    processed_dd_data,
                    "Order",
                    order,
                    direction,
                    window,
                    other + "_Gene",
                    chrom_col="Chromosome_" + other,
                    index_col="Index_Val_" + other,
                ).merge(cleaned_with_te_vals)

                # NOTE repeated somewhat, this time for the 'other' subgenome
                # Rename the TE val column for the previously added
                # args.bias_str subgenome (A or B). This way we can keep track
                # of whether the values are the A or B values
                cleaned_with_te_vals.rename(
                    columns={te_column_str: te_column_str + "_" + other},
                    inplace=True,
                )

                # This is the number of unique gene pairs in the set,
                # Useful for graphing purposes
                number_unique_subgenome_tissue_specific_genes = len(
                    cleaned_with_te_vals[args.bias_str + "_Gene"].unique()
                )

                gen_boxplot(
                    cleaned_with_te_vals,
                    args.bias_str,
                    other,
                    order,
                    window,
                    direction,
                    te_column_str,
                    args.num_tissues,
                    long_genome_name,
                    p_val_dict_corrected[te_column_str],
                    number_unique_subgenome_tissue_specific_genes,
                    args.output_dir,
                )

    # End analysis of special genes
    # -----------------------------------------------------------------------
    # Begin analysis of regular genes
    # TODO I should refactor this out because it only needs to be done once,
    # not twice, once for the A set and once for the B set

    # Do the loops to create graphs for the general set
    for window in windows_i_want:
        for direction in directions:
            for order in orders_i_want:

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

                # TODO edit this table to have a better docstring
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
