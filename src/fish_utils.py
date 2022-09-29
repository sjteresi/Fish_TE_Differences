"""
Sundry helper functions
"""
import pandas as pd
import logging as logger
from functools import partial
from configparser import ConfigParser

from transposon.import_filtered_genes import import_filtered_genes
from transposon.density_data import DensityData
from transposon.gene_data import GeneData


def parse_analysis_config(config_path):
    """Return parameters for running fish TE density calculations."""

    parser = ConfigParser(
        converters={"list": lambda x: [i.strip() for i in x.split(",")]}
    )
    parser.read(config_path)
    orders = list(map(str, parser.getlist("analysis_parameters", "orders")))
    windows = list(map(int, parser.getlist("analysis_parameters", "windows")))
    directions = list(map(str, parser.getlist("analysis_parameters", "directions")))
    return windows, directions, orders


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
