"""
Filter the TE annotations to the appropriate format for TE Density algorithm
"""

__author__ = "Scott Teresi"

import pandas as pd
import argparse
import os
import logging
import coloredlogs


def check_nulls(my_df, logger):
    """Check the TE dataframe for ANY null values in ANY rows

    Args:
        my_df (pandas.core.DataFrame): Pandas dataframe of TE values from TE
            annotation
    """
    Bool = my_df.isnull().values.any()
    if Bool:
        logger.critical("You have null values in your dataframe!")
        logger.critical("Here are the null values in the output:")
        null_columns = my_df.columns[my_df.isnull().any()]
        print((my_df[my_df.isnull().any(axis=1)][null_columns].head()))


def write_cleaned_transposons(te_pandaframe, output_dir, genome_name, logger):
    file_name = os.path.join(output_dir, ("Cleaned_" + genome_name + "_EDTA_TEs.tsv"))

    logger.info("Writing cleaned TE file to: %s" % file_name)
    te_pandaframe.to_csv(file_name, sep="\t", header=True, index=False)


def import_transposons(tes_input_path, te_annot_renamer, genome_name, logger):
    """Import TE file and read as a dataframe in Pandas

    Args:
        tes_input_path (str): string of the file path to the TE annotation

        logger (logging obj): The object to call logs and info
    """
    col_names = [
        "Chromosome",
        "Software",
        "Feature",
        "Start",
        "Stop",
        "Score",
        "Strand",
        "Phase",
        "Attribute",
    ]

    TE_Data = pd.read_csv(
        tes_input_path,
        sep="\t+",
        header=None,
        engine="python",
        names=col_names,
        comment="#",
        dtype={
            "Start": "float64",
            "Stop": "float64",
            "Chromosome": str,
            "Strand": str,
            "Attribute": str,
        },
    )

    # Drop extraneous columns
    TE_Data.drop(columns=["Score", "Software", "Phase", "Feature"], inplace=True)

    # Create Order and SuperFamily column from Attribute column
    # Because that column contains the detailed TE information
    # Then remove old Attribute column
    TE_Data["Attribute"] = TE_Data["Attribute"].str.extract(r"Classification=(.*?);")

    TE_Data[["Order", "SuperFamily"]] = TE_Data.Attribute.str.split("/", expand=True)
    TE_Data.drop(columns=["Attribute"], inplace=True)
    TE_Data.Order = TE_Data.Order.astype(str)
    TE_Data.SuperFamily = TE_Data.SuperFamily.astype(str)

    # Call renamer
    TE_Data = te_annot_renamer(TE_Data)

    # Declare data types
    TE_Data["Length"] = TE_Data.Stop - TE_Data.Start + 1
    check_nulls(TE_Data, logger)

    TE_Data.sort_values(by=["Chromosome", "Start"], inplace=True)

    # MAGIC I only want the first 7
    # MAGIC chromosome/scaffold names
    # NOTE this is kind of bad form to be doing things this way, alternatively
    # could provide some sort of input arg to handle the chromosomes I want but
    # this is easiest to hard-code it in.
    # NOTE chromosomes_i_want is a string variable that is set per genome and
    # used to subset the dataframe by chromosome ID
    if genome_name == "562":
        # NOTE prepend a chromosome name to the integers of chromosome IDs
        TE_Data["Chromosome"] = "562_scaffold_" + TE_Data["Chromosome"]
        chromosomes_i_want = ["562_scaffold_" + str(i) for i in range(8)]

    if genome_name == "2339":
        # NOTE prepend a chromosome name to the integers of chromosome IDs
        TE_Data["Chromosome"] = "2339_scaffold_" + TE_Data["Chromosome"]
        chromosomes_i_want = ["2339_scaffold_" + str(i) for i in range(8)]

    if genome_name == "H4":
        # NOTE chromosome IDs already in a good format, just define range
        chromosomes_i_want = ["Fvb" + str(i) for i in range(8)]

    if genome_name == "502":
        # NOTE prepend a chromosome name to the integers of chromosome IDs
        TE_Data["Chromosome"] = "502_scaffold_" + TE_Data["Chromosome"]
        chromosomes_i_want = ["502_scaffold_" + str(i) for i in range(8)]

    TE_Data = TE_Data.loc[TE_Data["Chromosome"].isin(chromosomes_i_want)]
    diagnostic(TE_Data)
    return TE_Data


def diagnostic(TE_Data):
    print(sorted(TE_Data["Order"].unique().tolist()))
    print(sorted(TE_Data["SuperFamily"].unique().tolist()))


def te_annot_renamer(TE_Data):
    U = "Unknown_Order"
    master_order = {
        "Unknown": U,
        "MITE": "TIR",
        "pararetrovirus": "pararetrovirus",
        "DNA": "TIR",
    }

    U = "Unknown_Superfam"
    master_superfamily = {
        # EDTA/Wicker et al 2007 renames to common name:
        "RLC": "Copia",
        "RLG": "Gypsy",
        "RLB": "Bel_Pao",
        "RLR": "Retrovirus",
        "RLE": "ERV",
        "RYD": "DIRS",
        "RYN": "Ngaro",
        "RYV": "VIPER",
        "RPP": "Penelope",
        "RIR": "R2",
        "RIT": "RTE",
        "RIJ": "Jockey",
        "RIL": "L1",
        "RII": "I",
        "RST": "tRNA",
        "RSL": "7SL",
        "RSS": "5S",
        "DTT": "Tc1-Mariner",
        "DTA": "hAT",
        "DTM": "Mutator",
        "DTE": "Merlin",
        "DTR": "Transib",
        "DTP": "P",
        "DTB": "PiggyBac",
        "DTH": "PIF-Harbinger",
        "DTC": "CACTA",
        "DYC": "Crypton",
        "DHH": "Helitron",
        "DMM": "Maverick",
        # Custom changes
        "unknown": U,
        "Unknown": U,
        "None": U,
        "EnSpm_CACTA": "CACTA",
        "MuDR_Mutator": "Mutator",
        "PIF_Harbinger": "PIF-Harbinger",
    }

    TE_Data.SuperFamily.fillna(
        value="Unknown_Superfam", inplace=True
    )  # replace None w U

    # Invoke dictionary to fix names
    TE_Data.Order.replace(master_order, inplace=True)
    TE_Data.SuperFamily.replace(master_superfamily, inplace=True)

    # Rename the superfamily value for pararetros as pararetrovirus
    TE_Data.loc[TE_Data.Order == "pararetrovirus", "SuperFamily"] = "pararetrovirus"

    # Rename unknown LINE element superfamilies to Unknown_LINE_Superfam to
    # distinguish between other unknowns
    TE_Data.loc[
        (TE_Data.Order == "LINE") & (TE_Data["SuperFamily"] == "Unknown_Superfam"),
        "SuperFamily",
    ] = "Unknown_LINE_Superfam"

    # Remove 'Maverick' TEs from plant genome
    TE_Data = TE_Data.drop(TE_Data[TE_Data.Order == "Maverick"].index)

    # Rename unknown LTR element superfamilies to Unknown_LTR_Superfam to
    # distinguish between other unknowns
    TE_Data.loc[
        (TE_Data["Order"] == "LTR") & (TE_Data["SuperFamily"] == "Unknown_Superfam"),
        "SuperFamily",
    ] = "Unknown_LTR_Superfam"

    # Rename unknown TIR element superfamilies to Unknown_TIR_Superfam to
    # distinguish between other unknowns
    TE_Data.loc[
        (TE_Data.Order == "TIR") & (TE_Data["SuperFamily"] == "Unknown_Superfam"),
        "SuperFamily",
    ] = "Unknown_TIR_Superfam"

    # Rename both values for Helitron elements, so that 'Helitron' is
    # both the Order and SuperFamily value
    # Some Helitron elements were labeled 'DNA' in the Order location, this is
    # technically correct but I prefer to differentiate the TIR DNA elements
    # from DNA elements as a whole
    TE_Data.loc[
        (TE_Data["Order"] == "TIR") & (TE_Data["SuperFamily"] == "Helitron"),
        ["Order", "SuperFamily"],
    ] = "Helitron"
    # If the Order is Helitron and the SuperFamily is unknown make the
    # superfamily 'Helitron'
    TE_Data.loc[
        (TE_Data["Order"] == "Helitron")
        & (TE_Data["SuperFamily"] == "Unknown_Superfam"),
        "SuperFamily",
    ] = "Helitron"

    # For TEs that are unknown for both Order AND SuperFamily we will call
    # those 'Completely_Unknown'
    TE_Data.loc[
        (TE_Data["Order"] == "Unknown_Order")
        & (TE_Data["SuperFamily"] == "Unknown_Superfam"),
        ["Order", "SuperFamily"],
    ] = "Completely_Unknown"

    return TE_Data


def diagnostic_cleaner_helper(TE_Data):
    print()
    print(TE_Data.Order.unique())
    print(TE_Data.SuperFamily.unique())
    print()

    # To see unique for a given type:
    # print(TE_Data.loc[TE_Data['Order'] == 'LINE'].SuperFamily.unique())
    return None


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Reformat TE annotation file")
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    output_default = os.path.join(dir_main, "../results")
    parser.add_argument(
        "TE_input_file", type=str, help="Parent path of TE annotation file"
    )
    parser.add_argument("genome_name", type=str, help="Name for the genome")

    parser.add_argument(
        "--output_dir",
        "-o",
        type=str,
        default=output_default,
        help="Parent directory to output results",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.TE_input_file = os.path.abspath(args.TE_input_file)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Execute
    cleaned_transposons = import_transposons(
        args.TE_input_file, te_annot_renamer, args.genome_name, logger
    )
    write_cleaned_transposons(
        cleaned_transposons, args.output_dir, args.genome_name, logger
    )
