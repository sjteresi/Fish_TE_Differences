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

    # Drop extraneous columns
    col_to_use = [
        "Chromosome",
        "Start",
        "Stop",
        "Attribute",
    ]

    TE_Data = pd.read_csv(
        tes_input_path,
        sep="\t+",
        header=None,
        engine="python",
        names=col_names,
        usecols=col_to_use,
        comment="#",
        dtype={
            "Start": "float64",
            "Stop": "float64",
            "Chromosome": str,
            "Attribute": str,
        },
    )
    # MAGIC
    # This is specifically for these fish genomes. I
    # don't care about the million scaffolds
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

    # Filter on my chromosome whitelist
    TE_Data = TE_Data.loc[TE_Data["Chromosome"].isin(chromosomes_i_want)]

    # Create Order and SuperFamily column from Attribute column
    # Because that column contains the detailed TE information
    # Then remove old Attribute column
    TE_Data["Attribute"] = TE_Data["Attribute"].str.extract(r"Class=(.*?);")
    # NOTE fix early on in the pipeline because the SS genome creator decided
    # to use their delimiter character in a TE name, not present in the
    # other two genomes. Normally this fix would be in the rename section but
    # the str split command below fails if there are more than one '/' in the
    # field.
    TE_Data.loc[TE_Data["Attribute"] == "DNA/EN/SPM", "Attribute"] = "DNA/EnSpm"
    TE_Data[["Order", "SuperFamily"]] = TE_Data.Attribute.str.split("/", expand=True)
    TE_Data.drop(columns=["Attribute"], inplace=True)
    TE_Data.Order = TE_Data.Order.astype(str)
    TE_Data.SuperFamily = TE_Data.SuperFamily.astype(str)

    # NOTE this re-categorizes many TEs
    TE_Data = te_annot_renamer(TE_Data)

    TE_Data["Length"] = TE_Data.Stop - TE_Data.Start + 1
    check_nulls(TE_Data, logger)
    TE_Data.sort_values(by=["Chromosome", "Start"], inplace=True)

    return TE_Data


def te_annot_renamer(TE_Data):
    U = "Unknown_Order"
    master_order = {
        "Unknown": U,
        # Rename the Order value for DNA orders to TIR to better conform to the
        # Wicker naming system
        "DNA": "TIR",
        # Rename the RC Order to Helitron
        "RC": "Helitron",
    }

    U = "Unknown_Superfam"
    master_superfamily = {
        # Custom changes
        "unknown": U,
        "Unknown": U,
        "None": U,
        "U": U,
        # "EnSpm_CACTA": "CACTA",
        # "MuDR_Mutator": "Mutator",
        # "PIF_Harbinger": "PIF-Harbinger",
        # rename all of the hAT subtypes
        "hAT-Ac": "hAT",
        "hAT-Blackjack": "hAT",
        "hAT-Charlie": "hAT",
        "hAT-Tag1": "hAT",
        "hAT-Tip100": "hAT",
        "hAT-hAT5": "hAT",
        "hAT-hobo": "hAT",
        "hAT-hAT6": "hAT",
        # Rename all of the TcMar subtypes
        "TcMar-ISRm11": "Tc1_Mariner",
        "Tc1-Mariner": "Tc1_Mariner",
        "TcMar-Tigger": "Tc1_Mariner",
        "TcMar": "Tc1_Mariner",
        "TcMar-Tc2": "Tc1_Mariner",
        "TcMar-Tc4": "Tc1_Mariner",
        "TcMar-Tc1": "Tc1_Mariner",
        "TcMar-Fot1": "Tc1_Mariner",
        # Rename all of the Crypton subtypes
        "Crypton-A": "Crypton",
        "Crypton-V": "Crypton",
        # Rename all of the Sola subtypes
        "Sola-1": "Sola",
        "Sola-2": "Sola",
        # Rename all of the L1 subtypes
        "L1-Tx1": "L1",
        # Rename all of the SINE subtypes
        "5S-Deu-L2": "5S",
        # Rename all of the Kolobok subtypes
        "Kolobok-T2": "Kolobok",
        # Rename all of the tRNA subtypes
        "tRNA-L2": "tRNA",
        "tRNA-RTE": "tRNA",
        "tRNA-Core": "tRNA",
        "tRNA-Core-RTE": "tRNA",
        # Rename all of the RTE subtypes
        "RTE-X": "RTE",
        "RTE-BovB": "RTE",
        # Rename all of the CMC subtypes
        "CMC-EnSpm": "CMC",
        "EnSpm": "CMC",
        "CMC-Chapaev-3": "CMC",
        "CACTA": "CMC",
        # Rename all of the MULE subtypes
        "MULE-MuDR": "MULE",
        "MULE-NOF": "MULE",
        # Rename all of the Academ subtypes
        "Academ-1": "Academ",
        # Rename all of the R2 subtypes
        "R2-Hero": "R2",
        # Rename all of the ERV subtypes
        "ERV1": "ERV",
        "ERVL": "ERV",
        "ERVK": "ERV",
        # Rename the PIF/Harbinger and IS2/3 elements to a superfamily of their own
        "PIF-ISL2EU": "PHIS",
        "PIF-Harbinger": "PHIS",
        "IS3EU": "PHIS",
        "PIF": "PHIS",
        # Rename the Proto2 elements to be a part of RTE
        # https://www.girinst.org/2009/vol9/issue7/Proto2-1_BF.html
        "Proto2": "RTE",
    }
    # Invoke dictionary to fix names
    TE_Data.Order.replace(master_order, inplace=True)
    TE_Data.SuperFamily.replace(master_superfamily, inplace=True)

    # The dictionary is the easiest way to rename things. But if you need to be
    # more creative or set the value for the Order column based on the
    # superfamily column, you can start using pandas loc notation.

    # If only the superfamily is None (unknown), and we DO know the order,
    # rewrite the superfamily to be unknown
    TE_Data.loc[
        (TE_Data["Order"] != "Unknown_Order") & (TE_Data["SuperFamily"] == "None"),
        ["SuperFamily"],
    ] = "Unknown_Superfam"

    # Rename the Order value for DIRS superfams to DIRS, they are not LTR
    TE_Data.loc[TE_Data.SuperFamily == "Ngaro", "Order"] = "DIRS"
    TE_Data.loc[TE_Data.SuperFamily == "DIRS", "Order"] = "DIRS"

    # Rename the Order value for Penelope superfam to DIRS, they are not LINE
    TE_Data.loc[TE_Data.SuperFamily == "Penelope", "Order"] = "PLE"

    # For TEs that are unknown for both Order AND SuperFamily we will call
    # those 'Completely_Unknown'
    TE_Data.loc[
        (TE_Data["Order"] == "Unknown_Order")
        & (TE_Data["SuperFamily"] == "Unknown_Superfam"),
        ["Order", "SuperFamily"],
    ] = "Completely_Unknown"

    # This section is for help. You can use it to check your groups
    # print("NOTE")
    # print(sorted(TE_Data["Order"].unique().tolist()))
    # print(sorted(TE_Data["SuperFamily"].unique().tolist()))
    # print(len(TE_Data["SuperFamily"].unique().tolist()))
    # print()

    return TE_Data


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Reformat TE annotation file")
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser.add_argument(
        "TE_input_file", type=str, help="Parent path of TE annotation file"
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
    args.TE_input_file = os.path.abspath(args.TE_input_file)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Execute
    cleaned_transposons = import_transposons(
        args.TE_input_file, te_annot_renamer, args.genome_name, logger
    )
    file_name = os.path.join(
        args.output_dir, ("Cleaned_" + args.genome_name + "_TEs.tsv")
    )
    cleaned_transposons.to_csv(file_name, sep="\t", header=True, index=False)
