# Contact Information:
| Role          | Name          | GitHub                                                  | Email              |
|---------------|---------------|---------------------------------------------------------|--------------------|
| Analyst: | Scott Teresi  | [Personal GitHub](https://github.com/sjteresi) | <teresisc@msu.edu> |
| PI:           | Patrick Edger | [Lab GitHub](https://github.com/EdgerLab)               | <edgerpat@msu.edu> |

# Goal:
There are 3 fish genomes that are polyploids. Each fish genome has an A and B subgenome. Investigate if there are any TE differences relative to genes amongst the A and B subgenomes.

## Aims:
* Generate a histogram of the *difference* in TE Density values of syntelogs
* Generate a boxplot of the TE Density values of syntelogs

# Methods:
TE Density version (TODO, commit: `4b4d7c091733d2fbe2e537fb7486c50816e55037`) (TODO CITE) was run independently with default options for each genome.
TE Density calculates the proportion of TE-occupied base-pairs relative to genes and a given window measurement size.
TE density was calculated for the combination of (TE superfamily identity ∥ TE order identity)×(upstream ∥ intragenic ∥ downstream), with respect to a window length and an individual gene's start and stop positions.
The output matrices, representing the TE density data for each pseudomolecule, are of size |identity|×|windows|×|genes|×|direction|, where identity is the set of either the TE superfamilies or orders, windows is the set of window lengths, and genes is the set individual gene names.
The direction is the relative location of the window to a gene’s start and stop position, where direction∈upstream, intra, downstream.
Gene and TE annotation files, the inputs to the software, were reformatted to conform to the requirements of TE Density.
Analyses were conducted and graphs were generated using Python 3.8.0.
Version-controlled code and documentation related to this analysis can be found at [https://github.com/sjteresi/Fish_TE_Differences](https://github.com/sjteresi/Fish_TE_Differences), see the requirements directory in the project GitHub repository for a more complete list of minor packages.
A two-sample t-test was calculated for the boxplots comparing A and B subgenome TE density values. 
TODO add in something about the tissues?

# Results:

- General Histogram Caption:
	* Histogram of differences in TE density values of syntelogs for Subgenome A and B. The density values for the TE type being shown are derived from the [read fig TE type legend] grouping of TEs, and values were collected in a [read fig window legend] bp window [read fig direction legend] of genes. TE density difference values were calculated by subtracting the B subgenome syntelog values from the A subgenome syntelog vaues. Negative values on the graph reflect a higher TE value for the B subgenome syntelog and positive values reflect a higher TE value for the A subgenome syntelog. Values were binned into groupings reflecting 10% increases or decreases in TE density. All but the last (most positive) bin is half-open. For example, the leftmost bin relfects an interval of [-1.0, -0.9) and the rightmost bin reflects an interval of [0.9, 1.0]. The number of syntelog pairs shown in the histogram is [read fig total plotted genes legend]; the number of syntelog pairs with any difference in TE density values is [read fig no. 0 differences legend]. NOTE, this was essentially copied from the TE Density publication.

-General Boxplot Notes:
	* Boxplot of TE Density values of syntelogs for Subgenome A and B. The figures are a subset of syntelogs that are biased (HEB test done by collaborators) towards one subgenome. We settled on using a 3 tissue cutoff. That means that the syntelog pairs had to have a consistent bias towards one subgenome in **at least** 3 tissues. For example, given the `LC_Ntissues3_Total_TE_Density_10000_Upstream_Sub_A_DensityBoxPlot.svg`, the genes plotted there are syntelogs on the A and B subgenomes, the A syntelog was observed to be the "dominant" syntelog in at least 3 tissues (there were 6 I believe). 

## Figures:
- *Luciobarbus capito*:
	* `LC_Total_TE_Density_10000_Upstream_DensityDifferences.svg`
	* `LC_Total_TE_Density_10000_Downstream_DensityDifferences.svg`
	* `LC_Ntissues3_Total_TE_Density_10000_Upstream_Sub_A_DensityBoxPlot.svg`
	* `LC_Ntissues3_Total_TE_Density_10000_Downstream_Sub_A_DensityBoxPlot.svg`
	* `LC_Ntissues3_Total_TE_Density_10000_Upstream_Sub_B_DensityBoxPlot.svg`
	* `LC_Ntissues3_Total_TE_Density_10000_Downstream_Sub_B_DensityBoxPlot.svg`

- *Spinibarbus sinensis*:
	* `SS_Total_TE_Density_10000_Upstream_DensityDifferences.svg`
	* `SS_Total_TE_Density_10000_Downstream_DensityDifferences.svg`
	* `SS_Ntissues3_Total_TE_Density_10000_Upstream_Sub_A_DensityBoxPlot.svg`
	* `SS_Ntissues3_Total_TE_Density_10000_Downstream_Sub_A_DensityBoxPlot.svg`
	* `SS_Ntissues3_Total_TE_Density_10000_Upstream_Sub_B_DensityBoxPlot.svg`
	* `SS_Ntissues3_Total_TE_Density_10000_Downstream_Sub_B_DensityBoxPlot.svg`

- *Procypris rabaudi*:
	* `PR_Total_TE_Density_10000_Upstream_DensityDifferences.svg`
	* `PR_Total_TE_Density_10000_Downstream_DensityDifferences.svg`
	* `PR_Ntissues3_Total_TE_Density_10000_Upstream_Sub_A_DensityBoxPlot.svg`
	* `PR_Ntissues3_Total_TE_Density_10000_Downstream_Sub_A_DensityBoxPlot.svg`
	* `PR_Ntissues3_Total_TE_Density_10000_Upstream_Sub_B_DensityBoxPlot.svg`
	* `PR_Ntissues3_Total_TE_Density_10000_Downstream_Sub_B_DensityBoxPlot.svg`

# Code:
The code is broken up into several difference scripts inside the `src/` directory. The scripts below can most easily be invoked from the `Makefile`.

## Generation of TE Density Data:
* `clean_fish_gene_annotation.py`: Filters a gene annotation for TE Density tool
	- Inputs: Gene annotation file, provided by collaborators
	- Outputs: `Cleaned_(genome_name)_Genes.tsv`, a reformatted gene annotation
* `clean_fish_TE_annotation.py`: Filters a TE annotation for TE Density tool
	- Inputs: TE annotation file, provided by collaborators
	- Outputs: `Cleaned_(genome_name)_TEs.tsv`, a reformatted TE annotation
* `TE_Density_LC.sb`: Runs TE Density on the LC (*Luciobarbus capito*) genome
* `TE_Density_PR.sb`: Runs TE Density on the PR (*Procypris rabaudi*) genome
* `TE_Density_SS.sb`: Runs TE Density on the SS (*Spinibarbus sinensis*) genome

## Analysis of Syntelogs
* `analyze_regular_syntelog_density`: Analaysis script. Reads all of the TE Density data and the regular syntelog data, combines them, and creates the histograms.
	- Inputs: Regular syntelog file provided by collaborators, and the set of .h5 files that are the output from TE Density.
	- Outputs: Histograms
* `analyze_special_syntelog_density`: Analysis script. Reads all of the TE Density data and the special syntelog data, combines them, and creates the boxplots.
	- Inputs: "Special" syntelog (these are subgenome biased genes) file provided by collaborators, and the set of .h5 files that are the output from TE Density.
	- Outputs: Boxplots

## Helper Code:
* `fish_utils.py`: Contains functions useful to the analysis of both the regular syntelogs and the "special" HEB biased syntelogs
* `sync_local_to_remote.sh`: Utility code for Scott
* `sync_remote_to_local.sh`: Utility code for Scott

## Requirements:
See `requirements/requirements.txt` for a complete list. The TE Density SBATCH scripts contain information required to run TE Density on MSU's HPCC.
