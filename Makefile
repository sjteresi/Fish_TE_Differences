# scripts for reproducing analysis
# __file__ Makefile
# __author__ Scott Teresi
#
ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DEV_DATA := $(ROOT_DIR)/data
DEV_RESULTS := $(ROOT_DIR)/results

# NB this only needs to be done once, filepaths are specific to my machine,
# make a virtualenv for the python packages and install the right ones.
#
# Made a virtualenv called Fish_TEs

clean_gene_annotations:
	# First, the Luciobarbus genome 
	python $(ROOT_DIR)/src/clean_fish_gene_annotation.py $(DEV_DATA)/Luciobarbus_capito/Luciobarbus_capito.v20220607.gff3 LC $(DEV_RESULTS)
	# Next, the Procpyris genome
	python $(ROOT_DIR)/src/clean_fish_gene_annotation.py $(DEV_DATA)/Procypris_rabaudi/Procypris_rabaudi.v20220607.gff3 PR $(DEV_RESULTS)
	# Finally, the Spinibarbus genome
	python $(ROOT_DIR)/src/clean_fish_gene_annotation.py $(DEV_DATA)/Spinibarbus_sinensis/Spinibarbus_sinensis.v20220607.gff3 SS $(DEV_RESULTS)

clean_TE_annotations:
	# First, the Luciobarbus genome 
	python $(ROOT_DIR)/src/clean_fish_TE_annotation.py $(DEV_DATA)/Luciobarbus_capito/Luciobarbus_capito.TE.gff LC $(DEV_RESULTS)
	python $(ROOT_DIR)/src/clean_fish_TE_annotation.py $(DEV_DATA)/Procypris_rabaudi/Procypris_rabaudi.TE.gff PR $(DEV_RESULTS)
	python $(ROOT_DIR)/src/clean_fish_TE_annotation.py $(DEV_DATA)/Spinibarbus_sinensis/Spinibarbus_sinensis.TE.gff SS $(DEV_RESULTS)

clean : clean_gene_annotations clean_TE_annotations


analyze_special_LC_A:
	python $(ROOT_DIR)/src/analyze_special_syntelog_density.py \
		$(DEV_DATA)/Kevin_Special_Gene_Set/Diploid_Polyploid_Syntelogs_HEB_For_Scotty.csv \
		3 \
		'Sub_A' \
		LC \
		$(DEV_RESULTS)/Cleaned_LC_Genes.tsv \
		$(DEV_RESULTS)/LC/ \
		$(ROOT_DIR)/config/config.ini \
		$(DEV_RESULTS)/LC/graphs/

analyze_special_LC_B:
	python $(ROOT_DIR)/src/analyze_special_syntelog_density.py \
		$(DEV_DATA)/Kevin_Special_Gene_Set/Diploid_Polyploid_Syntelogs_HEB_For_Scotty.csv \
		3 \
		'Sub_B' \
		LC \
		$(DEV_RESULTS)/Cleaned_LC_Genes.tsv \
		$(DEV_RESULTS)/LC/ \
		$(ROOT_DIR)/config/config.ini \
		$(DEV_RESULTS)/LC/graphs/


analyze_special_PR_A:
	python $(ROOT_DIR)/src/analyze_special_syntelog_density.py \
		$(DEV_DATA)/Kevin_Special_Gene_Set/Diploid_Polyploid_Syntelogs_HEB_For_Scotty.csv \
		3 \
		'Sub_A' \
		PR \
		$(DEV_RESULTS)/Cleaned_PR_Genes.tsv \
		$(DEV_RESULTS)/PR/ \
		$(ROOT_DIR)/config/config.ini \
		$(DEV_RESULTS)/PR/graphs/

analyze_special_PR_B:
	python $(ROOT_DIR)/src/analyze_special_syntelog_density.py \
		$(DEV_DATA)/Kevin_Special_Gene_Set/Diploid_Polyploid_Syntelogs_HEB_For_Scotty.csv \
		3 \
		'Sub_B' \
		PR \
		$(DEV_RESULTS)/Cleaned_PR_Genes.tsv \
		$(DEV_RESULTS)/PR/ \
		$(ROOT_DIR)/config/config.ini \
		$(DEV_RESULTS)/PR/graphs/

analyze_special_SS_A:
	python $(ROOT_DIR)/src/analyze_special_syntelog_density.py \
		$(DEV_DATA)/Kevin_Special_Gene_Set/Diploid_Polyploid_Syntelogs_HEB_For_Scotty.csv \
		3 \
		'Sub_A' \
		SS \
		$(DEV_RESULTS)/Cleaned_SS_Genes.tsv \
		$(DEV_RESULTS)/SS/ \
		$(ROOT_DIR)/config/config.ini \
		$(DEV_RESULTS)/SS/graphs/

analyze_special_SS_B:
	python $(ROOT_DIR)/src/analyze_special_syntelog_density.py \
		$(DEV_DATA)/Kevin_Special_Gene_Set/Diploid_Polyploid_Syntelogs_HEB_For_Scotty.csv \
		3 \
		'Sub_B' \
		SS \
		$(DEV_RESULTS)/Cleaned_SS_Genes.tsv \
		$(DEV_RESULTS)/SS/ \
		$(ROOT_DIR)/config/config.ini \
		$(DEV_RESULTS)/SS/graphs/

analyze_all_special: analyze_special_LC_A analyze_special_LC_B analyze_special_PR_A analyze_special_PR_B analyze_special_SS_A analyze_special_SS_B

regular:
	# LC for now
	python $(ROOT_DIR)/src/analyze_regular_syntelog_density.py \
		$(DEV_DATA)/syntelogs.tsv \
		LC \
		$(DEV_RESULTS)/Cleaned_LC_Genes.tsv \
		$(DEV_RESULTS)/LC/ \
		$(ROOT_DIR)/config/config.ini \
		$(DEV_RESULTS)/LC/graphs/
	# PR
	python $(ROOT_DIR)/src/analyze_regular_syntelog_density.py \
		$(DEV_DATA)/syntelogs.tsv \
		PR \
		$(DEV_RESULTS)/Cleaned_PR_Genes.tsv \
		$(DEV_RESULTS)/PR/ \
		$(ROOT_DIR)/config/config.ini \
		$(DEV_RESULTS)/PR/graphs/
	# SS
	python $(ROOT_DIR)/src/analyze_regular_syntelog_density.py \
		$(DEV_DATA)/syntelogs.tsv \
		SS \
		$(DEV_RESULTS)/Cleaned_SS_Genes.tsv \
		$(DEV_RESULTS)/SS/ \
		$(ROOT_DIR)/config/config.ini \
		$(DEV_RESULTS)/SS/graphs/
