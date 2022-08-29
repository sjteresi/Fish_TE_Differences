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

analyze_LC_syntelogs:
	mkdir -p $(DEV_RESULTS)/LC/graphs
	python $(ROOT_DIR)/src/analyze_syntelog_density.py \
		$(DEV_DATA)/syntelogs.tsv \
		LC \
		$(DEV_RESULTS)/Cleaned_LC_Genes.tsv \
		$(DEV_RESULTS)/LC/ \
		$(DEV_RESULTS)/LC/graphs/

analyze_PR_syntelogs:
	mkdir -p $(DEV_RESULTS)/PR/graphs
	python $(ROOT_DIR)/src/analyze_syntelog_density.py \
		$(DEV_DATA)/syntelogs.tsv \
		PR \
		$(DEV_RESULTS)/Cleaned_PR_Genes.tsv \
		$(DEV_RESULTS)/PR/ \
		$(DEV_RESULTS)/PR/graphs/

analyze_SS_syntelogs:
	mkdir -p $(DEV_RESULTS)/SS/graphs
	python $(ROOT_DIR)/src/analyze_syntelog_density.py \
		$(DEV_DATA)/syntelogs.tsv \
		SS \
		$(DEV_RESULTS)/Cleaned_SS_Genes.tsv \
		$(DEV_RESULTS)/SS/ \
		$(DEV_RESULTS)/SS/graphs/

