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
	# TODO flesh this out for all gene annotations
	# Let's do the Luciobarbus genome first
	python $(ROOT_DIR)/src/clean_fish_gene_annotation.py
