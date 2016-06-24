#!/bin/bash
## Define paths to software and reference files

RCODE=/home/gosia/R/shimizu_rna_seq/drought_analysis_2016
RWD=/home/Shared/data/seq/shimizu_rna_seq
ANALYSIS_DIR=Analysis_2016-06-22
DATA_DIR=Data_2016
ROUT=$OUT/Rout

mkdir -p $ROUT

##############################################################################
# Preparation of variables describing samples (new_samps_interpolation.xls, trees.xls)
##############################################################################


R33 CMD BATCH --no-save --no-restore "--args rwd='$RWD' rcode='$RCODE' data_dir='$DATA_DIR' analysis_dir='$ANALYSIS_DIR' out_dir='Plots_Samples' Interpolate_and_roll_function='Interpolate_and_roll.R'" $RCODE/1_PrepareSamples.R $ROUT/1_PrepareSamples.Rout

tail $ROUT/1_PrepareSamples.Rout

##############################################################################
# Prepare description of genes (genes_description_full.xls)
##############################################################################


R33 CMD BATCH --no-save --no-restore "--args rwd='$RWD' rcode='$RCODE' data_dir='$DATA_DIR'" $RCODE/1_PrepareGeneDescriptions.R $ROUT/1_PrepareGeneDescriptions.Rout

tail $ROUT/1_PrepareGeneDescriptions.Rout


##############################################################################
# Crude analysis - MDS plots
##############################################################################


R33 CMD BATCH --no-save --no-restore "--args rwd='$RWD' rcode='$RCODE' data_dir='$DATA_DIR' analysis_dir='$ANALYSIS_DIR' out_dir='Plots_MDS'" $RCODE/2_MDS.R $ROUT/2_MDS.Rout


##############################################################################
# Plots of expression for flowering genes
##############################################################################


R33 CMD BATCH --no-save --no-restore "--args rwd='$RWD' rcode='$RCODE' data_dir='$DATA_DIR' analysis_dir='$ANALYSIS_DIR' out_dir='Plots_of_flowering_genes'" $RCODE/2_ExprPlots.R $ROUT/2_ExprPlots.Rout


##############################################################################
# Clustering
##############################################################################







##############################################################################
# + preparation of variables describing samples
# + crude analysis - MDS plots
# + plots of expression for flowering genes
# + PLOTS OF SELECTED GENES
# + GLM - models with one covariate
# + clustering
# + GO analysis with topGO
# + Fisher's exact tests
# + RUV analysis 
# + Fisher exact test for MolEcol clusters: Sh.l <-> Sh.b 
# + CAMERA Gene set testing
##############################################################################






