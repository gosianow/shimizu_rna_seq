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
# MDS plots based on all genes, housekeeping genes, drought genes
# Coloring by tree, time, wp level
##############################################################################


R33 CMD BATCH --no-save --no-restore "--args rwd='$RWD' rcode='$RCODE' data_dir='$DATA_DIR' analysis_dir='$ANALYSIS_DIR' out_dir='Plots_MDS'" $RCODE/2_MDS.R $ROUT/2_MDS.Rout


##############################################################################
# Plots of expression for interesting gene sets: flowering genes, drought genes
# Line plots for all the samples and heatmaps for 970, 8266, 1099 and 1377 trees
##############################################################################


R33 CMD BATCH --no-save --no-restore "--args rwd='$RWD' rcode='$RCODE' data_dir='$DATA_DIR' analysis_dir='$ANALYSIS_DIR' out_dir='Plots_of_flowering_genes'" $RCODE/2_ExprPlots.R $ROUT/2_ExprPlots.Rout


##############################################################################
# DGE analysis with DESeq2
##############################################################################


de_model="between_trees_per_tree_no_time"

R33 CMD BATCH --no-save --no-restore "--args rwd='$RWD' rcode='$RCODE' data_dir='$DATA_DIR' analysis_dir='$ANALYSIS_DIR' out_dir='Results_DGE' de_model='$de_model'" $RCODE/3_DGE_deseq2.R $ROUT/3_DGE_deseq2_${de_model}.Rout

tail $ROUT/3_DGE_deseq2_${de_model}.Rout


de_model="between_trees_per_tree"

R33 CMD BATCH --no-save --no-restore "--args rwd='$RWD' rcode='$RCODE' data_dir='$DATA_DIR' analysis_dir='$ANALYSIS_DIR' out_dir='Results_DGE' de_model='$de_model'" $RCODE/3_DGE_deseq2.R $ROUT/3_DGE_deseq2_${de_model}.Rout

tail $ROUT/3_DGE_deseq2_${de_model}.Rout


de_model="interaction_any_time_any_tree"

R33 CMD BATCH --no-save --no-restore "--args rwd='$RWD' rcode='$RCODE' data_dir='$DATA_DIR' analysis_dir='$ANALYSIS_DIR' out_dir='Results_DGE' de_model='$de_model'" $RCODE/3_DGE_deseq2.R $ROUT/3_DGE_deseq2_${de_model}.Rout

tail $ROUT/3_DGE_deseq2_${de_model}.Rout


de_model="interaction_per_time_per_tree"

R33 CMD BATCH --no-save --no-restore "--args rwd='$RWD' rcode='$RCODE' data_dir='$DATA_DIR' analysis_dir='$ANALYSIS_DIR' out_dir='Results_DGE' de_model='$de_model'" $RCODE/3_DGE_deseq2.R $ROUT/3_DGE_deseq2_${de_model}.Rout

tail $ROUT/3_DGE_deseq2_${de_model}.Rout


de_model="between_trees_any_tree"

R33 CMD BATCH --no-save --no-restore "--args rwd='$RWD' rcode='$RCODE' data_dir='$DATA_DIR' analysis_dir='$ANALYSIS_DIR' out_dir='Results_DGE' de_model='$de_model'" $RCODE/3_DGE_deseq2.R $ROUT/3_DGE_deseq2_${de_model}.Rout

tail $ROUT/3_DGE_deseq2_${de_model}.Rout


de_model="between_time_all_trees"

R33 CMD BATCH --no-save --no-restore "--args rwd='$RWD' rcode='$RCODE' data_dir='$DATA_DIR' analysis_dir='$ANALYSIS_DIR' out_dir='Results_DGE' de_model='$de_model'" $RCODE/3_DGE_deseq2.R $ROUT/3_DGE_deseq2_${de_model}.Rout

tail $ROUT/3_DGE_deseq2_${de_model}.Rout


de_model="intercept"

R33 CMD BATCH --no-save --no-restore "--args rwd='$RWD' rcode='$RCODE' data_dir='$DATA_DIR' analysis_dir='$ANALYSIS_DIR' out_dir='Results_DGE' de_model='$de_model'" $RCODE/3_DGE_deseq2.R $ROUT/3_DGE_deseq2_${de_model}.Rout

tail $ROUT/3_DGE_deseq2_${de_model}.Rout








##############################################################################
# Clustering
##############################################################################












