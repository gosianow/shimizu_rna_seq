##############################################################################
## <<1_PrepareGeneDescriptions.R>>

# BioC 3.3
# Created 23 June 2016


##############################################################################
Sys.time()
##############################################################################

library(org.At.tair.db)
library(limma)

##############################################################################
# Test arguments
##############################################################################

rwd='/home/Shared/data/seq/shimizu_rna_seq'
rcode='/home/gosia/R/shimizu_rna_seq/drought_analysis_2016'
data_dir='Data_2016'


##############################################################################
# Read in the arguments
##############################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(args)

##############################################################################

setwd(rwd)


##############################################################################
# Prepare gene descriptions
##############################################################################


### Annotation S. leprosula <-> A. thaliana  gene matches 

AT.id <- read.table(paste0(data_dir, "/Data_2016-06-22/orig/Blastp_result_against_Athaliana.txt"), sep="\t", stringsAsFactors = FALSE, header = FALSE)
head(AT.id)

colnames(AT.id)[1] <- "sl_id"

AT.id$at_id <- strsplit2(AT.id[,2], "\\.")[ ,1]


### Get AT symbols

syms <- mget(AT.id[, "at_id"], org.At.tairSYMBOL)
syms <- sapply(syms, function(u) if(is.na(u[1])) NA else paste(u,collapse=","))
AT.id$at_symbol <- syms


### Get functional description from TAIR10

tfd <- read.table(paste0(data_dir, "/Data_2016-06-22/orig/TAIR10_functional_descriptions.txt"), sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)

m <- match(AT.id[, "at_id"], gsub("\\.[1-9]", "", tfd$Model_name))

AT.id$at_description_tair10 <- as.character(tfd$Short_description[m])

### Remove problematic symbols from description
AT.id$at_description_tair10 <- gsub("[^[:alnum:],()+/\\-]", " ", AT.id$at_description_tair10)
AT.id$at_symbol <- gsub("[^[:alnum:],()+/\\-]", " ", AT.id$at_symbol)


genes.description <- AT.id[, c("sl_id", "at_id", "at_symbol", "at_description_tair10")]


write.table(genes.description, paste0(data_dir, "/Data_2016-06-22/genes_description.xls"), sep="\t", quote = FALSE, row.names = FALSE)


### Flowering genes - Table S4 from MolEcol paper

Athaliana.flowering.genes <- read.table(paste0(data_dir,"/GeneControlSets/Flowering_Genes/gene_list_flowering_related_MolEcol_Athaliana_S4.csv"), sep=",", header=T, stringsAsFactors=FALSE, skip=1)

# Some clean up 
Athaliana.flowering.genes[,1] <- gsub(" ", "", Athaliana.flowering.genes[,1])
colnames(Athaliana.flowering.genes) <- c("at_id", "at_description_flowering")
Athaliana.flowering.genes$flowering <- TRUE

Athaliana.flowering.genes$at_description_flowering <- gsub("[^[:alnum:],()+/\\-]", " ", Athaliana.flowering.genes$at_description_flowering)



### Moderate drought responce genes 

UP.genes <- read.table(paste0(data_dir,"/GeneControlSets/Drought_Genes/mDr_Day10_drought_up_regulated_genes_Harb_etal.csv"), sep=";", header=T, stringsAsFactors=FALSE)

DOWN.genes <- read.table(paste0(data_dir,"/GeneControlSets/Drought_Genes/mDr_Day10_drought_down_regulated_genes_Harb_etal.csv"), sep=";", header=T, stringsAsFactors=FALSE)


Drought.genes <- rbind(UP.genes[, c("Gene", "mDr.Day10", "Function")], DOWN.genes[, c("Gene", "mDr.Day10", "Function")])

names(Drought.genes) <- c("at_id", "drought_regulation", "at_description_drought")

Drought.genes$at_description_drought <- gsub("[^[:alnum:],()+/\\-]", " ", Drought.genes$at_description_drought)



### Merge gene descriptions

genes.full.description <- merge(genes.description, Drought.genes, by="at_id", all.x=TRUE, sort = FALSE)


genes.full.description <- merge(genes.full.description, Athaliana.flowering.genes, by = "at_id", all.x = TRUE, sort = FALSE)

### IMPORTANT to take unique bcs there are duplicates after merge
genes.full.description <- unique(genes.full.description)


write.table(genes.full.description, paste0(data_dir, "/Data_2016-06-22/genes_description_full.xls"), sep="\t", quote=F, row.names=F)




genes.description <- read.table(paste0(data_dir, "/Data_2016-06-22/genes_description_full.xls"), sep="\t", header = TRUE, as.is = TRUE)
head(genes.description)









