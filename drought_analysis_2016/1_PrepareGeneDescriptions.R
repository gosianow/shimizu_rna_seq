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


### AT gene matches 

AT.id <- read.table(paste0(data_dir, "/Data_2016-06-22/orig/Blastp_result_against_Athaliana.txt"), sep="\t", stringsAsFactors = FALSE, header = FALSE)
head(AT.id)

AT.id$at_id <- strsplit2(AT.id[,2], "\\.")[ ,1]


### Get AT symbols

syms <- mget(AT.id[, "at_id"], org.At.tairSYMBOL)
syms <- sapply(syms, function(u) if(is.na(u[1])) NA else paste(u,collapse=","))
AT.id$at_symbol <- syms


### Get functional description from TAIR10

tfd <- read.table(paste0(data_dir, "/Data_2016-06-22/orig/TAIR10_functional_descriptions.txt"), sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)


t10id <- gsub("\\.[1-9]", "", tfd$Model_name)
m <- match(AT.id[, "at_id"], t10id)

AT.id$description <- as.character(tfd$Short_description[m])

### Remove problematic symbols from description
AT.id$description <- gsub("[^[:alnum:],()+/\\-]", " ", AT.id$description)
AT.id$at_symbol <- gsub("[^[:alnum:],()+/\\-]", " ", AT.id$at_symbol)


genes.description <- AT.id[, c("V1", "at_id", "at_symbol", "description")]

colnames(genes.description) <- c( "ID", "AT_ID", "AT_symbol", "Description")

write.table(genes.description, paste0(data_dir, "/Data_2016-06-22/genes_description.xls"), sep="\t", quote = FALSE, row.names = FALSE)


### Flowering genes - Table S4 from MolEcol paper

Athaliana.flowering.genes <- read.table(paste0(data_dir,"/GeneControlSets/Flowering_Genes/gene_list_flowering_related_MolEcol_Athaliana_S4.csv"), sep=",", header=T, stringsAsFactors=FALSE, skip=1)

# Some clean up 
Athaliana.flowering.genes[,1] <- gsub(" ", "", Athaliana.flowering.genes[,1])
colnames(Athaliana.flowering.genes) <- c("AT_ID", "AT_NAME")
Athaliana.flowering.genes$Flowering <- TRUE



### Drought responce genes 

UP.genes <- read.table(paste0(data_dir,"/GeneControlSets/Drought_Genes/mDr_Day10_drought_up_regulated_genes_Harb_etal.csv"), sep=";", header=T, stringsAsFactors=FALSE)

DOWN.genes <- read.table(paste0(data_dir,"/GeneControlSets/Drought_Genes/mDr_Day10_drought_down_regulated_genes_Harb_etal.csv"), sep=";", header=T, stringsAsFactors=FALSE)


Drought.genes <- rbind(UP.genes[, c("Gene", "mDr.Day10")], DOWN.genes[, c("Gene", "mDr.Day10")])

names(Drought.genes) <- c("AT_ID", "Drought_regulation")



### Merge gene descriptions

genes.full.description <- merge(genes.description, Drought.genes, by="AT_ID", all.x=TRUE, sort = FALSE)


genes.full.description <- merge(genes.full.description, Athaliana.flowering.genes, by = "AT_ID", all.x = TRUE, sort = FALSE)



write.table(genes.full.description, paste0(data_dir, "/Data_2016-06-22/genes_description_full.xls"), sep="\t", quote=F, row.names=F)




genes.description <- read.table(paste0(data_dir, "/Data_2016-06-22/genes_description_full.xls"), sep="\t", header = TRUE, as.is = TRUE)
head(genes.description)









