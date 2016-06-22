
AT.id <- read.table("Data/genes_descr_control/best_hit_blast_result_Sl_predicted_exons.txt", sep=",", stringsAsFactors=FALSE)
head(AT.id)



library(org.At.tair.db)

syms <- mget(AT.id[,2], org.At.tairSYMBOL)
syms <- sapply(syms, function(u) if(is.na(u[1])) NA else paste(u,collapse=","))
AT.id$at_symbol <- syms

tfd <- read.table("Data/genes_descr_control/TAIR10_functional_descriptions",sep="\t", header=TRUE, quote="", comment.char="", stringsAsFactors=FALSE)

t10id <- gsub("\\.[1-9]","",tfd$Model_name)
m <- match(AT.id[,2],t10id)

AT.id$description <- as.character(tfd$Short_description[m])


genes.description <- merge(data.frame(rownames(x)), AT.id, by=1, all.x=T)

genes.description <- genes.description[,c(1,2,6,7)]

colnames(genes.description) <- c( "ID", "AT_ID", "AT_symbol", "Description")

write.table(genes.description, "Data/genes_descr_control/genes_description_full.xls", sep="\t", quote=F, row.names=F)

write.table(genes.description, "Data/genes_descr_control/genes_description_full.txt", sep="\t", quote=F, row.names=F)


# genes.x <- rownames(x)
# genes.full.description <- data.frame(ID=genes.x)
# genes.full.description <- merge(genes.full.description, genes.description, by=1, all.x=TRUE)
# genes.full.description <- merge(genes.full.description, cbind(Flowering="FL", Athaliana.flowering.genes), by.x=2, by.y=2, all.x=TRUE)
# colnames(genes.full.description) <- c("AT_ID", "ID", "AT_symbol", "Description", "Flowering", "Description2")
# genes.full.description <- genes.full.description[,c("ID", "AT_ID","AT_symbol", "Description", "Flowering", "Description2")]
# write.table(genes.full.description, "Data/genes_descr_control/genes_description_full.csv", sep=";", quote=F, row.names=F)





