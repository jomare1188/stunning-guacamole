# load libraries

#install.packages(c("BiocManager","readr","dplyr", "magrittr","ggplot2"),
# lib = "./",
# repos = 'http://cran.us.r-project.org')

#BiocManager::install(c("tximport","DESeq2"), lib = "./")

library(readr, lib.loc =  "./")
library(dplyr, lib.loc =  "./")
library(magrittr,  lib.loc = "./")
library(tximport, lib.loc = "./")
library(DESeq2, lib.loc = "./")
library(ggplot2, lib.loc =  "./")

# set working directory
setwd("./")
# load metadata
sample_table <-read.table("metadata_complete.csv", sep = ",", header = T)
# load files paths
sample_files = paste0(pull(sample_table , "Sample_file"), "/quant.sf")
# name table columns
names(sample_files) = pull(sample_table, "Sample_file")
# relate genes to transcripts
tx2gene = read.table("txgen2.txt", sep = "\t", col.names =c("genid","transid"))
# import count data to tximport
count_data = tximport( files = sample_files,
          type = "salmon",
          tx2gene =  tx2gene,
          ignoreTxVersion = T)
# change "~ 1" to "~ condition"
# convert condition to a factor
sample_table$Condition <- as.factor((sample_table$Condition))
data <- DESeqDataSetFromTximport(txi = count_data,
                                 colData = sample_table,
                                 design = ~ Condition)


# 
data <- estimateSizeFactors(data)
normalizationFactors(data)
counts(data, normalized = TRUE)
vst <- varianceStabilizingTransformation(data)
boxplot(assay(vst))

pca <- plotPCA(vst, intgroup ='Condition') +
  theme_bw() +
  theme( text = element_text(size=22), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
 	 panel.grid.minor = element_blank(),
         axis.line = element_line(colour = "black") )

ggsave( pca, filename = "./PCA_NRG.png",units = "cm",width = 20*1.3, height = 20,dpi = 320)

plotPCA(vst, intgroup ='Condition') +
  theme_bw() +
  theme(text = element_text(size=22),
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.line = element_line(colour = "black"))


#
