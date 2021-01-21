# load libraries
library(readr)
library(dplyr)
library(magrittr)
library(tximport)
library(DESeq2)
# set working directory
setwd("~/Downloads/test_RNAseq/")
# load metadata
sample_table <-read.table("metadata.csv", sep = "\t", header = T)
# load files paths
sample_files = paste0(pull(sample_table, "sample_name"), "/quant.sf")
# name table columns
names(sample_files) = pull(sample_table, "sample_name")
# relate genes to transcripts
tx2gene = read.table("txgen.txt", sep = "\t", col.names =c("genid","transid"))
# import count data to tximport
count_data = tximport( files = sample_files,
          type = "salmon",
          tx2gene =  tx2gene,
          ignoreTxVersion = T)
# change "~ 1" to "~ condition"
data <- DESeqDataSetFromTximport(txi = count_data,
                                 colData = sample_table,
                                 design = ~ 1)


