library(janitor
        )
rm(list=ls())
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(pheatmap)
library(readr)
library(dplyr)
library(Rtsne)
library(umap)
library(viridisLite)
library(uwot)
library(Seurat)
library(Matrix)
library(dplyr)
library(stringr)
library(GEOquery)
library(tidyverse)
library(readr)
library(readxl)
library(limma)
library(DESeq2)
#got 32K genes, 

setwd("/Users/u0166682/Desktop/projects/temporal and spatial expression DEEs/zebrafish")
list.files()
wt.data<-read.delim(file = "wt.data.txt", header = TRUE, row.names=1)
wt.data.mat<-as.matrix(wt.data)
hist(wt.data.mat)
wt.meta<-read.delim(file = "wt.meta.txt", header = TRUE, row.names = 1)
colnames(wt.data)==row.names(wt.meta)
#want a normalised matrix to extract the values and plot them. do data exploration to start
#through deseq2 can normalise and extract different versions (way to go)


#removes decimals for deseq2
wt.data.round<-round(wt.data)
my.design<-wt.meta
ddsFullCountTable<- DESeqDataSetFromMatrix(
  countData = wt.data.round,
  colData = my.design,
  design = ~1     
)
#normalise matrix
dds<-DESeq(ddsFullCountTable)

vst<-vst(dds, blind = TRUE)
expression<- assay(vst)
head(expression)

#canj start extracting zebrafish genes
#start at getting normalised counts and create a log2(cpm) matrix. so have 3 diff expression matrixes.