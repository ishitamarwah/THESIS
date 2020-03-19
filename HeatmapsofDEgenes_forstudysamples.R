#================================================
# Created by - Dr. Ishita Marwah
# Created on - 01/01/2020
# This script is for generating heatmaps for DEGs from DE analyses performed using getDEGenesUsingEdgeR_forstudysamples.R
# code adapted from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4934518/
#================================================


#Clear the workspace
WS = c(ls())
rm(WS, list = WS)
options("scipen"=100, "digits"=2)


library(org.Hs.eg.db)
library(stringi)
library(biomaRt)
library(edgeR)
library(corrplot)
library(psych)
library(vctrs)
library(pillar)
library(ggplot2)
library(ggfortify)
library(plyr)
library(gplots)
library(RColorBrewer)
library(ggrepel)
library(dplyr)
library(GOexpress)

#for MAC

mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org", path="/biomart/martservice")
ensembl <- useDataset("hsapiens_gene_ensembl",mart)

#for linux workstation

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org", path="/biomart/martservice")
ensembl <- useDataset("hsapiens_gene_ensembl",mart)

dir.create("./Images", showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create("./Images/DE_Heatmaps", showWarnings = TRUE, recursive = FALSE, mode = "0777")

locCountFiles<-"C:/Users/im2214/Desktop/Patient Runs/Patient_Run_June_2019_glmQL" #see main getDEGenesUsingEdgeR script for details
setwd(locCountFiles)
locCountFiles<- paste0(getwd(),"/")

locRObj<-"./DET_RObjects/"

#------------------------------
#COMMON FUNCTIONS
#------------------------------

#------------------------------------------------
# Derive HGNC gene annotations for ensembl identifiers in a dataframe
#------------------------------------------------

getHGNCSymbolsForDETs <-function(df_dat)
{
  if (!"transcript_id" %in% colnames(df_dat))
  { 
    df_dat$transcript_id <-rownames(df_dat)
  }
  G_list <- getBM(filters= "ensembl_transcript_id_version", attributes= c("ensembl_transcript_id_version","hgnc_symbol"),values= unique(df_dat$transcript_id),mart= ensembl)
  names(G_list)[1]<-"transcript_id"
  # merge with passed df_dat
  annotated_df_dat=join(df_dat,G_list,by="transcript_id", match="first")
  #write.csv(annotated_df_dat,paste0(locCountFiles,fileName))
  return(annotated_df_dat)
}

#-------------------------------------------------
#Function to read in RDS files or R Objects saved from DE analyses
#-------------------------------------------------

readDET_RObjects<-function(sampleSource1,sampleSource2,sampleType1,sampleType2)
{
  fitObjecttrial<-readRDS(paste0(locRObj,"list_",sampleSource1,"_",sampleType1,"_Vs_",sampleSource2,"_",sampleType2,".rds"))
  return(fitObjecttrial)
}

#-------------------------------------------------
#This function  helps calculate the total number of significant transcripts from decideTests.DGELRT()
#which can then be used to extract these transcripts from topTags
#-------------------------------------------------

getNumSigTranscripts=function(dgelrt)
{
  status=decideTests.DGELRT(dgelrt)
  Down=summary(status)[1,1]
  Up=summary(status)[3,1]
  n= (Down + Up)
  
  return(n)
}

#-------------------------------------------------
#PLOTTING DE HEATMAPS
#-------------------------------------------------

# In order to make heatmaps of significant annotated DEGs (genes not transcripts) imported from IPA do the following:
# set sampleSource and SampleType and SampleType for DE comparison of interest, for example for EU vs EI BULK BAL DEG heatmap
# set sampleSource1 = "EU_BAL"; sampleSource2 ="EI_BAL"; sampleType1= "BULK"; sampleType2 = "BULK"

sampleSource1=""
sampleSource2=""
sampleType1=""
sampleType2=""

# read in the saved R object for the DE comparison of interest

fitObject=readDET_RObjects(sampleSource1,sampleSource2,sampleType1,sampleType2)
y=fitObject$y
y=calcNormFactors(y)


#calculate all significant topTags

top.TAG=topTags(fitObject$fitTreat1.2, n=getNumSigTranscripts(fitObject$fitTreat1.2))

#logcpm and manipulation and preparation for plotting

logCPM <- cpm (y, prior.count =2, log=TRUE)
rownames(logCPM) = y$genes$genes
colnames(logCPM) = paste(rownames(y$samples))
col.pan =colorpanel(100, "blue","white","red")
o = rownames(top.TAG$table)
logCPM = logCPM[o,]
logCPM= t(scale(t(logCPM)))
mergefile= as.data.frame(logCPM)


#Read in the annotated csv file for the same comparison for IPA, for example:

IPAdegFile=read.csv("../IPA output/20190810_unfiltered/Analysis-ready datasets from IPA/EU_EI_BAL_BULK_fc1.2.csv") # this is 
# an example, make sure to read in your file of interest

#Extract only gene symbols and ensembl id (which will be used to match to mergefile)

rownames(IPAdegFile)=IPAdegFile$ID
IPAdegFile=select(IPAdegFile,Symbol)

#merge mergefile with IPAdegFile

mergefile=merge(IPAdegFile,mergefile, by="row.names")
rownames(mergefile)=mergefile$Symbol
mergefile=mergefile[,-(1:2)]

logCPM=as.matrix(mergefile)


main=paste0("All ",nrow(unique(mergefile)), " significant DEGs ",sampleSource1,"_",
            sampleType1,"_VS_",sampleSource2,"_",sampleType2)

heatmap.2(logCPM, col=col.pan,Rowv =TRUE, scale="none", trace= "none", dendrogram ="both", 
          main=main,cexRow =1, cexCol=1.4,density.info ="none",margin=c(16,16),lhei=c(2,10),lwid=c(2,6))

#Saving the heatmap image

newdirectory=paste0("./Images/DE_Heatmaps/",sampleSource1,"_VS_",sampleSource2,"/")
dir.create(newdirectory, showWarnings = TRUE, recursive = FALSE, mode = "0777")
ImageName=paste0(newdirectory,"Heatmap_all_significant_DEGs_",sampleSource1,"_",sampleType1,"_VS_",sampleSource2,"_",sampleType2,".png")
png(filename = ImageName,  width = 900, height = 650)
par(cex.main=0.65)
heatmap.2(logCPM, col=col.pan,Rowv =TRUE, scale="none", trace= "none", dendrogram ="both", 
          main=main,cexRow =1, cexCol=1.1,density.info ="none",margin=c(16,16),lhei=c(2,10),lwid=c(2,6))
dev.off()

#=========
