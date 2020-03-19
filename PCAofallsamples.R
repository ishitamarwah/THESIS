#===============================================
# Created by - Dr. Ishita Marwah
# Created on - 01/01/2020
# SCRIPT for combining counts from different sources to plot PCAs
#===============================================

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
library(ggplot2)
library(ggfortify)
library(plyr)
library(gplots)
library(RColorBrewer)
library(dplyr)
library(readODS)

locCountFiles<- paste0(getwd(),"/") #see main getDEGenesUsingEdgeR script for details

locRObj<-"./DET_RObjects/"

dir.create("./Images", showWarnings = TRUE, recursive = FALSE, mode = "0777")

#--------------------------------------------------------------------------------
# Funtion to read in raw counts for futher combination and manipulation
#--------------------------------------------------------------------------------

readRawCounts<-function(sampleSource,sampleType)
{
  rawcountfile<-read.csv(paste0(locCountFiles,sampleSource,"_DECounts/",sampleSource,"_transcript_count_matrix_",sampleType,".csv"))
  rownames(rawcountfile)=rawcountfile$transcript_id 
  tag=paste0(sampleSource,"_",sampleType,"_")
  colnames(rawcountfile) <- paste0(tag, colnames(rawcountfile))
  rawcountfile=rawcountfile[,-1]
  if (!"transcript_id" %in% colnames(rawcountfile))
  { 
    rawcountfile$transcript_id <-rownames(rawcountfile)
  }
  rawcountfile=select(rawcountfile, transcript_id, everything())
  rownames(rawcountfile)=NULL
  return(rawcountfile)
}

#---------------------------------------------------------------------------------
# Fucntions to combine counts from files that are read in
#---------------------------------------------------------------------------------

# Function to combine all sampleType counts for each BAL-related sampleSource

compileforSampleSource <- function(sampleSource)
{
  
  R1=readRawCounts(sampleSource,"BULK")
  R2=readRawCounts(sampleSource,"CD4")
  R3=readRawCounts(sampleSource,"CD8")
  R4=readRawCounts(sampleSource,"CD206")
  R1R2=join(R1,R2,by="transcript_id")
  R3R4=join(R3,R4,by="transcript_id")
  
  combined=join(R1R2,R3R4,by="transcript_id")
  
  return(combined)
}


# Function to combine all sampleType counts for each BLOOD-related sampleSource

compileforSampleSource_blood <- function(sampleSource)
{
  
  R1=readRawCounts(sampleSource,"BULK")
  R2=readRawCounts(sampleSource,"CD4")
  R3=readRawCounts(sampleSource,"CD8")
  R1R2=join(R1,R2,by="transcript_id")
  
  combined=join(R1R2,R3,by="transcript_id")
  
  return(combined)
}

#=====================================================================================================

# PCA for all 100 study samples/sequences used for DE analysis

#=====================================================================================================

# Compiling per sampleSource to get a masterfile with all counts

AllEIBAL=compileforSampleSource("EI_BAL")
AllEUBAL=compileforSampleSource("EU_BAL")
AllUUBAL=compileforSampleSource("UU_BAL")

AllEUBlood=readRawCounts("EU_BLOOD","CD4")
AllEIBlood=readRawCounts("EI_BLOOD","CD4")
AllUUBlood=readRawCounts("UU_BLOOD", "CD4")

A=join(AllUUBAL,AllUUBlood, by="transcript_id")
B=join(A,AllEIBAL,by="transcript_id")
C=join(B,AllEIBlood,by="transcript_id")
D=join(C,AllEUBAL,by="transcript_id")
E=join(D,AllEUBlood,by="transcript_id")

#convert masterfile dataframe into DGElist format

E[is.na(E)]=0
rownames(E)=E$transcript_id
E=E[,-1]
y=DGEList(E,remove.zeros=T)

y=calcNormFactors(y)

# Assigning Study Group Type
y$samples$sampleGroupType=rownames(y$samples)
y$samples$sampleGroupType=gsub("_T.*", "\\1", y$samples$sampleGroupType)
y$samples$sampleGroupType=as.factor(y$samples$sampleGroupType)

# PCA
mydata1=y$samples
mydata0=cpm(y,log=T)
mydata2=as.data.frame(t(mydata0))
de=merge(mydata1,mydata2, by="row.names", all="TRUE")
df= de[-c(1:10)]
rownames(de)=de$Row.names
pca <- prcomp(df)

COLOURS=c("orange","orange","orange","orange","brown","purple","purple","purple","purple","magenta",
          "green3","green3","green3","green3","turquoise4")

SHAPES=rep(1:5, times =3)

# Plot PCA
autoplot(pca, data=de) +
  geom_point(aes(shape=sampleGroupType, color=sampleGroupType, size=sampleGroupType), size=3.5, stroke = 2) +
  scale_shape_manual(values=SHAPES) + 
  scale_color_manual(values=COLOURS)

#========================================================================================================

