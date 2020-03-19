#================================================
# Created by - Dr. Ishita Marwah
# Created on - 01/01/2020
# This script reads in the raw transcript count csv files obtained from prepDE.py to perform
# differential expression using the edgeR (edgeR_3.22.5) package for RNASeq analysis. This was used for the entirety of
# the study except for the the TOPFISH validation analyses, which are covered by ScriptforTOPFISHvalidation.R
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
library(plyr)

#for MAC OS
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org", path="/biomart/martservice")
ensembl <- useDataset("hsapiens_gene_ensembl",mart)

#OR

#for linux workstation
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org", path="/biomart/martservice")
ensembl <- useDataset("hsapiens_gene_ensembl",mart)

# <- <- <- <- <- <- <- <- <- <- <- <- <- <- 
#
# Read Input gene and transcript counts
#
# <- <- <- <- <- <- <- <- <- <- <- <- <- <- 


locCountFiles<-"C:/Users/im2214/Desktop/Patient Runs/Patient_Run_June_2019_glmQL" # OR Use location where files should be output and
# transcript count files are located within indiviudal sub-folders labelled 'studygroup_BAL/BLOOD_DECounts'

setwd(locCountFiles)

locCountFiles<- paste0(getwd(),"/")

dir.create("./DET_RObjects", showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create("./FC1.2", showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create("./FC1.2//Data_ModuleAnalysis_fc1.2", showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create("./FC1.2/DET_Files_qlf_fc1.2", showWarnings = TRUE, recursive = FALSE, mode = "0777")


locRObj<-"./DET_RObjects/"

#------------------------------------------
# This function reads in transcript count files generated with prepDE.py. 
# sampleSource = samplegroup_BAL or samplegroup_BLOOD i.e. "EU_BAL" EI_BAL" or "UU_BAL" and similarly for blood, 
# sampleType = BULK or cell-type i.e CD206, CD4, CD8.
#------------------------------------------

readDECountFile<-function(sampleSource,sampleType)
{
		df_dat<-read.csv(paste0(locCountFiles,sampleSource,"_DECounts/",sampleSource,"_transcript_count_matrix_", 
		                        sampleType,".csv"),header = T, stringsAsFactors = F)
		return(df_dat)
}


#----------------------------------------------------
# This function is to get corresponding HGNC symbols for a data 
# frame containing ensembl transcript identifiers. 
#----------------------------------------------------

getHGNCSymbolsForDETs <-function(df_dat)
{
  if (!"transcript_id" %in% colnames(df_dat))
  { 
    df_dat$transcript_id <-rownames(df_dat)
    }
  G_list <- getBM(filters= "ensembl_transcript_id_version", attributes= c("ensembl_transcript_id_version","hgnc_symbol"),
                  values= unique(df_dat$transcript_id),mart= ensembl)
names(G_list)[1]<-"transcript_id"
# merge with passed df_dat
annotated_df_dat=join(df_dat,G_list,by="transcript_id", match="first")
#write.csv(annotated_df_dat,paste0(locCountFiles,fileName))
return(annotated_df_dat)
}

#----------------------------------------------------
# This function extracts original counts for all samples for a selected list of transcripts
#----------------------------------------------------

getDataForModuleAnalysis<-function(selectedDETs, transCounts)
{
  selectedDETs$transcript_id <-rownames(selectedDETs)
  transCounts$transcript_id <-rownames(transCounts)
  #trans_merge <- merge(selectedDETs, transCounts, by = 'transcript_id' ,all.x = TRUE)
  trans_merge=join(selectedDETs, transCounts,by="transcript_id", match="first")
  rownames(trans_merge)=NULL
  return(trans_merge)
} 

#------------------------------------------------------
# This function helps calculate the total number of significant transcripts 
# using decideTests.DGELRT() which can then be used to extract these transcripts from topTags
#------------------------------------------------------

getNumSigTranscripts=function(dgelrt)
{
  status=decideTests.DGELRT(dgelrt)
  Down=summary(status)[1,1]
  Up=summary(status)[3,1]
  n= (Down + Up)
  
  return(n)
}

#----------------------------------------------------
# This function helps perform differential expression analysis with edgeR and incorporated all the functions above
# Differential expression is performed alphabetically i.e. uses the alphabetically prior group as baseline - for example -
# comparison between EI and EU will by default give transcripts differentially expressed in EU relative to EI. A condition for 
# reverseComparison is worked into the function below to reverse the default comparison direction if desired
#----------------------------------------------------

getDET<-function(locCountFiles, sampleSource1,sampleSource2, sampleType1, sampleType2, reverseComparison=FALSE, removeTrizol=FALSE)
{
	# read transcript data
	transS1 <-readDECountFile(sampleSource1, sampleType1)
	transS2 <-readDECountFile(sampleSource2, sampleType2)
	
	names(transS1) <- paste0(names(transS1),paste0(".", sampleSource1,"_",sampleType1))
	names(transS2) <- paste0(names(transS2),paste0(".", sampleSource2,"_",sampleType2))
  colnames(transS1)[1]="transcript_id"
	colnames(transS2)[1]="transcript_id"
	
	trans_merge=join(transS1,transS2,by="transcript_id", match="first")
	
	rownames(trans_merge)<-trans_merge[,1]
	trans_merge<-trans_merge[,-1]
	
	trans_merge[is.na(trans_merge)] <- 0
	
	sampleGp<-unlist(lapply(strsplit(paste(names(trans_merge)),".",fixed=T),function(l) l[2]))
	
	y <- DGEList(counts= trans_merge, group= sampleGp, genes = row.names(trans_merge), remove.zeros=T)
	
	#Filter out lowly expressed genes
	
	keep= filterByExpr.default(y, min.count=3)
	y <- y[keep, , keep.lib.sizes=FALSE]
	
	outFileName<-paste0(sampleSource1,"_",sampleType1,"_Vs_",sampleSource2,"_",sampleType2)

	if(reverseComparison)
	{
	  sampleGp<-as.factor(sampleGp)
	  revSampleGp<-rev(sampleGp)
	  levels(revSampleGp)<-rev(levels(sampleGp))
	  design <- model.matrix(~ revSampleGp)
	}
	else{
	  design <- model.matrix(~ sampleGp)
	}
	
	
	y <- estimateDisp(y, design, robust=TRUE)
	
	fit <- glmQLFit(y, design, robust=TRUE)
	
	#all differntially expressed transcripts tested for absolute fold change more than 1.2
	
	fitTreat1.2<-glmTreat(fit,lfc = log2(1.2))

	# Save the R data objects with full calculations about differential expression as a list
	
	listDETRObjects<-list()
	
	listDETRObjects[["y"]]<- y
	listDETRObjects[["qlf"]]<- qlf
	listDETRObjects[["fitTreat1.2"]]<- fitTreat1.2
	
	
	# Filter fitTreat1.2 by topTags to get the significant transcripts with FDR < 0.05% 

	sig_transcripts_fitTreat_fc1.2=topTags(fitTreat1.2, n=getNumSigTranscripts(fitTreat1.2))
	
	# For reference, fetch the original transcript counts per sample for the significant transcripts identified above 
	
	df_dataForModuleAnalysis_fc1.2= getDataForModuleAnalysis(sig_transcripts_fitTreat_fc1.2$table, trans_merge)
	
	# save signficant transcript info as a csv file for IPA analyis 
	
	write.csv(df_dataForModuleAnalysis_fc1.2, paste0("./FC1.2/Data_ModuleAnalysis_fc1.2/transCount_fc1.2_",outFileName,".csv"))
	
	# For future exploratory purposes including for tmod modular analysis, also export all transcripts used for
	# DE analysis (after clearing edgeR filters) as csv files, after annotation with HGNC symbols. Also save this as
	# part of the listDETRObjects generated above, for future reference and manipulation if required
	
	df_Annot_1.2 = getHGNCSymbolsForDETs(fitTreat1.2$table)
	
	write.csv(df_Annot_1.2, paste0("./FC1.2/DET_Files_qlf_fc1.2/DET_qlf_fc1.2_",outFileName,".csv"))
	
	listDETRObjects[["DET_HGNCAnnot_qlf_fitTreat1.2"]]<- df_Annot_1.2
	
	saveRDS(listDETRObjects, paste0("./DET_RObjects/list_",outFileName,".rds"))
  
	return(df_Annot)
}

# <- <- <- <- <- <- <- <- <- <- <- <- <- <-  <- 
#
# END OF FUNCTIONS USED FOR THE FOLLOWING ANALYSES
#
# <- <- <- <- <- <- <- <- <- <- <- <- <- <-  <- 

#---------------------------------------------------------------
# Performing DE analyses for the study 
#---------------------------------------------------------------

# All direct EU versus EI comparisons

DET_BAL_CD206_EUvsEI<-getDET(locCountFiles, "EU_BAL", "EI_BAL", "CD206", "CD206")

DET_BAL_CD4_EUvsEI<-getDET(locCountFiles, "EU_BAL", "EI_BAL", "CD4", "CD4")

DET_BAL_CD8_EUvsEI<-getDET(locCountFiles, "EU_BAL", "EI_BAL", "CD8", "CD8")

DET_BAL_BULK_EUvsEI<-getDET(locCountFiles, "EU_BAL", "EI_BAL", "BULK", "BULK")

DET_Blood_CD4_EUvsEI<-getDET(locCountFiles, "EU_BLOOD", "EI_BLOOD", "CD4", "CD4")


#all EU and EI individual comparisons versus UU i.e UU as the comparator

DET_BAL_CD206_EUvsUU<-getDET(locCountFiles, "EU_BAL", "UU_BAL", "CD206", "CD206", TRUE)
DET_BAL_CD206_EIvsUU<-getDET(locCountFiles, "EI_BAL", "UU_BAL", "CD206", "CD206", TRUE)

DET_BAL_CD4_EUvsUU<-getDET(locCountFiles, "EU_BAL", "UU_BAL", "CD4", "CD4", TRUE)
DET_BAL_CD4_EIvsUU<-getDET(locCountFiles, "EI_BAL", "UU_BAL", "CD4", "CD4", TRUE)

DET_BAL_CD8_EUvsUU<-getDET(locCountFiles, "EU_BAL", "UU_BAL", "CD8", "CD8", TRUE)
DET_BAL_CD8_EIvsUU<-getDET(locCountFiles, "EI_BAL", "UU_BAL", "CD8", "CD8", TRUE)

DET_BAL_BULK_EUvsUU<-getDET(locCountFiles, "EU_BAL", "UU_BAL", "BULK", "BULK", TRUE)
DET_BAL_BULK_EIvsUU<-getDET(locCountFiles, "EI_BAL", "UU_BAL", "BULK", "BULK", TRUE)

DET_BLOOD_CD4_EIvsUU<-getDET(locCountFiles, "EI_BLOOD", "UU_BLOOD", "CD4", "CD4", TRUE)
DET_BLOOD_CD4_EUvsUU<-getDET(locCountFiles, "EU_BLOOD", "UU_BLOOD", "CD4", "CD4", TRUE)

