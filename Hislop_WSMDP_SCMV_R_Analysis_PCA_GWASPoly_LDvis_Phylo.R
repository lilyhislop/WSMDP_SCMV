#########################
#This code was written by Lillian Hislop
#with reference to code written by Jenyne Loraca(GWASPoly) and Matheus Bassegio (LD Visualization)
#2020.2.24-
#the purpose of this code is to run additional/complementary GWAS on the SCMV Panel conducted in 2019.
#By using GWASpoly package. Previous code was written using FarmCPU
#########################

#########################
###Establish Workspace###
#########################
#begin by establishing a new fresh work space
rm(list=ls())

#load in the libraries
library(GWASpoly) # imported with a .zip file. #for GWASPoly
library(tidyverse) #For manipulating and renaming dataframes
library(data.table) #for the LD Visualization and datatable wrangling
library("scrime")#for recodeSNPS
library(lme4)#for linear modeling
library(ggtree)#for phylo tree
library(summarytools)
# library(gdsfmt)
library(SNPRelate) #For PCA and Visualizations
library("compiler") #needed to make GAPIT work
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")
library("bigmemory") #to make a matrix big
library(rrBLUP)

#establish the working directory
setwd("C:/Users/CaptainsScreen/Documents/0 Grad School/0 Lab/Diversity Panel/SCMV")
source("RCode/PCAFigureCreation.R")
source("RCode/PhenowithGenoPrep.R")
source("RCode/hmpToNumeric.R")
source("RCode/GWASPolyVis.R")
source("RCode/GWASPolyRunner.R")
source("RCode/WritePhenoGenotoFile.R")


#########################
###Phenotype Data Import###
#########################
phenoFile <- "Data/WSMDP_SCMV_Data_DifferentRankingScales.csv"
pheno <- read.delim(phenoFile, header = TRUE, sep = ",", dec = ".")
phenoSubsetGeno<- PhenowithGenoPrep(phenoFile)
taxalistwcol <- read.delim("Data/WSMDP_SCMV_420taxaPheno_wGenoData.txt",header = FALSE, sep = '\n')


#########################
###SNP Relate Establishment###
#########################
filename <- "WSMDP_SCMV_SeqB"


#read in VCF. This file has been previously filtered for 90% SNP call rate and MinorAlleleFrequency of 0.025
filename
vcfpath <- paste("Data/SequenceData/",filename,".vcf",sep = "")
snpgdsVCF2GDS(vcfpath, "Data/test.gds", method = "biallelic.only")
snpgdsSummary("Data/test.gds")
SCMV <- snpgdsOpen("Data/test.gds")

#Lets prune anything with an LD over 0.98
SCMV_MAF_LD <- snpgdsLDpruning(SCMV, ld.threshold = 0.98, start.pos = "first", verbose = TRUE)
names(SCMV_MAF_LD)
head(SCMV_MAF_LD$chr1)
SCMV_MAF_LD.id <- unlist(SCMV_MAF_LD)

#########################
### PCA ###
#########################
#Read in the SCMV gdsobj with only the snp.id's found by LD pruning and only the sample.id's that we phenotypes
SCMV_PCA <- snpgdsPCA(SCMV, sample.id = taxalistwcol$V1, snp.id = SCMV_MAF_LD.id)
# SCMV_PCA <- snpgdsPCA(SCMV, sample.id = phenowcol$V1) 
#cut off the numbers at the end of the sample.id that don't mean things to humans
SCMV_PCA$sample.id<-gsub(SCMV_PCA$sample.id, pattern = ":.*", replacement = "")

pc.perc <- SCMV_PCA$varprop*100
head(round(pc.perc,2))


PCAFigureCreation(SCMV_PCA,pc.perc,phenoSubsetGeno,filename,"EndospermType")
PCAFigureCreation(SCMV_PCA,pc.perc,phenoSubsetGeno,filename,"Program")
PCAFigureCreation(SCMV_PCA,pc.perc,phenoSubsetGeno,filename,"Region")
PCAFigureCreation(SCMV_PCA,pc.perc,phenoSubsetGeno,filename,"SusceptibilityRating03")
PCAFigureCreation(SCMV_PCA,pc.perc,phenoSubsetGeno,filename,"SusceptibilityRating02")
PCAFigureCreation(SCMV_PCA,pc.perc,phenoSubsetGeno,filename,"SusceptibilityRating05")
PCAFigureCreation(SCMV_PCA,pc.perc,phenoSubsetGeno,filename,"PercentInfectedAllRounds")

#########################
### Close Snp Relate ###
#########################
adendum <- "_LD098_PostSNPRelate_420taxa"
outfile <- paste("Data/SequenceData/",filename, adendum ,sep="")
snpgdsGDS2PED(SCMV, outfile, sample.id = taxalistwcol$V1, snp.id = SCMV_MAF_LD.id)
snpgdsClose(SCMV)

 
#########################
### GWASPoly ###
#########################
#read in genetic info post MAF and LD pruning. Pruning done by SNPrelate package. Outputted previously as a plink format, converted to Hapmap by tassel and read back in
hmppath <- paste("Data/SequenceData/",filename,adendum,".hmp.txt",sep = "")
SCMV_geno <- fread(hmppath,skip = "rs#")
geno_scmv <- SCMV_geno

str(geno_scmv)
colnames(geno_scmv)<-gsub(colnames(geno_scmv), pattern = ":.*", replacement = "")
str(geno_scmv)


#########################
### GWASPoly Function ###
#########################
#trait to analyze
trait <- "SusceptibilityRating02"
GWASPolyRunner("FullGroup_420lines_PC3_Functioned", trait, geno_scmv,filename,adendum,phenoSubsetGeno)
GWASPolyRunner("FullGroup_420lines_PC3_Functioned", "SusceptibilityRating03", geno_scmv,filename,adendum,phenoSubsetGeno)
GWASPolyRunner("FullGroup_420lines_PC3_Functioned", "SusceptibilityRating05", geno_scmv,filename,adendum,phenoSubsetGeno)
GWASPolyRunner("FullGroup_420lines_PC3_Functioned", "PercentInfectedAllRounds", geno_scmv,filename,adendum,phenoSubsetGeno)





#############
#data visualization
#############
pheno$PotLabel <- as.factor(pheno[,1])
pheno$StandCountCleaned <- as.numeric(pheno$StandCountCleaned)
pheno$SusceptibilityRating02 <- as.factor(pheno$SusceptibilityRating02)
pheno$SusceptibilityRating03 <- as.factor(pheno$SusceptibilityRating03)
pheno$SusceptibilityRating05 <- as.factor(pheno$SusceptibilityRating05)

summary(pheno)
summary(pheno$PercentInfectedAllRounds)
dfSummary(pheno)
freq(pheno)
descr(pheno)#descr shows that plant height has outlier with kurtosis
ctable(pheno)

plotData <- function(ToPlot, type, ByPlot=NULL){
  png(paste("Figures/Plot_Pheno_",ToPlot,"_",type,".png",sep=""))
  if(type == "hist"){
    pheno[,ToPlot]<-as.numeric(pheno[,ToPlot])
    hist(pheno[,ToPlot], main = paste(ToPlot,"Freq"), xlab = ToPlot)
}
  if(type == "plot"){plot(pheno[,ToPlot]~pheno[,ByPlot], main = paste(ToPlot,"by",ByPlot))}
  dev.off()
}

plotData("PercentInfectedAllRounds","hist")
plotData("SusceptibilityRating05","hist")
plotData("SusceptibilityRating03","hist")
plotData("SusceptibilityRating02","hist")


png("Figures/Plot_Susceptible.png")
barplot(prop.table(table(pheno$Susceptible)), xlab = "Susceptibility Scores",ylab = "Frequency of Occurance",ylim = c(0,1), main = "Susceptibility Scores within Panel")
dev.off()
png("Figures/Plot_Suscep_by_Pot.png")
plot(x =  pheno$PotLabel, y = pheno$Susceptible, main = "Susceptibility by Pot label")
dev.off()
png("Figures/Plot_StandCount_Hist.png")
hist(pheno$StandCountCleaned, main = "Stand Count Freq", xlab = "Stand Count")
dev.off()
png("Figures/Plot_Suscep_by_Region.png")
plot(pheno$Susceptible ~ pheno$Region, main = "Susceptibility by Region of Origin")
dev.off()
png("Figures/Plot_Suscep_by_Program.png")
plot(pheno$Susceptible ~ pheno$Program, main = "Susceptibility by Program of Origin")
dev.off()
png("Figures/Plot_Suscep_by_Endo.png")
plot(pheno$Susceptible ~ pheno$EndospermType, main = "Susceptibility by Endosperm Type")
dev.off()
png("Figures/Plot_Endo.png")
plot(pheno$EndospermType,ylab = "Count", main = "Representation of Endosperm Types in Panel")
dev.off()
png("Figures/Plot_Suscp_by_Date.png")
plot(pheno$Susceptible ~ pheno$DatePlanted, main = "Susceptibility by Planting Date")
dev.off()

#how many are resistant?
length(which(pheno$Susceptible == 0)) #47 not susceptible, or 47 resistant
resistantvar <- pheno$Line[which(pheno$Susceptible == 0 )]
print(resistantvar)
write.table(resistantvar, "resistant_varieties.txt", append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = TRUE)

pheno$EndospermType[which(pheno$Susceptible == 0 )]
resistantendo <- table(pheno$EndospermType[which(pheno$Susceptible == 0)])
print(resistantendo)
write.table(resistantendo,"resistant_endosperm_types.txt", append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = TRUE)

endocount <- table(pheno$EndospermType)
print(endocount)
write.table(endocount,"total_endosperm_types.txt", append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = TRUE)


# #########################
# ### GWASPoly, Full Group, No removes ###
# #########################
# #Start with the full group, nothing removed
# GWASPolyRunVersion <- "FullGroup_420lines_PC3"
# 
# #read in hmp data, output numeric format that GWASPoly likes
# SCMVPanel_nwithpos <- hmpToNumeric(geno_scmv)
# outfiles <- WritePhenoGenoToFile(GWASPolyRunVersion,trait,phenoSubsetGeno,SCMVPanel_nwithpos,filename,adendum)
# 
# #now we run the GWASpoly with the files in the proper format
# data <- read.GWASpoly(ploidy=2, 
#                       pheno.file = outfiles[1],
#                       geno.file=outfiles[2],
#                       format="numeric",
#                       n.traits=1,
#                       delim=",")
# 
# data2 <- set.K(data)
# params <- set.params(fixed=NULL, fixed.type=NULL,n.PC = 3, MAF = 0.025) #no fixed effects, MAF should do nothing as it's already been filtered
# 
# data3 <- GWASpoly(data2,models=c("additive","general","1-dom-alt","1-dom-ref"),traits=trait, params=params)
# 
# #visualize the gwas results
# GWASPolyVis(GWASPolyRunVersion, trait, data3)
# 
# #########################
# ### GWASPoly, Full Group, Remove Bad Stands ###
# #########################
# # GWASPolyRunVersion <- "BadStandsRemoved_420lines_PC3"
# # 
# # #subset the pheno data
# # phenoSubsetGeno_removes <- phenoSubsetGeno[!phenoSubsetGeno$Remove == "Remove",]
# # 
# # #remove the things that are labeled "Remove"
# # toRemove <- as.character(phenoSubsetGeno$GenoName[phenoSubsetGeno$Remove == "Remove"])
# # intersect(colnames(geno_scmv), toRemove)
# # geno_scmv_removes = dplyr::select(geno_scmv, -toRemove)
# # dim(geno_scmv)
# # dim(geno_scmv_removes)
# # 
# # 
# # #read in hmp data, output numeric format that GWASPoly likes
# # SCMVPanel_nwithpos <- hmpToNumeric(geno_scmv_removes)
# # 
# # outfiles <- WritePhenoGenoToFile(GWASPolyRunVersion,trait,phenoSubsetGeno_removes,SCMVPanel_nwithpos,filename,adendum)
# # 
# # 
# # #now we run the GWASpoly with the files in the proper format
# # data <- read.GWASpoly(ploidy=2, 
# #                       pheno.file = outfiles[1],
# #                       geno.file=outfiles[2],
# #                       format="numeric",
# #                       n.traits=1,
# #                       delim=",")
# # 
# # data2 <- set.K(data)
# # params <- set.params(fixed=NULL, fixed.type=NULL,n.PC = 3, MAF = 0.025)
# # 
# # 
# # data3 <- GWASpoly(data2,models=c("additive","general","1-dom-alt","1-dom-ref","diplo-general"),traits="Susceptible", params=params)
# # 
# # #visualize the gwas results
# # GWASPolyVis(GWASPolyRunVersion, trait, data3)
# 
# #########################
# ### GWASPoly, Full Group, Remove Bad Stands, 2&3 compressed, 1 removed###
# #########################
# GWASPolyRunVersion <- "BadStandsRemoved_Binary01_420lines_PC3"
# 
# #make the necessary phenotype/genotype files, required by the GWASPoly manual
# phenoAlpha_ss_removes <- phenoAlpha_ss[!phenoAlpha_ss$Remove == "Remove",]
# phenoAlpha_ss_removes_023 <- phenoAlpha_ss[!phenoAlpha_ss_removes$Susceptible == "1",]
# 
# #make all the susceptibility scores that are equal to 2 or 3 into 1s
# #this will make the measurement binary. That which is susceptible (1) and that which is resistant (0)
# phenoAlpha_ss_removes_023$Susceptible[which(phenoAlpha_ss_removes_023$Susceptible == "2")] = 3
# phenoAlpha_ss_removes_023$Susceptible[which(phenoAlpha_ss_removes_023$Susceptible == "3")] = 1
# 
# phenoAlpha_ss_removes_01 <- phenoAlpha_ss_removes_023
# phenoAlpha_ss_removes_01_fullpheno <- phenoAlpha_ss_removes_01[,c('GenoName','Susceptible', "EndospermType", "Region", "Program")]
# phenoAlpha_ss_removes_01_trait_Only <- phenoAlpha_ss_removes_01[,c('GenoName',trait)]
# 
# #read out the phenotype info into a format that GWASpoly will like
# phenooutfilefull <- paste("Data/SCMV_Panel_GWASpoly_Pheno_",GWASPolyRunVersion ,".csv",sep = "")
# write.table(phenoAlpha_ss_removes_01_fullpheno,
#             append = FALSE,
#             file = phenooutfilefull,
#             sep = ",",
#             dec = ".",
#             row.names = FALSE,
#             col.names = TRUE)
# 
# #read out the phenotype info into a format that GWASpoly will like. This time with only susceptibility info incase the other info is confusing
# phenooutfiletrait <- paste("Data/SCMV_Panel_GWASpoly_Pheno_only_",trait,"_",GWASPolyRunVersion ,".csv",sep = "")
# write.table(phenoAlpha_ss_removes_01_trait_Only,
#             append = FALSE,
#             file = phenooutfiletrait,
#             sep = ",",
#             dec = ".",
#             row.names = FALSE,
#             col.names = TRUE)
# 
# 
# #remove the things that are labeled "Remove"
# toRemove <- c(as.character(phenoAlpha_ss$GenoName[phenoAlpha_ss$Remove == "Remove"]),
#               as.character(phenoAlpha_ss$GenoName[phenoAlpha_ss$Susceptible == "1"]))
# geno_scmv_removes <- dplyr::select(geno_scmv, -c(toRemove))
# str(geno_scmv_removes)
# 
# #read in hmp data, output numeric format that GWASPoly likes
# SCMVPanel_nwithpos <- hmpToNumeric(geno_scmv_removes)
# 
# #export as CSV
# genooutput <- paste("Data/SequenceData/",filename,adendum,"_",GWASPolyRunVersion,"_numericFormat.csv",sep = "")
# write.table(SCMVPanel_nwithpos,
#             append = FALSE,
#             file = genooutput,
#             sep = ",",
#             dec = ".",
#             row.names = FALSE,
#             col.names = TRUE)
# 
# 
# #now we run the GWASpoly with the files in the proper format
# data <- read.GWASpoly(ploidy=2, 
#                       pheno.file = phenooutfilefull,
#                       geno.file=genooutput,
#                       format="numeric",
#                       n.traits=1,
#                       delim=",")
# 
# data2 <- set.K(data)
# # params <- set.params(fixed=NULL, fixed.type=NULL)
# params <- set.params(fixed=NULL, fixed.type=NULL,n.PC = 3, MAF = 0.025)
# 
# #fixed = defines the covariates in the pheno file (by column name)
# # fixed.type = numeric (covariates) or factor
# # in this case, we have none
# 
# data3 <- GWASpoly(data2,models=c("additive","general","1-dom-alt","1-dom-ref"),traits="Susceptible", params=params)
# 
# #visualize the gwas results
# GWASPolyVis(GWASPolyRunVersion, trait, data3)
# 
# #########################
# ### GWASPoly, Full Group, Remove Bad Stands,SCMV1 as fixed effect ###
# #########################
# GWASPolyRunVersion <- "BadStandsRemoved_420lines_PC3_SCMV1Fixed"
# 
# phenoAlpha_ss_removes <- phenoAlpha_ss[!phenoAlpha_ss$Remove == "Remove",]
# phenoAlpha_ss_removes_fullpheno <- phenoAlpha_ss_removes[,c('GenoName','Susceptible', "EndospermType", "Region", "Program")]
# phenoAlpha_ss_removes_trait_Only <- phenoAlpha_ss_removes[,c('GenoName',trait)]
# 
# #remove the things that are labeled "Remove"
# toRemove <- as.character(phenoAlpha_ss$GenoName[phenoAlpha_ss$Remove == "Remove"])
# intersect(colnames(geno_scmv), toRemove)
# geno_scmv_removes = dplyr::select(geno_scmv, -toRemove)
# dim(geno_scmv)
# dim(geno_scmv_removes)
# 
# #read in hmp data, output numeric format that GWASPoly likes
# SCMVPanel_nwithpos <- hmpToNumeric(geno_scmv_removes)
# 
# #figure out what the haplotype is at the SCMV1 gene in IL793a (it has Pa405 as parent)
# SCMVPanel_nwithpos_chr6 <- SCMVPanel_nwithpos[SCMVPanel_nwithpos$Chrom == "6"]
# #know the PAV for SCMV1 is around 19.4Mb. Find the closest postion to that
# SCMV1Pos <- as.integer(which(SCMVPanel_nwithpos_chr6$Position == 19369936))
# IL793apos <- as.integer(which(colnames(SCMVPanel_nwithpos_chr6) == "IL793a"))
# SCMV1status <- SCMVPanel_nwithpos_chr6[60,183]
# 
# withSCMV1pos <- which(SCMVPanel_nwithpos_chr6[60,] == as.integer(SCMV1status))
# withSCMV1names <- colnames(SCMVPanel_nwithpos_chr6)[withSCMV1pos]
# 
# phenoAlpha_ss_removes_fullpheno$SCMV1 = FALSE
# phenoAlpha_ss_removes_trait_Only$SCMV1 = FALSE
# phenoAlpha_ss_removes_fullpheno$SCMV1[phenoAlpha_ss_removes_fullpheno$GenoName %in% c(withSCMV1names)]=TRUE
# phenoAlpha_ss_removes_trait_Only$SCMV1[phenoAlpha_ss_removes_trait_Only$GenoName %in% c(withSCMV1names)]=TRUE
# 
# #read out the phenotype info into a format that GWASpoly will like
# phenooutfilefull <- paste("Data/SCMV_Panel_GWASpoly_Pheno_",GWASPolyRunVersion ,".csv",sep = "")
# write.table(phenoAlpha_ss_removes_fullpheno,
#             append = FALSE,
#             file = phenooutfilefull,
#             sep = ",",
#             dec = ".",
#             row.names = FALSE,
#             col.names = TRUE)
# 
# #read out the phenotype info into a format that GWASpoly will like. This time with only susceptibility info incase the other info is confusing
# phenooutfiletrait <- paste("Data/SCMV_Panel_GWASpoly_Pheno_only_",trait,"_",GWASPolyRunVersion ,".csv",sep = "")
# write.table(phenoAlpha_ss_removes_trait_Only,
#             append = FALSE,
#             file = phenooutfiletrait,
#             sep = ",",
#             dec = ".",
#             row.names = FALSE,
#             col.names = TRUE)
# 
# 
# #export as CSV
# genooutput <- paste("Data/SequenceData/",filename,adendum,"_",GWASPolyRunVersion,"_numericFormat.csv",sep = "")
# write.table(SCMVPanel_nwithpos,
#             append = FALSE,
#             file = genooutput,
#             sep = ",",
#             dec = ".",
#             row.names = FALSE,
#             col.names = TRUE)
# 
# 
# #now we run the GWASpoly with the files in the proper format
# data <- read.GWASpoly(ploidy=2, 
#                       pheno.file = phenooutfilefull,
#                       geno.file=genooutput,
#                       format="numeric",
#                       n.traits=1,
#                       delim=",")
# 
# data2 <- set.K(data)
# 
# #set presence of SCMV1 as a fixed effect
# para <- set.params(fixed=c("SCMV1"), fixed.type=c("factor"),n.PC = 3, MAF = 0.025)
# data3 <- GWASpoly(data2,models=c("additive","general","1-dom-alt","1-dom-ref"),traits=c("Susceptible"), params=para)
# 
# #visualize the gwas results
# GWASPolyVis(GWASPolyRunVersion, trait, data3)
# #########################
# ### GWASPoly, Only those with Scmv1, Remove Bad Stands,SCMV1 as fixed effect ###
# #########################
# GWASPolyRunVersion <- "BadStandsSCMV1Removed_420lines_PC3"
# 
# phenoAlpha_ss_removes <- phenoAlpha_ss[!phenoAlpha_ss$Remove == "Remove",]
# phenoAlpha_ss_removes_fullpheno <- phenoAlpha_ss_removes[,c('GenoName','Susceptible', "EndospermType", "Region", "Program")]
# phenoAlpha_ss_removes_trait_Only <- phenoAlpha_ss_removes[,c('GenoName',trait)]
# 
# #remove the things that are labeled "Remove"
# toRemove <- as.character(phenoAlpha_ss$GenoName[phenoAlpha_ss$Remove == "Remove"])
# intersect(colnames(geno_scmv), toRemove)
# geno_scmv_removes = dplyr::select(geno_scmv, -toRemove)
# dim(geno_scmv)
# dim(geno_scmv_removes)
# 
# #read in hmp data, output numeric format that GWASPoly likes
# SCMVPanel_nwithpos <- hmpToNumeric(geno_scmv_removes)
# 
# #figure out what the haplotype is at the SCMV1 gene in IL793a (it has Pa405 as parent)
# SCMVPanel_nwithpos_chr6 <- SCMVPanel_nwithpos[SCMVPanel_nwithpos$Chrom == "6"]
# #know the PAV for SCMV1 is around 19.4Mb. Find the closest postion to that
# SCMV1Pos <- as.integer(which(SCMVPanel_nwithpos_chr6$Position == 19369936))
# IL793apos <- as.integer(which(colnames(SCMVPanel_nwithpos_chr6) == "IL793a"))
# SCMV1status <- SCMVPanel_nwithpos_chr6[60,183]
# 
# withSCMV1pos <- which(SCMVPanel_nwithpos_chr6[60,] == as.integer(SCMV1status))
# withSCMV1names <- colnames(SCMVPanel_nwithpos_chr6)[withSCMV1pos]
# 
# phenoAlpha_ss_removes_fullpheno$SCMV1 = FALSE
# phenoAlpha_ss_removes_trait_Only$SCMV1 = FALSE
# phenoAlpha_ss_removes_fullpheno$SCMV1[phenoAlpha_ss_removes_fullpheno$GenoName %in% c(withSCMV1names)]=TRUE
# phenoAlpha_ss_removes_trait_Only$SCMV1[phenoAlpha_ss_removes_trait_Only$GenoName %in% c(withSCMV1names)]=TRUE
# 
# #remove all those without SCMV1
# #remove the things that are labeled "Remove"
# toRemove <- as.character(phenoAlpha_ss$GenoName[!phenoAlpha_ss_removes_fullpheno$SCMV1])
# intersect(colnames(geno_scmv), toRemove)
# geno_scmv_removes = dplyr::select(geno_scmv, -toRemove)
# dim(geno_scmv)
# dim(geno_scmv_removes)
# 
# 
# #read in hmp data, output numeric format that GWASPoly likes
# SCMVPanel_nwithpos <- hmpToNumeric(geno_scmv_removes)
# 
# #read out the phenotype info into a format that GWASpoly will like
# phenooutfilefull <- paste("Data/SCMV_Panel_GWASpoly_Pheno_",GWASPolyRunVersion ,".csv",sep = "")
# write.table(phenoAlpha_ss_removes_fullpheno,
#             append = FALSE,
#             file = phenooutfilefull,
#             sep = ",",
#             dec = ".",
#             row.names = FALSE,
#             col.names = TRUE)
# 
# #read out the phenotype info into a format that GWASpoly will like. This time with only susceptibility info incase the other info is confusing
# phenooutfiletrait <- paste("Data/SCMV_Panel_GWASpoly_Pheno_only_",trait,"_",GWASPolyRunVersion ,".csv",sep = "")
# write.table(phenoAlpha_ss_removes_trait_Only,
#             append = FALSE,
#             file = phenooutfiletrait,
#             sep = ",",
#             dec = ".",
#             row.names = FALSE,
#             col.names = TRUE)
# 
# 
# #export as CSV
# genooutput <- paste("Data/SequenceData/",filename,adendum,"_",GWASPolyRunVersion,"_numericFormat.csv",sep = "")
# write.table(SCMVPanel_nwithpos,
#             append = FALSE,
#             file = genooutput,
#             sep = ",",
#             dec = ".",
#             row.names = FALSE,
#             col.names = TRUE)
# 
# 
# #now we run the GWASpoly with the files in the proper format
# data <- read.GWASpoly(ploidy=2, 
#                       pheno.file = phenooutfilefull,
#                       geno.file=genooutput,
#                       format="numeric",
#                       n.traits=1,
#                       delim=",")
# 
# data2 <- set.K(data)
# 
# #set presence of SCMV1 as a fixed effect
# para <-  set.params(fixed=NULL, fixed.type=NULL,n.PC = 3, MAF = 0.025)
# data3 <- GWASpoly(data2,models=c("additive","general","1-dom-alt","1-dom-ref"),traits=c("Susceptible"), params=para)
# 
# #visualize the gwas results
# GWASPolyVis(GWASPolyRunVersion, trait, data3)
# 

#########################
### rrBlup version ###
#########################
GWAS(pheno = phenoAlpha_ss_removes_trait_Only,
     geno=SCMVPanel_nwithpos,
     K = NULL,
     n.PC = 3,
     min.MAF = 0.025)

# #########################
# ### Linear Model ###
# #########################
# 
# phenoAlpha_ss_removes$SCMV1 <- phenoAlpha_ss_removes_fullpheno$SCMV1
# 
# formula <- "Susceptible~SCMV1+DatePlanted+Program+StandCountCleaned"
# fit1 <- lm(formula,data=phenoAlpha_ss_removes)
# AIC1 <- extractAIC(fit1)
# AIC1
# summary(fit1)$r.square
# anova(fit1)
# 
# formula2 <- "Susceptible~SCMV1+DatePlanted+Program+BadStand+StandCountCleaned"
# fit2 <- lm(formula2,data=phenoAlpha_ss_removes)
# AIC2 <- extractAIC(fit2)
# AIC2
# summary(fit2)$r.square
# anova(fit2)
# 
# formula3 <- "Susceptible~SCMV1+DatePlanted+Program+EndospermType+BadStand+StandCountCleaned"
# fit3 <- lm(formula3,data=phenoAlpha_ss_removes)
# AIC3 <- extractAIC(fit3)
# AIC3
# summary(fit3)$r.square
# anova(fit3)
# 
# formula4 <- "Susceptible~SCMV1+DatePlanted+Program"
# fit4 <- lm(formula4,data=phenoAlpha_ss_removes)
# AIC4 <- extractAIC(fit4)
# AIC4
# summary(fit4)$r.square
# anova(fit4)
# 
# formula5 <- "Susceptible~SCMV1+DatePlanted+Program+BadStand"
# fit5 <- lm(formula5,data=phenoAlpha_ss_removes)
# AIC5 <- extractAIC(fit5)
# AIC5
# summary(fit5)$r.square
# anova(fit5)
# 
# 
# #what if we got rid of the late planting
# formula6 <- "Susceptible~SCMV1+DatePlanted+Program"
# fit6 <- lm(formula6,data=phenoAlpha_ss_removes[which(phenoAlpha_ss_removes$DatePlanted != "5/22/2019"),])
# AIC6 <- extractAIC(fit4)
# AIC6
# summary(fit6)$r.square
# anova(fit4)
# 
# formula7 <- "Susceptible~Program + Line"
# fit7 <- lm(formula7,data=phenoAlpha_ss_removes)
# AIC7 <- extractAIC(fit7)
# AIC7
# summary(fit7)$r.square
# anova(fit7)
# 


#########################
### LD Visualization ###
#########################
#this is a funtion to calculate what proportion of the compared snps at a given distance apart are in what LD
get.percentiles=function(x, probs=0.95, log_transform=log_transform, numbins=numbins, removefirst500 = removefirst500){
  # 	cat("Numbins is",numbins,"and range of distances is",range(x$dist),"\n")
  mydists = x$dist
  if(removefirst500 == TRUE){
    x <- x[-which(mydists < 501),]
    mydists <- mydists[-which(mydists < 501)]}
  if(log_transform == TRUE){
    is_negative = mydists<0 & !is.na(mydists)
    mydists = abs(mydists)	# Have to abs-value b/c log doesn"t do negatives
    print(range(mydists))
    mydists[mydists==0]=1 # To avoid problems with log-transforming 0-values
    mydists = log10(mydists)	
    mydists[is_negative] = mydists[is_negative] * -1	# Restore to negative to indicate is before the gene
  }
  results=as.list(tapply(X=x$rsq, INDEX=cut(mydists, breaks=numbins), FUN=quantile, na.rm=T, probs=probs))
  results=do.call(what=rbind, args=results)
  return(results)
}

#this is a funtion to calculate what's the median LD the compared SNPS
get.median=function(x, probs=5, log_transform=log_transform, numbins=numbins, removefirst500 = removefirst500){
  mydists = x$dist
  if(removefirst500 == TRUE){
    x <- x[-which(mydists < 501),]
    mydists <- mydists[-which(mydists < 501)]}
  if(log_transform == TRUE){
    is_negative = mydists<0 & !is.na(mydists) 
    mydists = abs(mydists)	# Have to abs-value b/c log doesn"t do negatives
    print(range(mydists))
    mydists[mydists==0]=1 # To avoid problems with log-transforming 0-values
    mydists = log10(mydists)	
    mydists[is_negative] <- mydists[is_negative] * -1	# Restore to negative to indicate is before the gene
  }
  results=as.list(tapply(X=x$rsq, INDEX=cut(mydists, breaks=numbins), FUN=median, na.rm=T))
  results=do.call(what=rbind, args=results)
  return(results)
}

#Import the name of the file. This matrix was calculated on Tassel using the Linkage Disequilibrium tool and the output was saved as a text file
  #since the resulting text file was 8GB, it was pruned using AWK on linux to only include the columns for distance and r2
# LDMatrixfileName <- "WiDivMerged+all_sweet_v120160307_rawcalls_min526_maf0025_maxhet001_mincall040_LD_FullMatrix_r2_distonly"
LDMatrixfileName <- "WiDivMerged+all_sweet_v120160307_rawcalls_min526_maf0025_LD_FullMatrix_locidistr2only"
fullpath <- paste("Data/",LDMatrixfileName,".txt",sep = "")

#Read in the .txt LD matric file
readin <- fread(fullpath, sep=' ')

#Uncomment if all the chromosomes want to be visualized individually, or just the whole geno.
  #based on the numbers of high quality markers I have (only 12k) doing each chromosome is a jittery mess
# chrmiterant <- c("All",6)
chrmiterant <- c("All")
for(j in 1:length(chrmiterant)){
  chrmnum = chrmiterant[j]
  data <- readin
  
  #Look at only one chromosome at a time
  if(chrmnum != "All"){
  data <- data[data$Locus1 == chrmnum]}
  
  #uncomment these if the Matrix file has more than two columns. Eliminate any columns not r2 or dist
  # data <- data[,c(13,14)]
  data <- data[,c(3,4)]
  
  #seperate the distances between chromosomes into 20 bins
  numbins = 20
  #change this if wanting to not log transfrom the data
  log_transform = TRUE
  absolute_distance = TRUE
  
  #change this if you want to eliminate any comparisons less than 500bp apart
  removefirst500 = FALSE
  names(data)=c("dist","rsq")
  data$dist <-as.numeric(data$dist)
  medians=get.median(data, log_transform=log_transform, numbins = numbins, removefirst500 = removefirst500)
  write.table(medians, file="medians.txt",quote=FALSE,sep="\t", col.names=TRUE, row.names=T)
  
  probs=c(0.5, 0.6, 0.7, 0.8, 0.9)
  percents=get.percentiles(data, probs=probs, log_transform=log_transform, numbins = numbins, removefirst500 = removefirst500)
  write.table(percents, file="percents.txt",quote=FALSE,sep="\t", col.names=TRUE, row.names=T)
  mean=tapply(X=data$rsq, INDEX=cut(data$dist, breaks=numbins), FUN=mean, na.rm=T)
  
  percents=read.table('percents.txt',sep='\t',header=T,row.names=1)
  medians=read.table('medians.txt',sep='\t',header=T,row.names=1)
  probs=c('50%','60%','70%','80%','90%')
  colnames(percents) <- probs
  ## Add all together
  row.names(percents)[20]
  
  #the rownames are the starting and end distance of each bin. convert that into a more usable format
  xvals = as.numeric(sub(rownames(percents), pattern="\\((.+),.+", repl="\\1"))
  xvals = c(xvals,8.49)
  mean.xvals <- NULL
  
  #find the middle between each bin border. This will be used for the Xlabel
  for (i in 1:c(length(xvals)-1)){
    mean.xvals <- c(mean.xvals,mean(c(xvals[i],xvals[i+1])))
  }
  
  #determine the file name to save the visualization
  if(removefirst500 == TRUE){
    outfile = paste("Figures/LD_20bins_less500bprmv_from",LDMatrixfileName,"_Chrm",chrmnum,".pdf",sep = "")}
  if(removefirst500 == FALSE){
    outfile = paste("Figures/LD_20bins_from",LDMatrixfileName,"_Chrm",chrmnum,".pdf",sep = "")
  }
  
  #plot and save the visualization
  {pdf(outfile, pointsize=8,family='serif',width = 5,height = 3.5)
    lwd=1
    ylim=c(0,1)
    xlim= c(0, range(mean.xvals)[2] * 1.1)	# Increase the upper xlim by a little bit
    par(cex=1,mar=c(5,5,1,1))
    if(log_transform == TRUE){
      xlabel = "Physical distance (log scale)"
    }else{
      xlabel = "Distance"
    }
    plot(NA, xlim=xlim, ylim=ylim, xlab=xlabel, ylab=bquote("Linkage disequilibrium"~(italic(r)^2)), xaxt="n")	#CHD added
    if(log_transform == TRUE){
      x_ticks = 0:8
      x_labels=c("1 bp","10 bp","0.1 kb","1 kb","10 kb","100 kb","1 Mb","10 Mb","100 Mb") #Jason original
    } else{
      x_ticks = seq(from=xlim[1], to=xlim[2], length.out=5)
      x_labels = x_ticks
    }
    axis(side=1, at=x_ticks, labels=x_labels)
    for(j in 1:length(probs)){
      lines(x=mean.xvals, y=percents[,j], lwd=lwd, col='black', lty=j)
    }
    text(x=sort(xvals)[8], y=c(0.11,0.30,0.46,0.68,0.9), labels=colnames(percents), pos=4, cex = 1)
    dev.off()
  }
}



############
##Phylogenetic Tree Visualization
###########
#tree file made in tassel using neighbor joining on the file outputed at the end of the SNPRelate Step
treeObj <- read.tree(paste(outfile,".nwk",sep=""))

#remove colons and stuff
treeObj_taxa_withcolons <- treeObj$tip.label
treeObj_taxa_wocolons<-gsub(treeObj_taxa_withcolons, pattern = ":.*", replacement = "") #remove the colons and excess characters from tassel

#reorder phenotype info to match the tree
phenoTreeOrder<-phenoAlpha_ss_un[match(treeObj_taxa_wocolons,phenoAlpha_ss_un$GenoName), ]
phenoTreeOrder$VarietyName<- phenoTreeOrder[,1]

#rename the tips of the tree to have the variety names, no the Geno names. I.e. no PI and AMES codes
treeObj$tip.label <- as.character(phenoTreeOrder$VarietyName)

#make sure all the lables agree and are in the same order
head(phenoTreeOrder$GenoName)
head(treeObj_taxa_wocolons)
head(treeObj_taxa_withcolons)
head(treeObj$tip.label)
tail(phenoTreeOrder$GenoName)
tail(treeObj_taxa_wocolons)
tail(treeObj_taxa_withcolons)
tail(treeObj$tip.label)


names <- as.matrix(phenoTreeOrder$VarietyName)[,1]
suscpTrait <- as.matrix(as.numeric(phenoTreeOrder$Susceptible))[,1]
endoTrait <- as.matrix(as.numeric(phenoTreeOrder$EndospermType))[,1]
names(suscpTrait) <- names
suscpTrait

fit <- phytools::fastAnc(treeObj, suscpTrait,vars=TRUE,CI=TRUE)
td <- data.frame(node = nodeid(treeObj, names),
                 Susceptible = suscpTrait)

nd <- data.frame(node = names(fit$ace), Susceptible = fit$ace)
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
treeSCMV <- full_join(treeObj, d, by = 'node')


#plot and export as files
plotname1 <- paste("Figures/Plot_",filename,adendum,"_Phylo_Tree_Suceptibility_ggtree_large.png",sep="")
png(plotname1, width = 3500, height = 3000)
ggtree(treeSCMV, aes(color=Susceptible), layout = 'circular', 
       ladderize = TRUE, continuous = TRUE, size=0.8) + 
  scale_color_gradientn(colours=c("red", 'orange', 'green', 'blue')) +
  geom_tiplab2(hjust = -.1, align = FALSE)  + theme(legend.position = c(.05, .85)) 
dev.off()

plotname2 <- paste("Figures/Plot_",filename,adendum,"_Phylo_Tree_Suceptibility_ggtree_nobranchlength.png",sep="")
png(plotname2, width = 2500, height = 2000)
ggtree(treeSCMV, aes(color=Susceptible), layout = 'circular', 
       ladderize = TRUE, continuous = TRUE, size=0.8,branch.length = 'none') + 
  scale_color_gradientn(colours=c("red", 'orange', 'green', 'blue')) +
  geom_tiplab2(hjust = -.1, align = FALSE)  + theme(legend.position = c(.0, .90))
dev.off()

plotname3 <- paste("Figures/Plot_",filename,adendum,"_Phylo_Tree_Suceptibility_ggtree_large.tiff",sep="")
tiff(plotname3, width = 3500, height = 3000)
ggtree(treeSCMV, aes(color=Susceptible), layout = 'circular', 
       ladderize = TRUE, continuous = TRUE, size=0.8) + 
  scale_color_gradientn(colours=c("red", 'orange', 'green', 'blue')) +
  geom_tiplab2(hjust = -.1, align = FALSE)  + theme(legend.position = c(.05, .85)) 
dev.off()

plotname4 <- paste("Figures/Plot_",filename,adendum,"Phylo_Tree_Suceptibility_ggtree_nobranchlength.tiff",sep="")
tiff(plotname4, width = 2500, height = 2000)
ggtree(treeSCMV, aes(color=Susceptible), layout = 'circular', 
       ladderize = TRUE, continuous = TRUE, size=0.8,branch.length = 'none') + 
  scale_color_gradientn(colours=c("red", 'orange', 'green', 'blue')) +
  geom_tiplab2(hjust = -.1, align = FALSE)  + theme(legend.position = c(.0, .90))
dev.off()

plotname5 <- paste("Figures/Plot_",filename,adendum,"Phylo_Tree_Suceptibility_ggtree_large.bmp",sep="")
bmp(plotname5, width = 3500, height = 3000)
ggtree(treeSCMV, aes(color=Susceptible), layout = 'circular', 
       ladderize = TRUE, continuous = TRUE, size=0.8) + 
  scale_color_gradientn(colours=c("red", 'orange', 'green', 'blue')) +
  geom_tiplab2(hjust = -.1, align = FALSE)  + theme(legend.position = c(.05, .85)) 
dev.off()


#plot with endosperm type not susceptibility
endoTrait[which(is.na(endoTrait))] <-0 
fit <- phytools::fastAnc(treeObj, endoTrait,vars=TRUE,CI=TRUE)
td <- data.frame(node = nodeid(treeObj, names),
                 EndospermType = endoTrait)

nd <- data.frame(node = names(fit$ace), EndospermType = fit$ace)
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
treeSCMV <- full_join(treeObj, d, by = 'node')


plotname6 <- paste("Figures/Plot_",filename,adendum,"Phylo_Tree_endosperm_ggtree_large.png",sep="")
png(plotname6, width = 3500, height = 3000)
ggtree(treeSCMV, aes(color=EndospermType), layout = 'circular', 
       ladderize = TRUE, continuous = TRUE, size=0.8) + 
  scale_color_gradientn(colours=c("black","red", "yellow",'orange', 'green', 'blue', 'purple')) +
  geom_tiplab2(hjust = -.1, align = FALSE)  + theme(legend.position = c(.05, .85)) 
dev.off()

plotname7 <- paste("Figures/Plot_",filename,adendum,"Phylo_Tree_endosperm_ggtree_nobranchlength.png",sep="")
png(plotname7, width = 2500, height = 2000)
ggtree(treeSCMV, aes(color=EndospermType), layout = 'circular', 
       ladderize = TRUE, continuous = TRUE, size=0.8,branch.length = 'none') + 
  scale_color_gradientn(colours=c("red", 'orange', 'green', 'blue')) +
  geom_tiplab2(hjust = -.1, align = FALSE)  + theme(legend.position = c(.0, .90))
dev.off()
#

plotname8 <- paste("Figures/Plot_",filename,adendum,"Phylo_Tree_endosperm_ggtree_large.tiff",sep="")
tiff(plotname8, width = 3500, height = 3000)
ggtree(treeSCMV, aes(color=EndospermType), layout = 'circular', 
       ladderize = TRUE, continuous = TRUE, size=0.8) + 
  scale_color_gradientn(colours=c("red", 'orange', 'green', 'blue')) +
  geom_tiplab2(hjust = -.1, align = FALSE)  + theme(legend.position = c(.05, .85)) 
dev.off()

plotname9 <- paste("Figures/Plot_",filename,adendum,"Phylo_Tree_endosperm_ggtree_nobranchlength.tiff",sep="")
tiff(plotname9, width = 2500, height = 2000)
ggtree(treeSCMV, aes(color=EndospermType), layout = 'circular', 
       ladderize = TRUE, continuous = TRUE, size=0.8,branch.length = 'none') + 
  scale_color_gradientn(colours=c("red", 'orange', 'green', 'blue')) +
  geom_tiplab2(hjust = -.1, align = FALSE)  + theme(legend.position = c(.0, .90))
dev.off()

plotname10 <- paste("Figures/Plot_",filename,adendum,"Phylo_Tree_endosperm_ggtree_large.bmp",sep="")
bmp(plotname10, width = 3500, height = 3000)
ggtree(treeSCMV, aes(color=EndospermType), layout = 'circular', 
       ladderize = TRUE, continuous = TRUE, size=0.8) + 
  scale_color_gradientn(colours=c("red", 'orange', 'green', 'blue')) +
  geom_tiplab2(hjust = -.1, align = FALSE)  + theme(legend.position = c(.05, .85)) 
dev.off()



#########################
###Control Analysis###
#########################
#how susceptible were the controls? if the controls were totally infected, that means our virus application method was effective
controls <- read.delim("Data/WSMDP_SCMV_Control_Data_Raw.csv", header = TRUE, sep = ",", dec = ".")
controls$Bench <- as.factor(controls$Bench)
controls$Check <- as.factor(controls$Check)
controls$ratio <- controls$Infected/controls$Stand
mean(controls$ratio)