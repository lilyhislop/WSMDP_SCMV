PhenowithGenoPrep <- function(phenoFile){
  #read in the phenotype data
  #this pheno data includes the name called by the GBS files, The colloquial variety names,
  #the date planted, the region and program, the stand count, the susceptibility score
  #use "pheno" if you want to consult all the whole phenotypic data set. use phenoAlpha_ss if you want only phenotypic data that has genotypic data associated with it
  pheno <- read.delim(phenoFile, header = TRUE, sep = ",", dec = ".")
  
  
  #Phenotype info cleanup
  #order the phenotype alphabetically by GBS name
  phenoAlpha <- pheno[order(pheno$GenoName),]
  
  #remove anything doesnt have a genoname, so doesnt have sequence data and can't be genetically analyzed
  phenoAlpha_ss <- phenoAlpha[!phenoAlpha$GenoName == "",]
  
  #remove any duplicates so I dont anaylze them twice
  phenoAlpha_ss <- phenoAlpha_ss[!duplicated(phenoAlpha_ss$GenoName),]
  
  phenoAlpha_ss_un <- phenoAlpha_ss[,c(1,4,26,28,29,30,22)]
  
  return(phenoAlpha_ss)}
