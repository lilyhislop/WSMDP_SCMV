GWASPolyRunner <- function(RunName, trait,geno_scmv,filename,adendum,phenoSubsetGeno){

  GWASPolyRunVersion <- RunName

#read in hmp data, output numeric format that GWASPoly likes
ifelse(!is.numeric(geno_scmv[15,30]),SCMVPanel_nwithpos <- hmpToNumeric(geno_scmv),SCMVPanel_nwithpos <- geno_scmv)
outfiles <- WritePhenoGenoToFile(GWASPolyRunVersion,trait,phenoSubsetGeno,SCMVPanel_nwithpos,filename,adendum)

#now we run the GWASpoly with the files in the proper format
data <- read.GWASpoly(ploidy=2, 
                      pheno.file = outfiles[1],
                      geno.file=outfiles[2],
                      format="numeric",
                      n.traits=4,
                      delim=",")

data2 <- set.K(data, LOCO = FALSE)
params <- set.params(fixed=NULL, fixed.type=NULL,n.PC = 3, MAF = 0.05) #no fixed effects, MAF should do nothing as it's already been filtered

data3 <- GWASpoly(data2,models=c("additive","general"),traits=trait, params=params)

#visualize the gwas results
GWASPolyVis(GWASPolyRunVersion, trait, data3,filename,adendum)
}