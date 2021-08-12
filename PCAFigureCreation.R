
PCAFigureCreation <- function(SCMV_PCA,pc.perc,phenoSubsetGeno,filename, Colourant ){
library(forcats)
progfile <- paste("Figures/Plot_PCA_",filename,"_by_",Colourant,"_test.png",sep = "")

if(Colourant == "Program"){
  #label with program
  phenoSubsetGeno$ProgramOthered <- phenoSubsetGeno$Program
  toOther <- c("FL","GA","MA","ME","MO","ND","NJ","OH","SC","SP","TN","USDA","")
  levels(phenoSubsetGeno$ProgramOthered) <- c(levels(phenoSubsetGeno$ProgramOthered),"Other")
  phenoSubsetGeno$ProgramOthered[phenoSubsetGeno$Program %in% toOther]<- "Other"
  colors <- c("slategray3", "violet", "royalblue4", "palegreen", "mediumpurple", "khaki", "red", "grey")
  symbols <- c(0,1,15,2,3,4,18,6)
  Colourant <- "ProgramOthered"
  x <- factor(phenoSubsetGeno$ProgramOthered)
  x2 <- fct_relevel(x,"Other",after = Inf)
  phenoSubsetGeno$ProgramOthered <- x2
}
colors <- c("slategray3", "violet", "royalblue4", "mediumpurple", "khaki", "grey", "palegreen","red")
symbols <- c(0,1,15,6,22,20,2,16)

if(Colourant == "SusceptibilityRating03"){
  colors <- c("royalblue4", "palegreen", "mediumpurple", "red")
symbols <- c(0,1,15,2)}

png(progfile,width = 750,height = 750)
#label PCA with region
tab <- data.frame(sample.id = SCMV_PCA$sample.id, 
                  pop = factor(phenoSubsetGeno[,Colourant])[match(phenoSubsetGeno$GenoName,SCMV_PCA$sample.id)],
                  EV1 = SCMV_PCA$eigenvect[,1],
                  EV2 = SCMV_PCA$eigenvect[,2],
                  stringsAsFactors = FALSE)
# plot(tab$EV2, tab$EV1,pch = symbols[as.integer(tab$pop)], col = colors[as.integer(tab$pop)], xlab = paste("PC 2 (",(round(pc.perc,2))[2],"%)" ,sep = ""), ylab = paste("PC 1 (",(round(pc.perc,2))[1],"%)" ,sep = ""), main = paste("PCA with", Colourant))
order <- levels(tab$pop)
plot(tab$EV2, tab$EV1,pch = symbols[as.integer(tab$pop)], col = colors[as.integer(tab$pop)], main = "", xlab = paste("PC 2 (",(round(pc.perc,1))[2],"%)" ,sep = ""), ylab = paste("PC 1 (",(round(pc.perc,1))[1],"%)" ,sep = ""))
legend("bottomleft", legend = order , pch = symbols, col = colors)
dev.off()
}



