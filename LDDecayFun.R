
###############################
# LDFunction
###############################
#write a function to calculate the LD decay for each chromosome
#for genofile, input the gdsobj you'd like to analyze
#for samples, input the sample.ids you like to investigate
#for method, input the method of LD analysis you'd like to use "composite", "r", "dprime", or "corr". default is r

LDDecay <- function(genofile,samples,method="r"){
  for(chrnum in 1:10){
    #which snps are within in the inputted chromosome
    ChrSnpID <- read.gdsn(index.gdsn(genofile,"snp.id"))[read.gdsn(index.gdsn(genofile, "snp.chromosome"))== chrnum]
    #get the data for all snps including position
    SNPdata <- snpgdsSNPList(genofile, sample.id = samples)
    #Get the data of the subset of snps denoted by chromosome
    SNPdatasub <- SNPdata[SNPdata$chromosome == chrnum,]
    print(paste("Chromosome ",chrnum," is ",length(ChrSnpID)," SNPs long"))
    
    #starting position on the chromosome
    Startpos <- SNPdata$position[head(ChrSnpID,n=1)]
    #ending position on the chromosome
    Endpos <- SNPdata$position[tail(ChrSnpID,n=1)]
    #total length of the chromosome
    Totalpos <- Endpos - Startpos
    #determine which snps to start the ld comparisons on. Start 10k away from first position and a million bp away from end
    #looking at 40 different starting location
    startpoint <- seq(Startpos+10000, Endpos-1000003,(Totalpos-1010003)/40)
    # startpoint <- seq(Startpos+10000, Endpos-1000003,(Totalpos-1010003)/5)
    
    #determine how far away the comparison snp should be
    #note! at the marker density used at time of writing, the shortest distances are esp wonky. 100:1000 seems to not work consistently
    spacing <- c(0,100,500,1000,2000,5000,10000,15000,20000,30000,45000,60000,100000,200000,500000,1000000,1500000)
    # spacing <- c(0,500,2000,10000,20000,45000,100000,1000000,1500000)#smaller version for debugging
    
    #establish a matrix where the LD data will go for several different starting snps and several different distances away from the starting snp
    #the n+1 row will be where we put the averages
    #this matrix will store the LD values
    LD <- matrix(nrow = length(startpoint)+1, ncol = length(spacing))
    colnames(LD) <- spacing
    #this matrix will store the records of which snps were analyzed
    LDSNPpos <- matrix(nrow = length(startpoint)+1, ncol = length(spacing))
    
    #forloops to look at all of these combinations of starting points and distances away
    for(i in 1:length(spacing)) #distance away from starting snp
    {
      for(j in 1:length(startpoint)) #starting snp
      {
        #only do this if the start point and the spacing wont go over the edge of the chromosome.
        #No point analying the LD for a snp that doesnt exist
        if((startpoint[j]+spacing[i])< Totalpos){
          
          #look at starting snp
          #this is looking for the snp that is closest to the starting position, since we don't have perfect snp coverage at all positions
          ####This might be were stuff is going weird?########
          snp1id <-SNPdata$snp.id[which.min(abs(SNPdatasub$position - startpoint[j]))]
          #look at snp spacing distance away from starting snp
          snp2id <- SNPdata$snp.id[which.min(abs(SNPdatasub$position - (startpoint[j]+spacing[i])))]
          
          #only bother calculating LD if theres a marker different marker the proper distance away to compre with
          #If the snps are the same because theres not enough density at that point, skip. if snp1==snp2 then LD will be 1
          if(snp1id != snp2id || spacing[i] < 100){
            #get snp1 and snp2 based on what the closes position is
            snp1 <- as.vector(snpgdsGetGeno(genofile, sample.id = samples, snp.id = snp1id))
            snp2 <- as.vector(snpgdsGetGeno(genofile, sample.id = samples, snp.id = snp2id))
            
            #whats the LD between those two snps
            temp <- snpgdsLDpair(snp1,snp2, method = method)
            
            #put those LD into a matrix
            LD[j,i] <- temp[1]
            #record which snps were compared and what distance apart they are supposed to represent
            LDSNPpos[j,i] <- paste(snp1id, ",", snp2id, ifelse(snp1id == snp2id, "SAME",""),"spacing",spacing[i])
          }
        }
      }
      #What the average LD? averaging among all the different starting points
      LD[length(startpoint)+1,i] <- mean(LD[1:length(startpoint),i],na.rm = TRUE)
    }
    
    
    get.percentiles=function(x, probs=0.95, log_transform=log_transform, numbins=numbins){
      # 	cat("Numbins is",numbins,"and range of distances is",range(x$dist),"\n")
      mydists = x$dist
      if(log_transform == TRUE){
        is_negative = mydists<0
        mydists = abs(mydists)	# Have to abs-value b/c log doesn"t do negatives
        print(range(mydists))
        mydists[mydists==0]=1 # To avoid problems with log-transforming 0-values
        mydists = log10(mydists)	
        mydists[is_negative] = mydists[is_negative] * -1	# Restore to negative to indicate is before the gene
      }
      results=tapply(X=x$rsq, INDEX=cut(mydists, breaks=numbins), FUN=quantile, na.rm=T, probs=probs)
      results=do.call(what=rbind, args=results)
      return(results)
    }
    
    get.median=function(x, probs=5, log_transform=log_transform, numbins=numbins){
      mydists = x$dist
      if(log_transform == TRUE){
        is_negative = mydists<0
        mydists = abs(mydists)	# Have to abs-value b/c log doesn"t do negatives
        print(range(mydists))
        mydists[mydists==0]=1 # To avoid problems with log-transforming 0-values
        mydists = log10(mydists)	
        mydists[is_negative] = mydists[is_negative] * -1	# Restore to negative to indicate is before the gene
      }
      results=tapply(X=x$rsq, INDEX=cut(mydists, breaks=numbins), FUN=median, na.rm=T)
      results=do.call(what=rbind, args=results)
      return(results)
    }
    
    
    
    #plot the LD for these parameters. Use a log scale to make the datapoints farther apart and more readable
    #I should probably use GG plot to make this prettier, but not a high priority at the moment
    if(chrnum == 1){
      plot(spacing,LD[nrow(LD),], log = "x", type = "l", col = chrnum, xaxt='n', main = paste("LD Decay using method",method),ylab = "LD", xlab = "Distance" )
      axis(side = 1, at=spacing, labels = as.character(spacing),las = 2)
      # returnlist <- list(LDSNPpos,LD[nrow(LD),],spacing)
    }
    if(chrnum != 1){
      lines(spacing,LD[nrow(LD),], log = "x", type = "l", col = chrnum)
      # returnlist <- list.append(LDSNPpos,LD[nrow(LD),],spacing)
    }
    #return 1 is the record of the snps analyzed
    #return 2 is the average LD at each distance
    #return 3 is the distances analyzed
    
  }
  legend(x=5e05,y=1,legend = c("Chr 1","Chr 2","Chr 3","Chr 4","Chr 5","Chr 6","Chr 7","Chr 8","Chr 9","Chr 10"),col = c(1:10),lty = 1)
  # return(returnlist)
}