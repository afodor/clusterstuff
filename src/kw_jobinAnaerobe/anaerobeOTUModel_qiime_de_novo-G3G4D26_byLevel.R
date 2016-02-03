##OTU model for QIIME de novo anaerobe data
##each variable separately
##Group 3 and 4, day 26 only
##2/1/16

rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")

setwd("C:\\Users\\kwinglee.cb3614tscr32wlt\\Documents\\Fodor\\JobinCollaboration\\C jejuni anaerobe v aerobe\\qiimeDeNovo\\by taxa level")

taxaLevels = c("phylum", "class", "order", "family", "genus")

for(taxa in taxaLevels) {
  paste(taxa)
  
  filename = paste("anaerobe_de_novo_", taxa, "_taxaAsCol_logNorm_G3G4D26.txt", sep="")
  table = read.table(filename, sep="\t", header=T)
  nc = ncol(table)
  cc = c(rep("character", 3), rep("numeric", 4), rep("character", 5), rep("numeric", nc-12))
  table = read.table(filename, sep="\t", header=T, colClasses=cc)
  
  ##output vectors
  names = vector()
  meanTaxa = vector()
  meanGroup3 = vector()
  meanGroup4 = vector()
  pGroup = vector()
  pCage = vector()
  index = 1
  pdf(paste("anaerobeStool_qiimeDeNovo_otuModel_boxplots_", taxa, "_G3G4D26.pdf", sep=""))
  
  for( i in 13:nc) {
    if(sum(table[,i] != 0 ) > nrow(table) / 4 ) {
      
      names[index] <- names(table)[i]
      bug <- table[,i]
      meanTaxa[index] <- mean(bug)
      meanGroup3[index] <- mean(bug[table$Group==3])
      meanGroup4[index] <- mean(bug[table$Group==4])
      group <- factor(table$Group)
      cage <- factor(table$Cage.) #cage and mouse are unique values (no repeats for different mice)
      
      myFrame <- data.frame(bug, group, cage)
      
      pGroup[index] = anova(lm(bug~group))$"Pr(>F)"[1]
      pCage[index] = anova(lm(bug~cage))$"Pr(>F)"[1]
      
      ##plots
      graphMain =  paste( names[index], 
                          ":\n pGroup=", format(pGroup[index], digits=3), 
                          "; pCage= ", format(pCage[index], digits=3), sep="")
      
      
      ##color by group
      col = ifelse(table$Group==3, "blue", "red")
      sh = 16
      
      par(mfrow=c(1,2), oma=c(.5,.5,5,.5), mar=c(4.8,4,1,1))
      ##group
      boxplot(bug~group, xlab="", ylab="relative abundance", main="group")
      points(bug~group, col=col, pch=sh)
      ##cage
      boxplot(bug~factor(cage), xlab="", ylab="relative abundance",  main="cage")
      points(bug~factor(cage), col=col, pch=sh)
      
      ##add title
      par(oma=c(0,0,0,0), mar=c(0,0,0,0), new=T, xpd=T, fig=c(0,1,0,1))
      plot(0,0,type="n", bty="n", xaxt="n", yaxt="n")
      legend("top", horiz=T, legend=graphMain, cex=.7, bty="n")
      
      index=index+1
      
    }
  }
  
  dFrame <- data.frame(names, meanTaxa, meanGroup3, meanGroup4, pGroup, pCage)
  dFrame$pAdjGroup <- p.adjust(dFrame$pGroup, method = "BH")
  dFrame$pAdjCage <- p.adjust(dFrame$pCage, method="BH")
  dFrame <- dFrame[order(dFrame$pGroup),]
  write.table(dFrame, file=paste("anaerobeStool_qiimeDeNovo_otuModel_pValues_", taxa, "_G3G4D26.txt", sep=""), sep="\t",row.names=FALSE, quote=F)
  dev.off()
}