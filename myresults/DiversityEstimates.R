setwd("~/UVM/EcologicalGenomics/EcologicalGenomics/myresults")

list.files()
SFS <- scan("XWS_outfold.sfs")
sumSFS <- sum(SFS)
sumSFS
pctPoly = 100*(1-(SFS[1]/sumSFS))
plotSFS <- SFS[-c(1,length(SFS))]
barplot(plotSFS, xlab="XWS Pop SFS")
div <- read.table("XWS_folded_allsites.thetas.idx.pestPG")
colnames(div)=c("window","chrome","wincenter","tW","tP","tF","tH","tL","tajD","fulif","fuliD","fayH","zengsE","numSites")
print(div)
div$tWpersite = div$tW/div$numSites
div$tPpersite = div$tP/div$numSites
tajDEst <- div$tPpersite-div$tWpersite
pdf("XWS_diversity_stats2.pdf")
par(mfrow=c(2,2))

hist(div$tWpersite,col="darkorchid",xlab="Theta_W",main="")
hist(div$tPpersite,col="darkorchid", xlab="Theta-Pi",main="")
hist(tajDEst,col="darkorchid",xlab="Tajima's D",main="")

summary(div)

barplot(plotSFS)
dev.off()

# Alison's experimenting to get all the graphs on one page, not part of the class work ggplot2 -----------------------------------------------------------------

library(cowplot)
p4 <- barplot(plotSFS, xlab="XWS Pop SFS")
p2 <- hist(div$tWpersite,col="darkorchid",xlab="Theta_W",main="")
p3 <- hist(tajDEst,col="darkorchid",xlab="Tajima's D",main="")
p1 <- hist(div$tPpersite,col="darkorchid", xlab="Theta-Pi",main="")
ggarrange(p4,p3,p2,p1,ncol=2,nrow=2)
                  
  