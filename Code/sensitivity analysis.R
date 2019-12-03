library(sensitivitymw)
# senmw does matching with multiple controls (same number of controls)
GammaSeq = seq(1,1.7,0.1)
GammaSeq = seq(1.5,1.6,0.01)
upperBound = matrix(0,length(GammaSeq),1)
tem = cbind(trt$Deathflag, con$Deathflag)
for(i in 1:length(GammaSeq)) {
  upperBound[i,1] =senmw(tem, gamma = GammaSeq[i], method = "t")$pval
}
out = cbind(GammaSeq,upperBound)
colnames(out) = c("Gamma","t-test")
round(out,3)

senmwCI(tem, gamma = 1.35, method = 't', one.sided = TRUE)

