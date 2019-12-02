library(sensitivitymw)
# senmw does matching with multiple controls (same number of controls)
GammaSeq = seq(1,1.7,0.01)
upperBound = matrix(0,length(GammaSeq),1)
tem = cbind(trtgroup$deathflag, ctrlgroup$deathflag)
for(i in 1:length(GammaSeq)) {
  upperBound[i,1] =senmw(tem, gamma = GammaSeq[i], method = "t")$pval
}
out = cbind(GammaSeq,upperBound)
colnames(out) = c("Gamma","t-test")
round(out,3)