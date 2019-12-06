library(sensitivitymw)
# senmw does matching with multiple controls (same number of controls)
GammaSeq = seq(1,1.7,0.1)
GammaSeq = seq(1.5,1.6,0.01)
upperBound = matrix(0,length(GammaSeq),3)
tem = cbind(trt$Deathflag, con$Deathflag)
for(i in 1:length(GammaSeq)) {
  upperBound[i,1] =senmw(tem, gamma = GammaSeq[i], method = "t")$pval
  upperBound[i,2] =senmw(tem, gamma = GammaSeq[i], method = "p")$pval
  upperBound[i,3] = senmw(tem, gamma = GammaSeq[i], method = "w")$pval
}
out = cbind(GammaSeq,upperBound)
colnames(out) = c("Gamma","t-test","trimmed mean test","weighted trimmed mean test")
round(out,3)

senmwCI(tem, gamma = 1.59, method = 't', one.sided = TRUE)

library(sensitivitymv)
amplify(1.59, c(4 : 7))
uniroot(function(x){amplify(1.59,x) - x},c(1.59+0.01,10))$root
