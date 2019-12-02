### Propensity Score

smahal = function(z, X) {
  X <- as.matrix(X)
  n <- dim(X)[1]
  rownames(X) <- 1:n
  k <- dim(X)[2]
  m <- sum(z)
  for (j in 1:k) X[, j] <- rank(X[, j])
  cv <- cov(X)
  vuntied <- var(1:n)
  rat <- sqrt(vuntied/diag(cv))
  cv <- diag(rat) %*% cv %*% diag(rat)
  out <- matrix(NA, m, n - m)
  Xc <- X[z == 0, ]
  Xt <- X[z == 1, ]
  rownames(out) <- rownames(X)[z == 1]
  colnames(out) <- rownames(X)[z == 0]
  library(MASS)
  icov <- ginv(cv)
  for (i in 1:m) out[i, ] <- mahalanobis(Xc, Xt[i, ], icov, inverted = T)
  out
}


addcaliper = function(dmat, z, logitp, calipersd = 0.2, penalty = 1000) {
  sd.logitp = sd(logitp)
  adif = abs(outer(logitp[z == 1], logitp[z == 0], "-"))
  adif = (adif - (calipersd * sd.logitp)) * (adif > (calipersd * sd.logitp))
  dmat = dmat + adif * penalty
}

matching <- function(x, paired=NA){
  flag = sum(is.na(paired))==0
  n = dim(x)[1]
  p = dim(x)[2]
  colnames(x)[p-1] = "Shock"
  index = sapply(1:n, function(i){
    return(sum(sapply((1:p)[-3], function(j){
      return(is.na(x[i,j]))
    })))
  })
  index = which(index==0)
  tmp = tmp1 = x[index,-c(3)]
  propscore.model = glm(Shock ~. -ID, family = binomial, x = TRUE, y = TRUE, data = tmp)
  
  X = propscore.model$x[, -1]
  A = propscore.model$y
  logitps = predict(propscore.model)
  names(logitps) = rownames(tmp)
  distmat = smahal(A, X)
  rownames(distmat) = rownames(tmp)[tmp$Shock == 1]
  colnames(distmat) = rownames(tmp)[tmp$Shock == 0]
  distmat_caliper = addcaliper(distmat, A, logitps, calipersd = 0.5)
  if(flag){
    U = union(paired[,1], paired[,2])
    index = which(rownames(tmp1) %in% U)
    tmp1 = tmp1[-index,]
    distmat_caliper = distmat_caliper[-which(rownames(distmat_caliper) %in% U),
                                      -which(colnames(distmat_caliper) %in% U)]
  }
  tmp1 = tmp1[,-(p-1)]
  
  noControls = 1
  matchvec = pairmatch(distmat_caliper, controls = noControls, data = tmp1)
  
  matchvec.num = as.numeric(substr(matchvec, start = 3, stop = 10))
  matchvec.num.notNA = matchvec.num[!is.na(matchvec.num)] #To remove individuals who didn't get matched
  matchID = unique(matchvec.num.notNA)
  I = length(matchID)
  matchedPairMat = matrix(0, I, 4)
  colnames(matchedPairMat) = c("SubjectID (Treated)", "SubjectID (Control)", "PS (Treated)",
                               "PS (Control)")
  treatedSubjID = rownames(tmp1)[tmp1$Shock == 1] 
  controlSubjID = rownames(tmp1)[tmp1$Shock == 0]
  for (i in 1:I) {
    subjectIDs = which(matchvec.num == matchID[i])
    subjectIDs = as.integer(rownames(tmp1)[subjectIDs])
    matchedPairMat[i, "SubjectID (Treated)"] = subjectIDs[subjectIDs %in% treatedSubjID]
    matchedPairMat[i, "SubjectID (Control)"] = subjectIDs[subjectIDs %in% controlSubjID]
    matchedPairMat[i, "PS (Treated)"] = round(logitps[which(names(logitps)
                                                            == matchedPairMat[i, "SubjectID (Treated)"])], 3)
    matchedPairMat[i, "PS (Control)"] = round(logitps[which(names(logitps)
                                                            == matchedPairMat[i, "SubjectID (Control)"])], 3)
  }
  knitr::kable(head(matchedPairMat), caption = "430 Matched Pairs")
  if(flag)
    return(rbind(paired,matchedPairMat))
  else
    return(matchedPairMat)
}
### create structured data
load("dat_v4.rdata")


system.time({
  I = matching(struc_dat[[1]])
  for(i in 2:12){
    I = matching(struc_dat[[i]], paired = I)
  }
})

dat$Deathflag = sapply(1:dim(dat)[1], function(i){
  if(is.na(dat$RADYEAR[i]))
    return(0)
  else
    return(1)
})

trt = dat[I[,1],]
con = dat[I[,2],]
Deathflag = sapply(1:dim(trt)[1], function(i){
  d1 = 0
  d2 = 0
  tmp = as.numeric(trt[i, 78:89]) - c(rep(1,11),0)
  s = which(tmp>0)[1]
  if(!is.na(trt$RADYEAR[i]) && trt$RADYEAR[i] <= 1994+2*s)
    d1 = 1
  if(!is.na(con$RADYEAR[i]) && con$RADYEAR[i] <= 1994+2*s)
    d2 = 1
  return(as.vector(c(d1,d2)))
})
trt$Deathflag = Deathflag[1,]
con$Deathflag = Deathflag[2,]
t.test(trt$Deathflag,con$Deathflag,paired = T)

### person year
years = as.matrix(dat[,57:68])
years = t(apply(years, 1, function(x){
  return(1-is.na(x))
}))
dat$years = apply(years,1,function(x){
  return(sum(2*x))
})
sum(years)
sum(trt$Deathflag)

###plot
library(tidyverse)
trt$tr = 1
con$tr = 0
plotdat = rbind(trt,con)
levels(plotdat$H6INPOV) = c("NO","YES")
levels(plotdat$H12INPOV) = c("NO","YES")
ggplot(plotdat, aes(factor(R1MSTAT), fill = factor(tr))) +
  geom_bar(position = position_dodge2(preserve = "single"))+
  coord_flip()+
  theme(legend.justification = "bottom")+
  labs(x = "Wave 1 Marital Status")
ggplot(plotdat, aes(factor(R6MSTAT), fill = factor(tr))) +
  geom_bar(position = position_dodge2(preserve = "single"))+
  theme(legend.justification = "bottom")+
  coord_flip()+
  labs(x = "Wave 6 Marital Status")
ggplot(plotdat, aes(factor(R12MSTAT), fill = factor(tr))) +
  geom_bar(position = position_dodge2(preserve = "single"))+
  coord_flip()+
  theme(legend.justification = "bottom")+
  labs(x = "Wave 12 Marital Status")
ggplot(plotdat, aes(factor(RAGENDER), fill = factor(tr))) +
  geom_bar(position = position_dodge2(preserve = "single"))+
  coord_flip()+
  theme(legend.justification = "bottom")+
  labs(x = "Gender")
ggplot(plotdat, aes(factor(RAEDEGRM), fill = factor(tr))) +
  geom_bar(position = position_dodge2(preserve = "single"))+
  coord_flip()+
  theme(legend.justification = "bottom")+
  labs(x = "Education Background")
ggplot(plotdat, aes(factor(RABYEAR), fill = factor(tr))) +
  geom_bar(position = position_dodge2(preserve = "single"))+
  coord_flip()+
  theme(legend.justification = "bottom")+
  labs(x = "Birth Year")
ggplot(plotdat, aes(factor(H6INPOV), fill = factor(tr))) +
  geom_bar(position = position_dodge2(preserve = "single"))+
  coord_flip()+
  theme(legend.justification = "bottom")+
  labs(x = "Wave 6 in Poverty")
ggplot(plotdat, aes(factor(H12INPOV), fill = factor(tr))) +
  geom_bar(position = position_dodge2(preserve = "single"))+
  coord_flip()+
  theme(legend.justification = "bottom")+
  labs(x = "Wave 12 in Poverty")
ggplot(plotdat, aes(factor(R1SHLT), fill = factor(tr))) +
  geom_bar(position = position_dodge2(preserve = "single"))+
  coord_flip()+
  theme(legend.justification = "bottom")+
  labs(x = "Wave 1 Self Report Health")
ggplot(plotdat, aes(factor(R6SHLT), fill = factor(tr))) +
  geom_bar(position = position_dodge2(preserve = "single"))+
  coord_flip()+
  theme(legend.justification = "bottom")+
  labs(x = "Wave 6 Self Report Health")
ggplot(plotdat, aes(factor(R12SHLT), fill = factor(tr))) +
  geom_bar(position = position_dodge2(preserve = "single"))+
  coord_flip()+
  theme(legend.justification = "bottom")+
  labs(x = "Wave 12 Self Report Health")
ggplot(plotdat, aes(factor(R1HLTHLM), fill = factor(tr))) +
  geom_bar(position = position_dodge2(preserve = "single"))+
  coord_flip()+
  theme(legend.justification = "bottom")+
  labs(x = "Wave 1 Health Influence Labor")
ggplot(plotdat, aes(factor(R6HLTHLM), fill = factor(tr))) +
  geom_bar(position = position_dodge2(preserve = "single"))+
  coord_flip()+
  theme(legend.justification = "bottom")+
  labs(x = "Wave 6 Health Influence Labor")
ggplot(plotdat, aes(factor(R12HLTHLM), fill = factor(tr))) +
  geom_bar(position = position_dodge2(preserve = "single"))+
  coord_flip()+
  theme(legend.justification = "bottom")+
  labs(x = "Wave 12 Health Influence Labor")


par(mfrow=c(2,2))
boxplot(con$H1ITOT, trt$H1ITOT,outline = F,ylim = c(0,1e5),names = c("Con","Trt"),ylab = "Wave 1 Income")
boxplot(con$H2ITOT, trt$H2ITOT,outline = F,ylim = c(0,1e5),names = c("Con","Trt"),ylab = "Wave 2 Income")
boxplot(con$H3ITOT, trt$H3ITOT,outline = F,ylim = c(0,1e5),names = c("Con","Trt"),ylab = "Wave 3 Income")
boxplot(con$H4ITOT, trt$H4ITOT,outline = F,ylim = c(0,1e5),names = c("Con","Trt"),ylab = "Wave 4 Income")
par(mfrow=c(2,2))
boxplot(con$H5ITOT, trt$H5ITOT,outline = F,ylim = c(0,1e5),names = c("Con","Trt"),ylab = "Wave 5 Income")
boxplot(con$H6ITOT, trt$H6ITOT,outline = F,ylim = c(0,1e5),names = c("Con","Trt"),ylab = "Wave 6 Income")
boxplot(con$H1ITOT, trt$H7ITOT,outline = F,ylim = c(0,1e5),names = c("Con","Trt"),ylab = "Wave 7 Income")
boxplot(con$H2ITOT, trt$H8ITOT,outline = F,ylim = c(0,1e5),names = c("Con","Trt"),ylab = "Wave 8 Income")
par(mfrow=c(2,2))
boxplot(con$H3ITOT, trt$H9ITOT,outline = F,ylim = c(0,1e5),names = c("Con","Trt"),ylab = "Wave 9 Income")
boxplot(con$H4ITOT, trt$H10ITOT,outline = F,ylim = c(0,1e5),names = c("Con","Trt"),ylab = "Wave 10 Income")
boxplot(con$H5ITOT, trt$H11ITOT,outline = F,ylim = c(0,1e5),names = c("Con","Trt"),ylab = "Wave 11 Income")
boxplot(con$H6ITOT, trt$H12ITOT,outline = F,ylim = c(0,1e5),names = c("Con","Trt"),ylab = "Wave 12 Income")



shock = as.matrix(dat[,78:89])
shock = apply(shock,1, function(x){
  x = as.numeric(x)
  return(as.numeric(sum(x,na.rm = T)>0))
})
sum(shock)
trt = dat[which(shock>0),]
con = dat[which(shock==0),]
t.test(trt$Deathflag,con$Deathflag)
