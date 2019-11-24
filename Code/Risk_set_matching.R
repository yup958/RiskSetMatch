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

### create structured data
struc_dat = list("1","2","3","4","5","6","7","8","9","10","11","12")
struc_dat[[1]] = dat[ , which(colnames(dat)%in%c("ID","R1AGEY_B","RAGENDER",
                                                 "H1ITOT","R1SHLT","R1HLTHLM","SHOCK2",
                                                 "R1MSTAT","RAEDEGRM","RADYEAR"))]
struc_dat[[2]] = dat[ , which(colnames(dat)%in%c("ID","R2AGEY_B","RAGENDER","SHOCK3",
                                                 "H2ITOT","R2SHLT","R2HLTHLM",
                                                 "R2MSTAT","RAEDEGRM","RADYEAR"))]
struc_dat[[3]] = dat[ , which(colnames(dat)%in%c("ID","R3AGEY_B","RAGENDER","SHOCK4",
                                                 "H3ITOT","R3SHLT","R3HLTHLM",
                                                 "R3MSTAT","RAEDEGRM","RADYEAR"))]
struc_dat[[4]] = dat[ , which(colnames(dat)%in%c("ID","R4AGEY_B","RAGENDER","SHOCK5",
                                                 "H4ITOT","R4SHLT","R4HLTHLM",
                                                 "R4MSTAT","RAEDEGRM","RADYEAR"))]
struc_dat[[5]] = dat[ , which(colnames(dat)%in%c("ID","R5AGEY_B","RAGENDER","SHOCK6",
                                                 "H5ITOT","R5SHLT","R5HLTHLM",
                                                 "R5MSTAT","RAEDEGRM","RADYEAR"))]
struc_dat[[6]] = dat[ , which(colnames(dat)%in%c("ID","R6AGEY_B","RAGENDER","H6INPOV","SHOCK7",
                                                 "H6ITOT","R6SHLT","R6HLTHLM",
                                                 "R6MSTAT","RAEDEGRM","RADYEAR"))]
struc_dat[[7]] = dat[ , which(colnames(dat)%in%c("ID","R7AGEY_B","RAGENDER","H7INPOV","SHOCK8",
                                                 "H7ITOT","R7SHLT","R7HLTHLM",
                                                 "R7MSTAT","RAEDEGRM","RADYEAR"))]
struc_dat[[8]] = dat[ , which(colnames(dat)%in%c("ID","R8AGEY_B","RAGENDER","H8INPOV","SHOCK9",
                                                 "H8ITOT","R8SHLT","R8HLTHLM",
                                                 "R8MSTAT","RAEDEGRM","RADYEAR"))]
struc_dat[[9]] = dat[ , which(colnames(dat)%in%c("ID","R9AGEY_B","RAGENDER","H9INPOV","SHOCK10",
                                                 "H9ITOT","R9SHLT","R9HLTHLM",
                                                 "R9MSTAT","RAEDEGRM","RADYEAR"))]
struc_dat[[10]] = dat[ , which(colnames(dat)%in%c("ID","R10AGEY_B","RAGENDER","H10INPOV","SHOCK11",
                                                 "H10ITOT","R10SHLT","R10HLTHLM",
                                                 "R10MSTAT","RAEDEGRM","RADYEAR"))]
struc_dat[[11]] = dat[ , which(colnames(dat)%in%c("ID","R11AGEY_B","RAGENDER","H11INPOV","SHOCK12",
                                                 "H11ITOT","R11SHLT","R11HLTHLM",
                                                 "R11MSTAT","RAEDEGRM","RADYEAR"))]
struc_dat[[12]] = dat[ , which(colnames(dat)%in%c("ID","R12AGEY_B","RAGENDER","H12INPOV","SHOCK13",
                                                 "H12ITOT","R12SHLT","R12HLTHLM",
                                                 "R12MSTAT","RAEDEGRM","RADYEAR"))]

### Test on 1st year
X = struc_dat[[12]]
tmp = apply(X, 1, function(x){
  return((is.na(x[9])))
})
X = X[-which(tmp>0),]
tmp = apply(dat[ ,18:30], 1, function(x){
  return(sum(is.na(x)))
})
which(tmp == 10)


