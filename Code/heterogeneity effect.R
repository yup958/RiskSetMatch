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
  return(as.vector(c(d1,d2,s)))
})
trt$Deathflag = Deathflag[1,]
con$Deathflag = Deathflag[2,]
trt$pairyear = Deathflag[3,]
con$pairyear = Deathflag[3,]

agetrt = sapply(1:dim(trt)[1], function(i){
  return (trt[i, 17 + trt$pairyear[i]])
})
trt$age = agetrt

agecon = sapply(1:dim(con)[1], function(i){
  return (con[i, 17 + con$pairyear[i]])
})
con$age = agecon

healthtrt = sapply(1:dim(trt)[1], function(i){
  return (trt[i, 30 + trt$pairyear[i]])
})
trt$health = healthtrt

healthcon = sapply(1:dim(con)[1], function(i){
  return (con[i, 30 + con$pairyear[i]])
})
con$health = healthcon


#Estimate the ATE by difference-in-means
tauhat_dr = mean(trt$Deathflag) - mean(con$Deathflag)
tauhat_dr


#Estimate the ATE by difference-in-means
tauhat_dr = mean(trt$Deathflag) - mean(con$Deathflag)
tauhat_dr

tauhat_dr_female = mean(trt$Deathflag[which(trt$RAGENDER == '2.Female')]) - 
  mean(con$Deathflag[which(con$RAGENDER == '2.Female')])
tauhat_dr_female 

tauhat_dr_male = mean(trt$Deathflag[which(trt$RAGENDER == '1.Male')]) - 
  mean(con$Deathflag[which(con$RAGENDER == '1.Male')])
tauhat_dr_male 

min(c(trt$age,con$age))#50
max(c(trt$age,con$age))#83

tauhat_dr_age = mean(trt$Deathflag[which(trt$age == 50)]) - 
  mean(con$Deathflag[which(con$age == 50)])
for (age in 51:83){
  tauhat_dr_age = c(tauhat_dr_age, mean(trt$Deathflag[which(trt$age == age)]) - 
                      mean(con$Deathflag[which(con$age == age)]))
}
tauhat_dr_age
plot(50:83, tauhat_dr_age,
     xlab = "Age", ylab = "Treatment Effect", 
     main="ATE vs Age")
lines(50:83,rep(tauhat_dr,34),type = "l", col="blue")


healthcondition = c('1.Excellent', '2.Very good', '3.Good', '4.Fair', '5.Poor')
tauhat_dr_health = numeric(5)
for (i in 1:5){
  tauhat_dr_health[i] = mean(trt$Deathflag[which(trt$health == healthcondition[i])]) - 
    mean(con$Deathflag[which(con$health == healthcondition[i])])
}
plot(1:5, tauhat_dr_health)


