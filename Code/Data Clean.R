setwd("~/Documents/UW/Slides/992/Project/")
library(foreign)

#dat = read.spss("randhrs1992_2016v1_SPSS/randhrs1992_2016v1.sav",to.data.frame = T)
load("hrs.RData")
dat = dat[which(dat$RACOHBYR == "3.Hrs"),]
tmp = as.matrix(dat[, 743:753])

### Reserve HRS sample cohort that have responded to H2-H12
s = apply(tmp,1,function(x){
  t = sapply(x, function(y){
    if(y == "1.Resp,alive"){
      return(1)
    }else if(y == "0.NonResp"){
      return(0)
    }
  })
  return(sum(t))
})
dat = dat[which(s >0),]

### Extract information of respond death year
tmp_label1 = c("RADYEAR", "RABYEAR")


### Extract information of change in wealth and total income
tmp_label2 = c("H2ATOTAC","H3ATOTAC","H4ATOTAC",
              "H5ATOTAC","H6ATOTAC","H7ATOTAC",
              "H8ATOTAC","H9ATOTAC","H10ATOTAC",
              "H11ATOTAC","H12ATOTAC","H13ATOTAC",
              "H1ATOTA","H2ATOTA","H3ATOTA",
              "H4ATOTA","H5ATOTA","H6ATOTA","H7ATOTA","H8ATOTA","H9ATOTA","H10ATOTA",
              "H11ATOTA","H12ATOTA","H13ATOTA")

### Extract information of covariates
tmp_label3 = c("R1AGEY_B","R2AGEY_B","R3AGEY_B","R4AGEY_B","R5AGEY_B","R6AGEY_B",
               "R7AGEY_B","R8AGEY_B","R9AGEY_B","R10AGEY_B","R11AGEY_B","R12AGEY_B",
               "R13AGEY_B","RAGENDER","H1ITOT","H2ITOT","H3ITOT","H4ITOT","H5ITOT","H6ITOT","H7ITOT","H8ITOT","H9ITOT",
               "H10ITOT","H11ITOT","H12ITOT","H13ITOT","H6INPOV","H7INPOV","H8INPOV",
               "H9INPOV","H10INPOV","H11INPOV","H12INPOV","H13INPOV","R1SHLT","R2SHLT",
               "R3SHLT","R4SHLT","R5SHLT","R6SHLT","R7SHLT","R8SHLT","R9SHLT","R10SHLT",
               "R11SHLT","R12SHLT","R13SHLT","R1HLTHLM","R2HLTHLM","R3HLTHLM","R4HLTHLM",
               "R5HLTHLM","R6HLTHLM","R7HLTHLM","R8HLTHLM","R9HLTHLM","R10HLTHLM",
               "R11HLTHLM","R12HLTHLM","R13HLTHLM","R1MSTAT","R2MSTAT","R3MSTAT",
               "R4MSTAT","R5MSTAT","R6MSTAT","R7MSTAT","R8MSTAT","R9MSTAT","R10MSTAT",
               "R11MSTAT","R12MSTAT","R13MSTAT","RAEDEGRM")
label = c(tmp_label1,tmp_label2,tmp_label3)
dat = dat[,which(colnames(dat) %in% label)]


### Wealth Shock Indicator
dat$DeathAge = sapply(1:dim(dat)[1],function(x){
  if(is.na(dat$RADYEAR[x])){
    return(NA)
  } else {
    return(dat$RADYEAR[x] - dat$RABYEAR[x])
  }
})

dat$SHOCK2 = sapply(1:dim(dat)[1],function(x){
  if(is.na(dat$H2ATOTAC[x]) | is.na(dat$H1ATOTA[x])){
    return(NA)
  } else if((dat$H2ATOTAC[x] < 0) && (abs(dat$H2ATOTAC[x]) > 0.75*dat$H1ATOTA[x])){
    return(1)
  } else{
    return(0)
  }
})

dat$SHOCK3 = sapply(1:dim(dat)[1],function(x){
  if(is.na(dat$H3ATOTAC[x]) | is.na(dat$H2ATOTA[x])){
    return(NA)
  } else if((dat$H3ATOTAC[x] < 0) && (abs(dat$H3ATOTAC[x]) > 0.75*dat$H2ATOTA[x])){
    return(1)
  } else{
    return(0)
  }
})

dat$SHOCK4 = sapply(1:dim(dat)[1],function(x){
  if(is.na(dat$H4ATOTAC[x]) | is.na(dat$H3ATOTA[x])){
    return(NA)
  } else if((dat$H4ATOTAC[x] < 0) && (abs(dat$H4ATOTAC[x]) > 0.75*dat$H3ATOTA[x])){
    return(1)
  } else{
    return(0)
  }
})

dat$SHOCK5 = sapply(1:dim(dat)[1],function(x){
  if(is.na(dat$H5ATOTAC[x]) | is.na(dat$H4ATOTA[x])){
    return(NA)
  } else if((dat$H5ATOTAC[x] < 0) && (abs(dat$H5ATOTAC[x]) > 0.75*dat$H4ATOTA[x])){
    return(1)
  } else{
    return(0)
  }
})

dat$SHOCK6 = sapply(1:dim(dat)[1],function(x){
  if(is.na(dat$H6ATOTAC[x]) | is.na(dat$H5ATOTA[x])){
    return(NA)
  } else if((dat$H6ATOTAC[x] < 0) && (abs(dat$H6ATOTAC[x]) > 0.75*dat$H5ATOTA[x])){
    return(1)
  } else{
    return(0)
  }
})

dat$SHOCK7 = sapply(1:dim(dat)[1],function(x){
  if(is.na(dat$H7ATOTAC[x]) | is.na(dat$H6ATOTA[x])){
    return(NA)
  } else if((dat$H7ATOTAC[x] < 0) && (abs(dat$H7ATOTAC[x]) > 0.75*dat$H6ATOTA[x])){
    return(1)
  } else{
    return(0)
  }
})

dat$SHOCK8 = sapply(1:dim(dat)[1],function(x){
  if(is.na(dat$H8ATOTAC[x]) | is.na(dat$H7ATOTA[x])){
    return(NA)
  } else if((dat$H8ATOTAC[x] < 0) && (abs(dat$H8ATOTAC[x]) > 0.75*dat$H7ATOTA[x])){
    return(1)
  } else{
    return(0)
  }
})

dat$SHOCK9 = sapply(1:dim(dat)[1],function(x){
  if(is.na(dat$H9ATOTAC[x]) | is.na(dat$H8ATOTA[x])){
    return(NA)
  } else if((dat$H9ATOTAC[x] < 0) && (abs(dat$H9ATOTAC[x]) > 0.75*dat$H8ATOTA[x])){
    return(1)
  } else{
    return(0)
  }
})

dat$SHOCK10 = sapply(1:dim(dat)[1],function(x){
  if(is.na(dat$H10ATOTAC[x]) | is.na(dat$H9ATOTA[x])){
    return(NA)
  } else if((dat$H10ATOTAC[x] < 0) && (abs(dat$H10ATOTAC[x]) > 0.75*dat$H9ATOTA[x])){
    return(1)
  } else{
    return(0)
  }
})

dat$SHOCK11 = sapply(1:dim(dat)[1],function(x){
  if(is.na(dat$H11ATOTAC[x]) | is.na(dat$H10ATOTA[x])){
    return(NA)
  } else if((dat$H11ATOTAC[x] < 0) && (abs(dat$H11ATOTAC[x]) > 0.75*dat$H10ATOTA[x])){
    return(1)
  } else{
    return(0)
  }
})

dat$SHOCK12 = sapply(1:dim(dat)[1],function(x){
  if(is.na(dat$H12ATOTAC[x]) | is.na(dat$H11ATOTA[x])){
    return(NA)
  } else if((dat$H12ATOTAC[x] < 0) && (abs(dat$H12ATOTAC[x]) > 0.75*dat$H11ATOTA[x])){
    return(1)
  } else{
    return(0)
  }
})

dat$SHOCK13 = sapply(1:dim(dat)[1],function(x){
  if(is.na(dat$H13ATOTAC[x]) | is.na(dat$H12ATOTA[x])){
    return(NA)
  } else if((dat$H13ATOTAC[x] < 0) && (abs(dat$H13ATOTAC[x]) > 0.75*dat$H12ATOTA[x])){
    return(1)
  } else{
    return(0)
  }
})
### Missing Value
for (i in 1:dim(dat)[1]) {
  state = 0
  for (j in 1:12) {
    if(state == 1){
      dat[i,102+j] = 1
    } else if(is.na(dat[i,102+j])){
      dat[i,102+j] = NA
    } else if(dat[i, 102+j] == 1){
      state = 1
    } else{next}
  }
}
dat = dat[,-c(57:81)]
dat = dat[-which(is.na(dat$RAGENDER)),]