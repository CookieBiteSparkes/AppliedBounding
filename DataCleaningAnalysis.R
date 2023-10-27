library(igraph)
library(dplyr)
library(geepack)
library(ergm)
library(gee)
library(intergraph)
library(boot)


##################################################
# Read in data and clean - Forgive the ugly code!
##################################################

raw = read.csv("C:/Users/shane/OneDrive/Desktop/03-sns.csv")
raw = filter(raw, school != 112)
raw$keep1 = 0

idMap = function(x, d){
  
  d2 = filter(d, photoid == x)
  return(d2$id)
  
}

cleanD = function(x){
  
  b = data.frame(x)
  n = nrow(b)
  m = ncol(b)
  for(i in 1:n){
    
    for(j in 9:(m-1)){
      if(!is.na(b[i, j]) && !(b[i, j] %in% b$photoid)){
        b[i, j] = NA
      }
      if(!is.na(b[i, j]) && (b[i, j] %in% b$photoid)){
        b[i, j] = idMap(b[i, j], b)
      }
    }
    
    
  }
  
  return(b)
}

##Schools

sch1 = filter(raw, school == 111)
sch1$id = 1:nrow(sch1)
sch2 = filter(raw, school == 113)
sch2$id = 1:nrow(sch2)
sch3 = filter(raw, school == 114)
sch3$id = 1:nrow(sch3)
sch4 = filter(raw, school == 115)
sch4$id = 1:nrow(sch4)

#### We need to subset to students who non-missing friendship identifications

##School 1 Waves
sch11 = select(sch1, c(1:3, 4, 8, 12, 16, 20, 24:42, 100:101))
sch11 = cleanD(sch11)
sch12 = select(sch1, c(1:3, 5, 9, 13, 17, 21, 43:61, 100:101))
sch12 = cleanD(sch12)
sch13 = select(sch1, c(1:3, 6, 10, 14, 18, 22, 62:80, 100:101))
sch13 = cleanD(sch13)
sch14 = select(sch1, c(1:3, 7, 11, 15, 19, 23, 81:101))
sch14 = cleanD(sch14)

##School 2 Waves
sch21 = select(sch2, c(1:3, 4, 8, 12, 16, 20, 24:42, 100:101))
sch21 = cleanD(sch21)
sch22 = select(sch2, c(1:3, 5, 9, 13, 17, 21, 43:61, 100:101))
sch22 = cleanD(sch22)
sch23 = select(sch2, c(1:3, 6, 10, 14, 18, 22, 62:80, 100:101))
sch23 = cleanD(sch23)
sch24 = select(sch2, c(1:3, 7, 11, 15, 19, 23, 81:101))
sch24 = cleanD(sch24)

##School 3 Waves
sch31 = select(sch3, c(1:3, 4, 8, 12, 16, 20, 24:42, 100:101))
sch31 = cleanD(sch31)
sch32 = select(sch3, c(1:3, 5, 9, 13, 17, 21, 43:61, 100:101))
sch32 = cleanD(sch32)
sch33 = select(sch3, c(1:3, 6, 10, 14, 18, 22, 62:80, 100:101))
sch33 = cleanD(sch33)
sch34 = select(sch3, c(1:3, 7, 11, 15, 19, 23, 81:101))
sch34 = cleanD(sch34)

##School 4 Waves
sch41 = select(sch4, c(1:3, 4, 8, 12, 16, 20, 24:42, 100:101))
sch41 = cleanD(sch41)
sch42 = select(sch4, c(1:3, 5, 9, 13, 17, 21, 43:61, 100:101))
sch42 = cleanD(sch42)
sch43 = select(sch4, c(1:3, 6, 10, 14, 18, 22, 62:80, 100:101))
sch43 = cleanD(sch43)
sch44 = select(sch4, c(1:3, 7, 11, 15, 19, 23, 81:101))
sch44 = cleanD(sch44)


##Network datasets

getNetData = function(x){
  
  b = data.frame(x)
  b$sch_friend120 = b$id
  data11 = reshape(b, direction="long",  varying = c(9:27, 30), sep="")
  dyad2  = data11[!is.na(data11$sch_friend),]  
  dyad1  = as.data.frame(dyad2)
  tdf = distinct(select(dyad1, photoid, contains("eversmk")), photoid, .keep_all = T)
  tdf = tdf[order(tdf$photoid),]
  colnames(tdf)[2] = "eversmk"
  vars   = c("id", "sch_friend")
  dyad   = as.matrix(dyad1[vars])
  net_data  = graph_from_edgelist(dyad)
  net_data  = simplify(net_data)
  V(net_data)$smk = tdf$eversmk
  
  
  
  return(net_data)
  
}


getMeasures = function(x, smk, id){
  
  adj_mat = as_adjacency_matrix(x, sparse = FALSE)
  degO = degree(x, mode = "out")
  degI = degree(x, mode = "in")
  deg = degree(x)
  deg_c = stdz(deg)
  betC = betweenness(x)
  outcls  = closeness(x, mode = "out", normalized = TRUE)    # out closeness
  incls  = closeness(x, mode = "in", normalized = TRUE)     # in closeness
  constr  = constraint(x)            
  
  sym_mat_s = sna::symmetrize(adj_mat, rule="strong")
  net_data_sym_s = graph_from_adjacency_matrix(sym_mat_s, mode="undirected")
  recip   = degree(net_data_sym_s, mode = "out")

  getReachO = function(v){
    length(subcomponent(x, v, mode = "out"))
  }
  getReachI = function(v){
    length(subcomponent(x, v, mode = "in"))
  }
  getReachA = function(v){
    length(subcomponent(x, v, mode = "all"))
  }
  betBrok = function(x, mode){
    bet = betweenness(x, directed = F)
    df = data.frame(bet)
    deg = c()
    K = c()
    ests = c()
    if(mode ==1){
      deg = degree(x)
      K =  unlist(lapply(id, getReachA))
      df = mutate(df, br = if_else(bet == 0, 0, (2*bet + length(deg)-1)/deg))
      ests = df$br
    } else {
      deg = degree(x, mode = "out")
      K =  unlist(lapply(id, getReachO))
      df = mutate(df, br = if_else(bet == 0, 0, (bet + K)/deg))
      ests = df$br
    }
    
  }
  
  broke = betBrok(x, 1)
  fsmk = matrix(smk, nrow = length(smk), 1)
  f_smk = adj_mat%*%fsmk
  measures = data.frame(degI, degO, deg, deg_c, betC, recip, broke, f_smk)
  
  return(measures)
  
}
na.replace = function(x, k){
  
  ifelse(is.na(x), k, x)
}

n11 = getNetData(sch11)
n11M = getMeasures(n11, na.replace(V(n11)$smk, 0), 1:nrow(sch1))
colnames(n11M) = c("indeg1", "outdeg1", "deg1", "degc1", "bet1", "recip1", "broke1", "fsmk1")
n12 = getNetData(sch12)
n12M = getMeasures(n12, na.replace(V(n12)$smk,0), 1:nrow(sch1))
colnames(n12M) = c("indeg2", "outdeg2", "deg2", "degc2", "bet2", "recip2", "broke2", "fsmk2")
n13 = getNetData(sch13)
n13M = getMeasures(n13,  na.replace(V(n13)$smk, 0), 1:nrow(sch1))
colnames(n13M) = c("indeg3", "outdeg3", "deg3", "degc3", "bet3", "recip3", "broke3", "fsmk3")
n14 = getNetData(sch14)
n14M = getMeasures(n14,  na.replace(V(n14)$smk,0) , 1:nrow(sch1))
colnames(n14M) = c("indeg4", "outdeg4", "deg4", "degc4", "bet4", "recip4", "broke4", "fsmk4")


sch1m = cbind(sch1[,1:23], n11M, n12M, n13M, n14M)
sch1m = mutate(sch1m, keep = if_else(is.na(female1) & is.na(female2) & is.na(female3) &
                                       is.na(female4), 0, 1))
sch1m = filter(sch1m, keep ==1)
sch1m$female1 = na.replace(sch1m$female1, 0)
sch1m$female2 = na.replace(sch1m$female2, 0)
sch1m$female3 = na.replace(sch1m$female3, 0)
sch1m$female4 = na.replace(sch1m$female4, 0)
sch1m2 = mutate(sch1m, female = female1 + female2 + female3 + female4)
sch1m2 = mutate(sch1m2, female = if_else(female >= 1, 1, 0))
sch1m2 = select(sch1m2, -c(female1, female2, female3, female4, keep))

sch1f = select(sch1m2, photoid, female, school, hispanic, grades1, grades2, grades3, grades4,
               eversmk1, eversmk2, eversmk3, eversmk4,
               everdrk1, everdrk2, everdrk3, everdrk4,
               home1, home2, home3, home4,
               indeg1, indeg2, indeg3, indeg4,
               outdeg1, outdeg2, outdeg3, outdeg4,
               deg1, deg2, deg3, deg4,
               degc1, degc2, degc3, degc4,
               recip1, recip2, recip3, recip4,
               bet1, bet2, bet3, bet4,
               broke1, broke2, broke3, broke4,
               fsmk1, fsmk2, fsmk3, fsmk4)

sch1f = reshape(sch1f, varying = 5:ncol(sch1f), sep = "", direction = "long")
sch1f = sch1f[order(sch1f$id),]
sch1f  = sch1f[!is.na(sch1f$eversmk),]
sch1f  = sch1f[!is.na(sch1f$everdrk),]
sch1f  = sch1f[!is.na(sch1f$home),]
sch1f  = sch1f[!is.na(sch1f$grades),]
sch1f  = sch1f[!is.na(sch1f$hispanic),]

sch1ff = data.frame(sch1f)


##School 2


n21 = getNetData(sch21)
n21M = getMeasures(n21, na.replace(V(n21)$smk,0),  1:nrow(sch2))
colnames(n21M) = c("indeg1", "outdeg1", "deg1", "degc1", "bet1", "recip1", "broke1", "fsmk1")
n22 = getNetData(sch22)
n22M = getMeasures(n22, na.replace(V(n22)$smk, 0), 1:nrow(sch2))
colnames(n22M) = c("indeg2", "outdeg2", "deg2", "degc2", "bet2", "recip2", "broke2", "fsmk2")
n23 = getNetData(sch23)
n23M = getMeasures(n23, na.replace(V(n23)$smk, 0), 1:nrow(sch2))
colnames(n23M) = c("indeg3", "outdeg3", "deg3", "degc3", "bet3", "recip3", "broke3", "fsmk3")
n24 = getNetData(sch24)
n24M = getMeasures(n24, na.replace(V(n24)$smk, 0), 1:nrow(sch2))
colnames(n24M) = c("indeg4", "outdeg4", "deg4", "degc4", "bet4", "recip4", "broke4", "fsmk4")


sch2m = cbind(sch2[,1:23], n21M, n22M, n23M, n24M)
sch2m = mutate(sch2m, keep = if_else(is.na(female1) & is.na(female2) & is.na(female3) &
                                       is.na(female4), 0, 1))
sch2m = filter(sch2m, keep ==1)
sch2m$female1 = na.replace(sch2m$female1, 0)
sch2m$female2 = na.replace(sch2m$female2, 0)
sch2m$female3 = na.replace(sch2m$female3, 0)
sch2m$female4 = na.replace(sch2m$female4, 0)
sch2m2 = mutate(sch2m, female = female1 + female2 + female3 + female4)
sch2m2 = mutate(sch2m2, female = if_else(female >= 1, 1, 0))
sch2m2 = select(sch2m2, -c(female1, female2, female3, female4, keep))

sch2f = select(sch2m2, photoid, female, school, hispanic, grades1, grades2, grades3, grades4,
               eversmk1, eversmk2, eversmk3, eversmk4,
               everdrk1, everdrk2, everdrk3, everdrk4,
               home1, home2, home3, home4,
               indeg1, indeg2, indeg3, indeg4,
               outdeg1, outdeg2, outdeg3, outdeg4,
               deg1, deg2, deg3, deg4,
               degc1, degc2, degc3, degc4,
               recip1, recip2, recip3, recip4,
               bet1, bet2, bet3, bet4,
               broke1, broke2, broke3, broke4,
               fsmk1, fsmk2, fsmk3, fsmk4)

sch2f = reshape(sch2f, varying = 5:ncol(sch2f), sep = "", direction = "long")
sch2f = sch2f[order(sch2f$id),]
sch2f  = sch2f[!is.na(sch2f$eversmk),]
sch2f  = sch2f[!is.na(sch2f$everdrk),]
sch2f  = sch2f[!is.na(sch2f$home),]
sch2f  = sch2f[!is.na(sch2f$grades),]
sch2f  = sch2f[!is.na(sch2f$hispanic),]

sch2ff = data.frame(sch2f)


## School 3


n31 = getNetData(sch31)
n31M = getMeasures(n31, na.replace(V(n31)$smk, 0), 1:nrow(sch3))
colnames(n31M) = c("indeg1", "outdeg1", "deg1", "degc1", "bet1", "recip1", "broke1", "fsmk1")
n32 = getNetData(sch32)
n32M = getMeasures(n32, na.replace(V(n32)$smk,0), 1:nrow(sch3))
colnames(n32M) = c("indeg2", "outdeg2", "deg2", "degc2", "bet2", "recip2", "broke2", "fsmk2")
n33 = getNetData(sch33)
n33M = getMeasures(n33, na.replace(V(n33)$smk,0), 1:nrow(sch3))
colnames(n33M) = c("indeg3", "outdeg3", "deg3", "degc3", "bet3", "recip3", "broke3", "fsmk3")
n34 = getNetData(sch34)
n34M = getMeasures(n34, na.replace(V(n34)$smk, 0), 1:nrow(sch3))
colnames(n34M) = c("indeg4", "outdeg4", "deg4", "degc4", "bet4", "recip4", "broke4", "fsmk4")


sch3m = cbind(sch3[,1:23], n31M, n32M, n33M, n34M)
sch3m = mutate(sch3m, keep = if_else(is.na(female1) & is.na(female2) & is.na(female3) &
                                       is.na(female4), 0, 1))
sch3m = filter(sch3m, keep ==1)
sch3m$female1 = na.replace(sch3m$female1, 0)
sch3m$female2 = na.replace(sch3m$female2, 0)
sch3m$female3 = na.replace(sch3m$female3, 0)
sch3m$female4 = na.replace(sch3m$female4, 0)
sch3m2 = mutate(sch3m, female = female1 + female2 + female3 + female4)
sch3m2 = mutate(sch3m2, female = if_else(female >= 1, 1, 0))
sch3m2 = select(sch3m2, -c(female1, female2, female3, female4, keep))

sch3f = select(sch3m2, photoid, female, school, hispanic, grades1, grades2, grades3, grades4,
               eversmk1, eversmk2, eversmk3, eversmk4,
               everdrk1, everdrk2, everdrk3, everdrk4,
               home1, home2, home3, home4,
               indeg1, indeg2, indeg3, indeg4,
               outdeg1, outdeg2, outdeg3, outdeg4,
               deg1, deg2, deg3, deg4,
               degc1, degc2, degc3, degc4,
               recip1, recip2, recip3, recip4,
               bet1, bet2, bet3, bet4,
               broke1, broke2, broke3, broke4,
               fsmk1, fsmk2, fsmk3, fsmk4)

sch3f = reshape(sch3f, varying = 5:ncol(sch3f), sep = "", direction = "long")
sch3f = sch3f[order(sch3f$id),]
sch3f  = sch3f[!is.na(sch3f$eversmk),]
sch3f  = sch3f[!is.na(sch3f$everdrk),]
sch3f  = sch3f[!is.na(sch3f$home),]
sch3f  = sch3f[!is.na(sch3f$grades),]
sch3f  = sch3f[!is.na(sch3f$hispanic),]

sch3ff = data.frame(sch3f)



## School 4

n41 = getNetData(sch41)
n41M = getMeasures(n41, na.replace(V(n41)$smk, 0), 1:nrow(sch4))
colnames(n41M) = c("indeg1", "outdeg1", "deg1", "degc1", "bet1", "recip1", "broke1", "fsmk1")
n42 = getNetData(sch42)
n42M = getMeasures(n42, na.replace(V(n42)$smk,0), 1:nrow(sch4))
colnames(n42M) = c("indeg2", "outdeg2", "deg2", "degc2", "bet2", "recip2", "broke2", "fsmk2")
n43 = getNetData(sch43)
n43M = getMeasures(n43, na.replace(V(n43)$smk, 0), 1:nrow(sch4))
colnames(n43M) = c("indeg3", "outdeg3", "deg3", "degc3", "bet3", "recip3", "broke3", "fsmk3")
n44 = getNetData(sch44)
n44M = getMeasures(n44, na.replace(V(n44)$smk, 0), 1:nrow(sch4))
colnames(n44M) = c("indeg4", "outdeg4", "deg4", "degc4", "bet4", "recip4", "broke4", "fsmk4")

na.replace = function(x, k){
  
  ifelse(is.na(x), k, x)
}
sch4m = cbind(sch4[,1:23], n41M, n42M, n43M, n44M)
sch4m = mutate(sch4m, keep = if_else(is.na(female1) & is.na(female2) & is.na(female3) &
                                       is.na(female4), 0, 1))
sch4m = filter(sch4m, keep ==1)
sch4m$female1 = na.replace(sch4m$female1, 0)
sch4m$female2 = na.replace(sch4m$female2, 0)
sch4m$female3 = na.replace(sch4m$female3, 0)
sch4m$female4 = na.replace(sch4m$female4, 0)
sch4m2 = mutate(sch4m, female = female1 + female2 + female3 + female4)
sch4m2 = mutate(sch4m2, female = if_else(female >= 1, 1, 0))
sch4m2 = select(sch4m2, -c(female1, female2, female3, female4, keep))

sch4f = select(sch4m2, photoid, female, school, hispanic, grades1, grades2, grades3, grades4,
               eversmk1, eversmk2, eversmk3, eversmk4,
               everdrk1, everdrk2, everdrk3, everdrk4,
               home1, home2, home3, home4,
               indeg1, indeg2, indeg3, indeg4,
               outdeg1, outdeg2, outdeg3, outdeg4,
               deg1, deg2, deg3, deg4,
               degc1, degc2, degc3, degc4,
               recip1, recip2, recip3, recip4,
               bet1, bet2, bet3, bet4,
               broke1, broke2, broke3, broke4,
               fsmk1, fsmk2, fsmk3, fsmk4)

sch4f = reshape(sch4f, varying = 5:ncol(sch4f), sep = "", direction = "long")
sch4f = sch4f[order(sch4f$id),]
sch4f  = sch4f[!is.na(sch4f$eversmk),]
sch4f  = sch4f[!is.na(sch4f$everdrk),]
sch4f  = sch4f[!is.na(sch4f$home),]
sch4f  = sch4f[!is.na(sch4f$grades),]
sch4f  = sch4f[!is.na(sch4f$hispanic),]

sch4ff = data.frame(sch4f)

final_data = rbind(sch1ff, sch2ff, sch3ff, sch4ff)

#############################################################
#### Construct dataset for analysis 
##############################################################

## School 1:4

pre_f = rbind(cbind(sch14[,1:23], n11M, n12M, n13M, n14M),
              cbind(sch24[,1:23], n21M, n22M, n23M, n24M),
              cbind(sch34[,1:23], n31M, n32M, n33M, n34M),
              cbind(sch44[,1:23], n41M, n42M, n43M, n44M))
dfCN = select(pre_f, school, hispanic, female4, grades4, eversmk4, everdrk4, home4, indeg4,
              outdeg4, deg4, bet4, recip4, broke4, fsmk4)
dfCN = mutate(dfCN, grades = (ifelse(grades4 > 3.78, 1, 0)))


net1 = n14
net2= n24
net3= n34
net4= n44

netW1 = asNetwork(net1)
netW2 = asNetwork(net2)
netW3 = asNetwork(net3)
netW4 = asNetwork(net4)

dfCN$id = 1:nrow(dfCN)
dfN = na.omit(dfCN)
removed = setdiff(1:nrow(dfCN), dfN$id)
dfCN = mutate(dfCN, grades = ifelse(grades4 > 3.78, 1, 0))
dfCN = mutate(dfCN, own = ifelse(home4 == 1, 1, 0))



adj1 = as.matrix(as_adjacency_matrix(n14))
diag(adj1) = NA
adj1 = t(adj1)
links1 = as.vector(adj1)

adj2 = as.matrix(as_adjacency_matrix(n24))
diag(adj2) = NA
adj2 = t(adj2)
links2 = as.vector(adj2)

adj3 = as.matrix(as_adjacency_matrix(n34))
diag(adj3) = NA
adj3 = t(adj3)
links3 = as.vector(adj3)

adj4 = as.matrix(as_adjacency_matrix(n44))
diag(adj4) = NA
adj4 = t(adj4)
links4 = as.vector(adj4)

## Make school datasets
schoolData = function(data, schoolID, links, num, start){
  dfC = filter(data, school == schoolID)
  n = nrow(dfC)
  
  
  H_g = mat.or.vec(n,n)
  H_smk = mat.or.vec(n,n)
  H_drk = mat.or.vec(n,n)
  H_hisp = mat.or.vec(n,n)
  H_own = mat.or.vec(n,n)
  H_grades = mat.or.vec(n,n)
  
  
  for(s in 1:n){
    for(t in 1:n){
      if(!is.na(dfC$female4[s]) & !is.na(dfC$female4[t])){
        if(dfC$female4[s] == dfC$female4[t]){
          H_g[s, t] = 1
        }
      } else{
        H_g[s, t] = NA
      }
      
      if(!is.na(dfC$eversmk4[s]) & !is.na(dfC$eversmk4[t])){
        if(dfC$eversmk4[s] == 1 & dfC$eversmk4[t] == 1){
          H_smk[s, t] = 1
        }
      } else{
        H_smk[s, t] = NA
      }
      
      
      if(!is.na(dfC$everdrk4[s]) & !is.na(dfC$everdrk4[t])){
        if(dfC$everdrk4[s] == 1 & dfC$everdrk4[t] == 1){
          H_drk[s, t] = 1
        }
      } else{
        H_drk[s, t] = NA
      }
      
      if(!is.na(dfC$hispanic[s]) & !is.na(dfC$hispanic[t])){
        if(dfC$hispanic[s] == dfC$hispanic[t]){
          H_hisp[s, t] = 1
        }
      } else{
        H_hisp[s, t] = NA
      }
      
      if(!is.na(dfC$grades[s]) & !is.na(dfC$grades[t])){
        if(dfC$grades[s] == dfC$grades[t]){
          H_grades[s, t] = 1
        }
      } else{
        H_grades[s, t] = NA
      }
      
      if(!is.na(dfC$own[s]) & !is.na(dfC$own[t])){
        if(dfC$own[s] == dfC$own[t]){
          H_own[s, t] = 1
        }
      } else{
        H_own[s, t] = NA
      }
      
      
    }
  }
  dfC$nodeID = 1:n
  dft = data.frame(Links = links,
                   DrinkR = as.vector((matrix(dfC$everdrk4, nrow = nrow(dfC), ncol = nrow(dfC)))),
                   HispanicR = as.vector((matrix(dfC$hispanic, nrow = nrow(dfC),ncol = nrow(dfC)))),
                   SmokeR = as.vector((matrix(dfC$eversmk4,nrow = nrow(dfC), ncol = nrow(dfC)))),
                   GenderR = as.vector((matrix(dfC$female4,nrow = nrow(dfC), ncol = nrow(dfC)))),
                   OwnR = as.vector((matrix(dfC$own,nrow = nrow(dfC), ncol = nrow(dfC)))),
                   GradesR = as.vector((matrix(dfC$grades,nrow = nrow(dfC), ncol = nrow(dfC)))),
                   DrinkS = as.vector(t(matrix(dfC$everdrk4, nrow = nrow(dfC), ncol = nrow(dfC)))),
                   HispanicS = as.vector(t(matrix(dfC$hispanic, nrow = nrow(dfC),ncol = nrow(dfC)))),
                   SmokeS = as.vector(t(matrix(dfC$eversmk4,nrow = nrow(dfC), ncol = nrow(dfC)))),
                   GenderS = as.vector(t(matrix(dfC$female4,nrow = nrow(dfC), ncol = nrow(dfC)))),
                   OwnS = as.vector(t(matrix(dfC$own,nrow = nrow(dfC), ncol = nrow(dfC)))),
                   GradesS = as.vector(t(matrix(dfC$grades,nrow = nrow(dfC), ncol = nrow(dfC)))))
  
  dft = mutate(dft, DrinkH = as.vector(t(H_drk)))
  dft = mutate(dft, HispanicH = as.vector(H_hisp))
  dft = mutate(dft, SmokeH = as.vector(H_smk))
  dft = mutate(dft, GenderH = as.vector(H_g))
  dft = mutate(dft, OwnH = as.vector(H_own))
  dft = mutate(dft, GradesH = as.vector(H_grades))
  
  id = (start+1):(start + n)
  dft$id = as.vector((matrix(id, nrow = nrow(dfC), ncol = nrow(dfC))))
  dft$school = rep(num, n*n)
  dft = arrange(dft, id)
  
  dft$H2 = ifelse(dft$HispanicS ==1 & dft$HispanicR ==1, 1, 0)
 
  return(dft) 
}

df1 = schoolData(dfCN, schoolID=111, links1, num=1, start = 0)
df2 = schoolData(dfCN, schoolID=113, links2, num=2, start = nrow(df1))
df3 = schoolData(dfCN, schoolID=114, links3, num=3, start = nrow(df1) + nrow(df2))
df4 = schoolData(dfCN, schoolID=115, links4, num=4, start = nrow(df1) + nrow(df2) + nrow(df3))

df11 =na.omit(df1)
df22 =na.omit(df2)
df33 =na.omit(df3)
df44 =na.omit(df4)
df = rbind(df11, df22, df33, df44)

df$school = factor(df$school)
df$bid = 1:nrow(df)
N=nrow(df)

#########################################################################
# Propensity scores and exch. correlation estimation
###########################################################################

### Get propensity scores
pmodelGH = glm(GenderH*HispanicH ~ SmokeS + SmokeR + SmokeH +
               DrinkR + DrinkS + DrinkH + OwnR + OwnS + OwnH +
               GradesR + GradesS + GradesH + school, family = binomial, data = df)
pscoresGH = predict(pmodelGH, type = "response")
df$pscoresGH = pscoresGH
df$strataGH = factor(ntile(pscoresGH, ceiling(N/(log(N)))))

pmodelS = glm(SmokeH ~ HispanicS + HispanicR + HispanicH + GenderS + GenderR + GenderH +
                DrinkR + DrinkS + DrinkH + OwnR + OwnS + OwnH +
                GradesR + GradesS + GradesH + school, family = binomial, data = df)
pscoresS = predict(pmodelS, type = "response")
df$pscoresS = pscoresS
df$strataS = factor(ntile(pscoresS, ceiling(N/(log(N)))))


df = arrange(df, strataGH)
GH = gee(Links ~ HispanicS + HispanicR + HispanicH + GenderS + GenderR + GenderH +
           DrinkR + DrinkS + DrinkH +
           SmokeS + SmokeR + SmokeH + OwnR + OwnS + OwnH +
           GradesR + GradesS + GradesH + school,family = binomial(link="logit"), corstr = "exchangeable",id = strataGH, df)

df = arrange(df, strataS)
Sm = gee(Links ~ HispanicS + HispanicR + HispanicH + GenderS + GenderR + GenderH +
           DrinkR + DrinkS + DrinkH +
           SmokeS + SmokeR + SmokeH + OwnR + OwnS + OwnH +
           GradesR + GradesS + GradesH + school,family = binomial(link="logit"), corstr = "exchangeable",id = strataS, df)


df = arrange(df, id)
Eg = gee(Links ~ HispanicS + HispanicR + HispanicH + GenderS + GenderR + GenderH +
           DrinkR + DrinkS + DrinkH +
           SmokeS + SmokeR + SmokeH + OwnR + OwnS + OwnH +
           GradesR + GradesS + GradesH + school,family = binomial(link="logit"), corstr = "exchangeable",id = id, df)

ExCoef = c(GH$working.correlation[1,2], Sm$working.correlation[1,2], Eg$working.correlation[1,2])
Ec = max(ExCoef)


### Main model, boots, and CIs

mainModel = glm(Links ~ SmokeS + SmokeR + SmokeH + DrinkR + DrinkS + DrinkH +
                  HispanicS + HispanicR + HispanicH + GenderS + GenderR + GenderH +
                  OwnR + OwnS + OwnH +
                  GradesR + GradesS + GradesH + school, family = binomial, data=df)
mSum = summary(mainModel)
L1 = mSum$coefficients[,1] - 1.96*sqrt(1 + Ec*sqrt(N-1))*mSum$coefficients[,2]
U1 = mSum$coefficients[,1] + 1.96*sqrt(1 + Ec*sqrt(N-1))*mSum$coefficients[,2]

## Hoeffding boot

fitM = function(data, indices){
  
  
  fM = function(data){
    model = glm(Links ~ HispanicS + HispanicR + HispanicH + GenderS + GenderR + GenderH +
                  DrinkR + DrinkS + DrinkH +
                  SmokeS + SmokeR + SmokeH + OwnR + OwnS + OwnH +
                  GradesR + GradesS + GradesH + school, family = binomial(link="logit"), data=data)
    return(coef(model))
  }
  dfo = data[indices,]
  
  return(fM(dfo))
}


boots = boot(df, fitM, R =5000, parallel = "snow")
Rs = apply(rbind(boots$t, matrix(boots$t0, nrow = 1)), 2, max) - apply(rbind(boots$t, matrix(boots$t0, nrow = 1)), 2, min)
L2 = boots$t0 - Rs*sqrt(log(2/.05)/6)
U2 = boots$t0 + Rs*sqrt(log(2/.05)/6)