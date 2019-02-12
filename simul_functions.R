# Appendix 1. R code and functions to apply the cookbook to analyse biological responses 
#             to multiple stressors 
#
# Analysing the impact of multiple stressors on aquatic biology: cookbook with applications in R 
# C. K. Feld, P. Segurado & C. Gutierrez-Canovas
#
# Functions to simulate multi-stressor data and test the influence of sample size and
# gradient length on stressor hierarchy, goodness-of-fit and SES estimation
#
# Code written by C. K. Feld, P. Segurado & C. Gutierrez-Canovas 
#
# Problems with the code, send email to GutierrezCanovasC@cardiff.ac.uk

# Loading required libraries
# The required packages should be installed

library(MuMIn) # Multi-Model Inference
library(randomForestSRC) # RF
library(ggRandomForests) # RF
library(gbm) # BRT
library(dismo) # BRT

# sim.multi.str() produces artificial multi-stressor data
# This function computes the hypervolume to estimate how each community fills
# the functional space
#
# Inputs:
# n: number of sites
# ses: a vector indicanting the SES of each variable
# ac: a number indicating the SES of the model error
# mod.type: the model type. It could be "antagonistic","synergistic","opposing","additive","complex"
#
# plot.int: logical value indicating wether plotting the model responses along stressor 1 and 2 gradients
# 
# Output:
# a list containing:
# sim.dat: simulated multi-stressor data
# mod: linear model accounting for the defined model
# var.names: response and predictor names

sim.multi.str<-function(n=100, ses, ac, 
                        mod.type=c("antagonistic","synergistic","opposing","additive","complex"), 
                        plot.int=F){
  
  # Creating 20 stressors
  s1<-rnorm(n)
  s2<-rnorm(n)
  s3<-rnorm(n)
  s4<-rnorm(n)
  s5<- s1 + 2*rnorm(n)
  s6<- s1 + 2*rnorm(n)
  s7<- s2 + 2*rnorm(n)
  s8<- s2 + 2*rnorm(n)
  s9<- s3 + 2*rnorm(n)
  s10<- s3 + 2*rnorm(n)
  s11<- s4 + 2*rnorm(n)
  s12<- s4 + 2*rnorm(n)
  s13<- rnorm(n)
  s14<- rnorm(n)
  s15<- rnorm(n)
  s16<- rnorm(n)
  s17<- rnorm(n)
  s18<- rnorm(n)
  s19<- rnorm(n)
  s20<- rnorm(n)
  
  err<-ac*rnorm(n) # global stressor effect
  
  # Defining model and response variable
  if (mod.type=="mixed") y =  (- ses[1]*s1 - ses[2]*s2 - ses[3]*s3 + ses[4]*s4 
                                      - ses[5]*s1*s2 + ses[6]*s1*s3 + ses[7]*s2*s4 
                                      + err)
  
  if (mod.type=="antagonistic") y =  - ses[1]*s1 - ses[2]*s2 + ses[3]*s1*s2 + err
  
  if (mod.type=="synergistic") y =  - ses[1]*s1 - ses[2]*s2 - ses[3]*s1*s2 + err
  
  if (mod.type=="opposing") y =  + ses[1]*s1 - ses[2]*s2 - ses[3]*s1*s2 + err
  
  if (mod.type=="additive") y =  - ses[1]*s1 - ses[2]*s2 + err 
  
  
  # dataset of predictors
  raw.predictors<-data.frame(s1,s2,s3,s4,
                             s5,s6,s7,s8,s9,s10,s11,s12,
                             s13,s14,s15,s16,s17,s18,s19,s20)
  
  sim.dat<-data.frame(y,scale(raw.predictors)) # standardising predictors
  
  # Running linear model
  if (mod.type=="mixed") mod<-lm (y ~ s1 + s2 + s3 + s4 + s1:s2 + s1:s3 + s2:s4, data=sim.dat)
  if (mod.type!="mixed") mod<-lm (y ~ s1 + s2 + s1:s2, data=sim.dat)
  
  if (mod.type=="mixed") var.names<-c("s1","s2","s3","s4","s1:s2","s1:s3","s2:s4")
  if (mod.type!="mixed") var.names<-c("s1","s2","s1:s2")
  
  if(plot.int==T) {
  
  # Setting graphical parameters
  par(mfrow=c(1,2),cex=1,cex.axis=1,cex.lab=1, pch=16)
  col.int<-rainbow(3)
  
  # Generating the whole gradient for the first stressor
  s1.seq<-seq(min(sim.dat$s1),max(sim.dat$s1),length.out = 1000)
  s2.seq<-seq(min(sim.dat$s2),max(sim.dat$s2),length.out = 1000)
  
  #plotting legend in the graph
  
   
  # Scatter plot 1
  plot(sim.dat$s1, sim.dat$y, xlab="stressor 1",ylab="index")
  
  # fitted values
  
  # 90th s2, mean s3
  lines(s1.seq,predict(mod,data.frame(s1=s1.seq,
                                      s2=as.numeric(quantile(sim.dat$s2,0.90)),
                                      s3=mean(sim.dat$s3),
                                      s4=mean(sim.dat$s4))),
                                      lty=1,lwd=3,col=col.int[3])
  # 50th s2, mean s3
  lines(s1.seq,predict(mod,data.frame(s1=s1.seq,
                                      s2=as.numeric(quantile(sim.dat$s2,0.50)),
                                      s3=mean(sim.dat$s3),
                                      s4=mean(sim.dat$s4))),
                                      lty=1,lwd=3,col=col.int[2])
  
  # 10th s2, mean s3
  lines(s1.seq,predict(mod,data.frame(s1=s1.seq,
                                      s2=as.numeric(quantile(sim.dat$s2,0.10)),
                                      s3=mean(sim.dat$s3),
                                      s4=mean(sim.dat$s4))),
                                      lty=1,lwd=3,col=col.int[1])
  
    ## Legend for plot 1
  
  if(mod.type!="opposing") legend( "bottomleft",
                                   legend = c("s2 low","s2 intermediate","s2 high"),
                                   bty = "n", col = col.int,
                                   lwd=3, lty = 1)
  
  if(mod.type=="opposing") legend("topleft",
                                  legend = c("s2 low","s2 intermediate","s2 high"),
                                  bty = "n", col = col.int,
                                  lwd=3, lty = 1)
  
  # Scatter plot 2
  plot(sim.dat$s2, sim.dat$y, xlab="stressor 2",ylab="index")
  
  # fitted values
  
  # 90th s1, mean s3
  lines(s2.seq,predict(mod,data.frame(s1=as.numeric(quantile(sim.dat$s1,0.90)),
                                      s2=s2.seq,
                                      s3=mean(sim.dat$s3),
                                      s4=mean(sim.dat$s4))),
                                      lty=1,lwd=3,col=col.int[3])
  # 50th s1, mean s3
  lines(s2.seq,predict(mod,data.frame(s1=as.numeric(quantile(sim.dat$s1,0.50)),
                                      s2=s2.seq,
                                      s3=mean(sim.dat$s3),
                                      s4=mean(sim.dat$s4))),
                                      lty=1,lwd=3,col=col.int[2])
  
  # 10th s1, mean s3
  lines(s2.seq,predict(mod,data.frame(s1=as.numeric(quantile(sim.dat$s1,0.10)),
                                      s2=s2.seq,
                                      s3=mean(sim.dat$s3),
                                      s4=mean(sim.dat$s4))),
                                      lty=1,lwd=3,col=col.int[1])
  
  legend( "bottomleft", legend = c("s1 low","s1 intermediate","s1 high"),
                                   bty = "n", col = col.int,
                                   lwd=3, lty = 1)
  
  
  
  }
  par(mfrow=c(1,1))
  return(list(sim.dat=sim.dat, mod=mod, var.names=var.names))
}

# qual_analysis() runs BRT and RF analyses
#
# Inputs:
# sim.set: simulated multi-stressor data obtained using sim.multi.str() 
# ntree: number of trees for RF analysis
# lr: learning rate for BRT analysis
# BRT_mod: logical value indicating wether to perform BRT analysis. 
# When sample size is > 100 this is set automatically to FALSE
#
# RF_mod: logical value indicating wether to perform RF analysis. 
# 
# Output:
#
# res_brt: BRT results
# res_rf: RF results
# res_lm_brt: linear model averaging results for the predictors selected using BRT
# res_lm_rf: linear model averaging results for the predictors selected using RF

qual_analysis<-function(sim.set, ntree, lr, BRT_mod, RF_mod){
  
  res.names<-c("s1","s2","s3","s4","int1","int2","int3","num.vars","max.cor","cv","r2","RMSE","ntree")
  res_rf<-res_brt<-data.frame(matrix(NA,1,length(res.names)))
  names(res_rf)<-names(res_brt)<-res.names
  
  lm.names<-c("s1","s2","s3","s4","s1:s2","s1:s3","s2:s4","r.cor","r2")
  res_lm_brt<-res_lm_rf<-data.frame(matrix(NA,1,length(lm.names)))
  names(res_lm_brt)<-names(res_lm_rf)<-lm.names
  
  sim.set$sim.dat->sim.dat
  sim.set$var.names->var.names
  
  max(abs(as.dist(cor(sim.dat[,-1]))))->res_rf[9]->res_brt[9] # storing max pairwise correlation
  
  p<- ncol(sim.dat)-1 # number of predictors
  
  if(BRT_mod==T){
    
    my.gbm.full <- gbm.step(data=sim.dat, gbm.x=c(2:ncol(sim.dat)), gbm.y="y", 
                            family="gaussian", tree.complexity = 2,
                            learning.rate=lr, bag.fraction=0.66)
    

    # Explained deviance based on CV
    my.gbm.full$self.statistics$mean.null->null.dev
    my.gbm.full$cv.statistics$deviance.mean->resid.dev
    
    1-resid.dev/null.dev->res_brt[11]
    
    # RMSE
    sqrt(1/nrow(sim.dat) * sum((sim.dat$y - my.gbm.full$fitted)^2))->res_brt[12]
    
    my.gbm.full$n.trees->res_brt[13] # number of trees
    
    res_brt[8] <- p # number of variables
    
    # Variable importance
    abs(rank(summary(my.gbm.full)[,2])-p-1)->rank.tmp
    summary(my.gbm.full)[,1]->names(rank.tmp)
    rank.tmp[paste("s",1:p,sep="")][1:4]->res_brt[1:4]
    
    sel.st <- names(rank.tmp)[1:4]
    
    # compute pairwise interactions among explanatory variables
    int.tmp <- gbm.interactions(my.gbm.full)$rank.list
    
    if(nrow(int.tmp)<3) int.length<-nrow(int.tmp) else int.length<-3 # In case there are less than 3 interactions
    
    
    for (a in 1:int.length) {
      if(int.tmp[a,1] > int.tmp[a,3]) paste(int.tmp[a,4],int.tmp[a,2],sep=":")->res_brt[a+4] # To ensure right order
      
      if(int.tmp[a,1] < int.tmp[a,3]) paste(int.tmp[a,2],int.tmp[a,4],sep=":")->res_brt[a+4]
      
    }
    
    sel.int <- paste(res_brt[5:7])
    
    lm_est(sel.st, sel.int, sim.dat, var.names)-> res_lm_brt
    
  } else res_brt<-NULL
  
  
  #### RF
  
  if (RF_mod==T) {
    
    my.rf <- rfsrc(y ~ ., mtry = ceiling(p/3), ntree = ntree, importance = "permute",  
                   data = sim.dat)
    
    # CV error
    res_rf[10] <- round(my.rf$err.rate[ntree],1)
    
    # RMSE
    sqrt(1/nrow(sim.dat) * sum((sim.dat$y - my.rf$predicted)^2))->res_rf[12]
    
    ntree->res_rf[13] # number of trees
    
    # variable importance (VIPM)
    abs(rank(my.rf$importance)-p-1)->rank.tmp
    rank.tmp[paste("s",1:p,sep="")][1:4]->res_rf[1:4] # storing variable rank for s1 to s4
    
    # selecting best stressor subset
    sel.st <- max.subtree(my.rf)$topvars # top variables
    res_rf[8] <- length(sel.st) # number of selected variales by max.subtree analysis
    
    sort(names(sort(rank.tmp)[1:4]))-> sel.st # best 4 candidates to be used for interaction detection and model averaging
    
    if(res_rf[8]>1) my.rf.interaction <- find.interaction(my.rf, xvar.names = sel.st, 
                                          importance= "permute", 
                                          method = "vimp", nrep = 3)
    
    if(res_rf[8]<2) int.length<-0 else nrow(my.rf.interaction)->int.length # number of interactions
    
    
    
    # Just in case there are less than 3 interactions
    
    if(int.length==0) rep(NA,3)->sel.int
    
    if(int.length==1) paste(sel.st[1:2],collapse =":")->sel.int->res_rf[5]
    
    if(int.length==2) names(sort(abs(my.rf.interaction[,5]),decreasing=T)[1:2])->sel.int->res_rf[5:6]
    
    if(int.length>2) names(sort(abs(my.rf.interaction[,5]),decreasing=T)[1:3])->sel.int->res_rf[5:7]
    
    # Model averaging
    
    if(res_rf[8] > 0 & int.length == 0) lm_est(sel.st=sel.st, sim.dat=sim.dat, var.names=var.names)-> res_lm_rf
    
    if(res_rf[8] > 0 & int.length > 0) lm_est(sel.st, sel.int, sim.dat, var.names)-> res_lm_rf
  
  } else res_rf<-NULL
  
  
  return(list(res_brt=res_brt, res_rf=res_rf, res_lm_brt=res_lm_brt, res_lm_rf=res_lm_rf))
}

# lm_est() produces the average lineal model from a global model
#
# Inputs:
#
# sel.st: selected single stressors
# sel.int: selected stressor interactions
# sim.dat: simulated data matrix
# var.names:  variable names
# 
# Output:
#
# SES for the variables and interactions of interest

lm_est<-function(sel.st=NULL, sel.int=NULL, sim.dat, var.names) {
  
  res<-rep(NA,length(var.names))
  
  # Creating the model formula with the selected variables (stressors + interactions)
  form.mod<-as.formula(paste("y", paste(c(sel.st,sel.int), collapse=" + "), sep=" ~ "))
  
  # Defining the global model
  mod<-lm (form.mod, data=sim.dat)
  
  options(na.action = "na.fail") # necessary to run dredge()
  
  # Ranking models from variable subsets from the global model
  mod.d<-dredge(mod, rank = "AIC",extra = 
                  c(R2=function(x) (summary(x)$adj.r.squared)))
  
  # Selecting models with AIC delta <=6
  mod.d.gm<-get.models(mod.d, subset=delta<=6)
  
  length(mod.d.gm)->n.mod # number of models
  weighted.mean(mod.d$R2[1:n.mod],mod.d$weight[1:n.mod])->r2.w # weighted r2
  
  # Model averaging when length(mod.d.gm) > 1, otherwise selecting the best model
  
  if (length(mod.d.gm) > 1) {
    
    # Model averaging, extracting coefficients
    model.avg(mod.d.gm,revised.var = TRUE)->mod.av
    coef.av<-mod.av$coefficients[1,-1]
    cor(predict(mod.av),sim.dat$y)->r.cor
    
    # Storing model coefficients from target variables and interactions
    
    for (j in var.names){
      res.tmp<-coef.av[j]
      if (is.na(res.tmp)==F) res[which(var.names==j)]<-res.tmp else res[which(var.names==j)]<-0
    }
  }
  
  if (length(mod.d.gm) == 1) {
    
    summary(mod.d.gm[[1]])$coef[,1]->coef.av # extracting the best model
    cor(predict(mod.d.gm[[1]]),sim.dat$y)->r.cor
    # Storing model coefficients from target variables and interactions
    
    for (j in var.names){
      res.tmp<-coef.av[j]
      if (is.na(res.tmp)==F) res[which(var.names==j)]<-res.tmp else res[which(var.names==j)]<-0
    }
  }
  
  return(c(res,r.cor,r2.w))
  
}

# sim_ntest() runs a number of simulations to assess the effect of varying sample size
# in multi-stressor analysis
#
# Inputs:
# sim.set: simulated multi-stressor data obtained using sim.multi.str()
# n.test: number of sites to subset the full dataset
# ntree: number of trees for RF analysis
# lr: learning rate for BRT analysis
# BRT_mod: logical value indicating wether to perform BRT analysis. 
# When sample size is > 100 this is set automatically to FALSE
# RF_mod: logical value indicating wether to perform RF analysis.
# runs: number of simulations
# 
# Output:
#
# res_brt: BRT results
# res_rf: RF results
# res_lm_brt: linear model averaging results for the predictors selected using BRT
# res_lm_rf: linear model averaging results for the predictors selected using RF
# ntest: number of sites to subset the full dataset
# runs: number of simulations


sim_ntest<-function(sim.set,n.test,ntree, lr, BRT_mod, RF_mod, runs){
  
  res.names<-c("s1","s2","s3","s4","int1","int2","int3","num.vars","max.cor","cv","r2","RMSE","ntree")
  sim_res_brt<-sim_res_rf<-data.frame(matrix(NA,runs,length(res.names)))
  colnames(sim_res_brt)<-colnames(sim_res_rf)<-res.names
  
  lm.names<-c("s1","s2","s3","s4","s1:s2","s1:s3","s2:s4","r.cor","r2")
  sim_res_lm_brt<-sim_res_lm_rf<-data.frame(matrix(NA,1,length(lm.names)))
  names(sim_res_lm_brt)<-names(sim_res_lm_rf)<-lm.names
  
  if (n.test<100) BRT_mod<-F else BRT_mod<-T
  
  for (i in 1:runs) {
    
    if (n.test==F) sim.set$sim.dat->sim.dat.ntest else sim.set$sim.dat[sample(1:nrow(sim.set$sim.dat),n.test),]->sim.dat.ntest 
    
    sim.set.ntest<-list(sim.dat=sim.dat.ntest,var.names=sim.set$var.names)
    
    qual_analysis(sim.set.ntest, ntree, lr, BRT_mod,RF_mod)->simul.res
    
    if (n.test>=100) simul.res$res_brt->sim_res_brt[i,]
    simul.res$res_rf->sim_res_rf[i,]
    if (n.test>=100)simul.res$res_lm_brt->sim_res_lm_brt[i,]
    simul.res$res_lm_rf->sim_res_lm_rf[i,]
  }
  return(list(res_brt=sim_res_brt, res_rf=sim_res_rf, res_lm_brt=sim_res_lm_brt, res_lm_rf=sim_res_lm_rf,
              n.test=n.test,runs=runs))
}

# sim_length() runs a number of simulations to assess the effect of varying sample size
# and gradient length in multi-stressor analysis
#
# Inputs:
# sim.set: simulated multi-stressor data obtained using sim.multi.str() 
# ntree: number of trees for RF analysis
# n.test: number of sites to subset the full dataset
# stressor: stressor used to subset the data
# st_leng: gradient length expressed as quantile. e.g. 0.25 represents 25% of the data range
# st_level: stressor intensity from which the gradient will start. O means
# lr: learning rate for BRT analysis
# BRT_mod: logical value indicating wether to perform BRT analysis. 
# When sample size is > 100 this is set automatically to FALSE
# RF_mod: logical value indicating wether to perform RF analysis.
# runs: number of simulations
# 
# Output:
#
# res_brt: BRT results
# res_rf: RF results
# res_lm_brt: linear model averaging results for the predictors selected using BRT
# res_lm_rf: linear model averaging results for the predictors selected using RF


sim_length<-function(sim.set,n.test,stressor,st_leng,st_level,ntree, lr, BRT_mod, RF_mod, runs){
  
  
  res.names<-c("s1","s2","s3","s4","int1","int2","int3","num.vars","max.cor","cv","r2","RMSE","ntree")
  sim_res_brt<-sim_res_rf<-data.frame(matrix(NA,runs,length(res.names)))
  colnames(sim_res_brt)<-colnames(sim_res_rf)<-res.names
  
  lm.names<-c("s1","s2","s3","s4","s1:s2","s1:s3","s2:s4","r.cor","r2")
  sim_res_lm_brt<-sim_res_lm_rf<-data.frame(matrix(NA,1,length(lm.names)))
  names(sim_res_lm_brt)<-names(sim_res_lm_rf)<-lm.names
  
  n.sites<-rep(NA,runs)
  sim.set$sim.dat->sim.dat
  
    for (i in 1:runs) {
    
    unlist(sim.dat[stressor])->st
    l_set<-as.numeric(which((st>=quantile(st,st_level)) & (st<quantile(st,st_level+st_leng))))
    
    sim.dat[l_set,]->sim.dat.l
    
    sim.dat.l[sample(1:nrow(sim.dat.l),n.test),]->sim.dat.ntest
    
    nrow(sim.dat.ntest)->n.sites[i]
    
    if (n.test>n.sites[i]) {
      warning("The number of sites is smaller than the n.test parameter. n.sites is used instead.")
      n.test<-n.sites[i]
    }
    
    
    if (n.sites[i]<100) BRT_mod<-F else BRT_mod<-T
    
    sim.set.l<-list(sim.dat=sim.dat.ntest,var.names=sim.set$var.names)
    
    qual_analysis(sim.set.l, ntree, lr, BRT_mod, RF_mod)->simul.res
    
    if (BRT_mod==T) simul.res$res_brt->sim_res_brt[i,]
    simul.res$res_rf->sim_res_rf[i,]
    if (BRT_mod==T) simul.res$res_lm_brt->sim_res_lm_brt[i,]
    simul.res$res_lm_rf->sim_res_lm_rf[i,]
  }
  return(list(res_brt=sim_res_brt, res_rf=sim_res_rf, res_lm_brt=sim_res_lm_brt, res_lm_rf=sim_res_lm_rf, n.sites=n.sites, 
              n.test=n.test,stressor=stressor,st_leng=st_leng,st_level=st_level,runs=runs))
}

# sim_length() runs a number of simulations to assess the effect of varying sample size
# and gradient length in multi-stressor analysis
#
# Inputs:
#
# res.sets: list of simulation results
# len_sim: vector indicating wich simulations were performed using sim_length()
# ses: Standardised Effect Sizes of the putative stressors, which were defined a priori
# BRT_mods: vector indicating wich simulations were performed using BRT
# RF_mods: vector indicating wich simulations were performed using RF
# 
# Output:
#
# brt_qual_sim_res: overall results for BRT simulations
# rf_qual_sim_res: overall results for BRT simulations
# lm_brt_qual_sim_res: overall results for LM simulations using BRT-based candidates
# lm_rf_qual_sim_res: overall results for BRT simulations using RF-based candidates
# brt_mean_rank: mean rank for the single putative stressors (s1-s4) in BRT analysis
# rf_mean_rank: mean rank for the single putative stressors (s1-s4) in RF analysis
# top_int_brt: frequency of occurrence of the interacions s1:S2,s1:s3 and 
#              s2:S4 within the top three interactions in BRT analysis
# top_int_rf: frequency of occurrence of the interacions s1:S2,s1:s3 and 
#             s2:S4 within the top three interactions in RF analysis
# range.trees: Minimum and maximum number of trees used in BRT analysis for 
#              each combinationof sample size and gradient length

sim_res<-function(res.sets, len_sim, ses, BRT_mods, RF_mods){
  
  res.names<-c("n.test","n.sites", "st_leng", "st_level","s1","s2","s3","s4","s1:s2","s1:s3","s2:s4","num.vars","max.cor","cv","r2","RMSE","ntree")
  brt_qual_sim_res<-rf_qual_sim_res<-data.frame(matrix(NA,1,length(res.names)))
  colnames(brt_qual_sim_res)<-colnames(rf_qual_sim_res)<-res.names
  
  mean_rank_names<-c("s1","s2","s3","s4")
  brt_mean_rank<-rf_mean_rank<-data.frame(matrix(NA,1,length(mean_rank_names)))
  colnames(brt_mean_rank)<-colnames(rf_mean_rank)<-mean_rank_names
  
  lm.names<-c("n.test","n.sites", "st_leng", "st_level","s1","s2","s3","s4","s1:s2","s1:s3","s2:s4","r.cor","r2")
  lm_brt_qual_sim_res<-lm_rf_qual_sim_res<-data.frame(matrix(NA,1,length(lm.names)))
  names(lm_brt_qual_sim_res)<-names(lm_rf_qual_sim_res)<-lm.names
  
  top_int_names<-c("s1:s2","s1:s3","s2:s4")
  top_int_brt<-top_int_rf<-data.frame(matrix(NA,1,length(top_int_names)))
  names(top_int_brt)<-names(top_int_rf)<-top_int_names
  
  range.trees.names<-c("min.trees","max.trees","less1000")
  range.trees<-data.frame(matrix(NA,1,length(range.trees.names)))
  colnames(range.trees)<-range.trees.names
  
    for (i in 1:length(res.sets)){
      
      res.sets[[i]]->dat
      
      dat$runs->runs
      
      # BRT
      if (BRT_mods[i]==T){
        
        dat$n.test->brt_qual_sim_res[i,1]->lm_brt_qual_sim_res[i,1]
        if (len_sim[i]==T & dat$n.test>100) mean(dat$n.sites)->brt_qual_sim_res[i,2]->lm_brt_qual_sim_res[i,2] else dat$n.test->brt_qual_sim_res[i,2]->lm_brt_qual_sim_res[i,2]
        if (len_sim[i]==T) dat$st_leng*100->brt_qual_sim_res[i,3]->lm_brt_qual_sim_res[i,3] else 100->brt_qual_sim_res[i,3]->lm_brt_qual_sim_res[i,3]
        if (len_sim[i]==T) dat$st_level->brt_qual_sim_res[i,4]->lm_brt_qual_sim_res[i,4] else 0->brt_qual_sim_res[i,4]->lm_brt_qual_sim_res[i,4]        

        if (dat$n.test>100){
          colMeans(dat$res_brt[,c(8,10:13)], na.rm = T)->brt_qual_sim_res[i,c(c(8,10:13)+4)]
          
          colMeans(dat$res_brt[,c(1:4)], na.rm = T)->brt_mean_rank[i,]
          
          range(dat$res_brt$ntree)->range.trees[i,1:2]
          
          length(which(dat$res_brt$ntree<1000))/runs->range.trees[i,3]
          
          length(which(dat$res_brt$s1==1))/runs->brt_qual_sim_res[i,5]
          length(which(dat$res_brt$s2==2))/runs->brt_qual_sim_res[i,6]
          length(which(dat$res_brt$s3==3))/runs->brt_qual_sim_res[i,7]
          length(which(dat$res_brt$s4==4))/runs->brt_qual_sim_res[i,8]
          
          length(which(dat$res_brt$int2=="s1:s2"))/runs->brt_qual_sim_res[i,9]
          length(which(dat$res_brt$int1=="s1:s3"))/runs->brt_qual_sim_res[i,10]
          length(which(dat$res_brt$int3=="s2:s4"))/runs->brt_qual_sim_res[i,11]
          
          length(which(unlist(dat$res_brt[,5:7])=="s1:s2"))/runs->top_int_brt[i,1]
          length(which(unlist(dat$res_brt[,5:7])=="s1:s3"))/runs->top_int_brt[i,2]
          length(which(unlist(dat$res_brt[,5:7])=="s2:s4"))/runs->top_int_brt[i,3]
          
          max(dat$res_rf$max.cor)->brt_qual_sim_res[i,13]
          
          # LM
          (abs(colMeans(dat$res_lm_brt[,1:7], na.rm = T))-abs(ses))*100/abs(ses)->lm_brt_qual_sim_res[i,5:11]
          colMeans(dat$res_lm_brt[,8:9], na.rm = T)->lm_brt_qual_sim_res[i,12:13]
          }
      }
      
      # RF
      if (RF_mods[i]==T){
        dat$n.test->rf_qual_sim_res[i,1]->lm_rf_qual_sim_res[i,1]
        if (len_sim[i]==T & dat$n.test>100) mean(dat$n.sites)->rf_qual_sim_res[i,2]->lm_rf_qual_sim_res[i,2] else dat$n.test->rf_qual_sim_res[i,2]->lm_rf_qual_sim_res[i,2]
        if (len_sim[i]==T) dat$st_leng*100->rf_qual_sim_res[i,3]->lm_rf_qual_sim_res[i,3] else 100->rf_qual_sim_res[i,3]->lm_rf_qual_sim_res[i,3]
        if (len_sim[i]==T) dat$st_level->rf_qual_sim_res[i,4]->lm_rf_qual_sim_res[i,4] else 0->rf_qual_sim_res[i,4]->lm_rf_qual_sim_res[i,4]
        

        colMeans(dat$res_rf[,c(8,10:13)], na.rm = T)->rf_qual_sim_res[i,c(c(8,10:13)+4)]
        
        colMeans(dat$res_rf[,c(1:4)], na.rm = T)->rf_mean_rank[i,]
        
        length(which(dat$res_rf$s1==1))/runs->rf_qual_sim_res[i,5]
        length(which(dat$res_rf$s2==2))/runs->rf_qual_sim_res[i,6]
        length(which(dat$res_rf$s3==3))/runs->rf_qual_sim_res[i,7]
        length(which(dat$res_rf$s4==4))/runs->rf_qual_sim_res[i,8]
        
        length(which(dat$res_rf$int2=="s1:s2"))/runs->rf_qual_sim_res[i,9]
        length(which(dat$res_rf$int1=="s1:s3"))/runs->rf_qual_sim_res[i,10]
        length(which(dat$res_rf$int3=="s2:s4"))/runs->rf_qual_sim_res[i,11]
        
        length(which(unlist(dat$res_rf[,5:7])=="s1:s2"))/runs->top_int_rf[i,1]
        length(which(unlist(dat$res_rf[,5:7])=="s1:s3"))/runs->top_int_rf[i,2]
        length(which(unlist(dat$res_rf[,5:7])=="s2:s4"))/runs->top_int_rf[i,3]
        
        max(dat$res_rf$max.cor)->rf_qual_sim_res[i,13]
        # LM
        (abs(colMeans(dat$res_lm_rf[,1:7], na.rm = T))-abs(ses))*100/abs(ses)->lm_rf_qual_sim_res[i,5:11]
        colMeans(dat$res_lm_rf[,8:9], na.rm = T)->lm_rf_qual_sim_res[i,12:13]
      }
      
    }
    
    return(list(brt_qual_sim_res=brt_qual_sim_res, rf_qual_sim_res=rf_qual_sim_res, 
           lm_brt_qual_sim_res=lm_brt_qual_sim_res, lm_rf_qual_sim_res=lm_rf_qual_sim_res,
           brt_mean_rank=brt_mean_rank,rf_mean_rank=rf_mean_rank,
           top_int_brt=top_int_brt,top_int_rf=top_int_rf,range.trees=range.trees))
}

# sim_plots() runs a number of simulations to assess the effect of varying sample size
# and gradient length in multi-stressor analysis
#
# Inputs:
#
# res.mat: one of the matrices resulting from a sim_res() analysis
# var.names: stressor names
# ses_test: logical (TRUE or FALSE)
#           Are you going to perform the SES estimation error assessment?
# resp.name: Name of the variable to show (it will be assign to ylab and main in the plots)
# leg.pos: position of the legend indicating gradient length. 
#          Type help(legend) for more info about this parameter
# 
# Output:
#
# A plot with showing the specified information

sim_plots<-function(res.mat, var.names, ses_test=F, resp.name="SES error percentage",
                     leg.pos="bottomright"){
  
  res.mat[,c("n.sites","st_leng",var.names)]->res.mat
  
  na.omit(res.mat)->res.mat
  
  par(cex=1.5,cex.axis=1.5,cex.lab=1.5,mar = c(5, 5, 3, 1))
  
  # Setting graphic parameters
  if(length(var.names)==1) par(mfrow=c(1,1))
  
  if(length(var.names)==2) par(mfrow=c(1,2))
  
  if(length(var.names)>4 & length(var.names)<9) par(mfrow=c(2,4))
  for (i in var.names) {
    
    # Create Line Chart
    
    # convert factor to numeric for convenience
    classes<-sort(unique(res.mat$st_leng))
    leng_classes <- length(classes)
    
    # get the range for the x and y axis
    xrange <- range(c(25,res.mat$n.sites))
    if(ses_test==T) yrange <- range(c(0, res.mat[,i])) else yrange <- range(res.mat[,i])
    
    # set up the plot
    plot(xrange, yrange, type="n", xlab="Number of sites",
         ylab= resp.name)
    if(ses_test==T) abline(h=0,lwd=2, lty=3)
    colors <- rainbow(leng_classes)
    linetype <- c(1:leng_classes)
    plotchar <- 17+c(1:leng_classes)
    
      # add lines
    for (j in 1:leng_classes) {
      sub.dat <- subset(res.mat, st_leng==classes[j])
       lines(sub.dat$n.sites, sub.dat[,i], type="b", lwd=3,
            lty=linetype[j], col=colors[j], pch=plotchar[j],cex=2)
    }
    

    
    # add a title
    title(paste(i))
    
    if(var.names[length(var.names)]==i & length(var.names)>3) {
      plot(1, type = "n", axes=FALSE, xlab="", ylab="")
      legend(leg.pos, bty="n", legend = classes, pch=plotchar, 
             lty=linetype, title="Gradient length", col=colors, 
             lwd=3, cex=1.25, pt.cex=2, horiz = F, y.intersp = 0.5,
             xjust=0.5)
    }
    if(var.names[1]==i & length(var.names)<3) {
      legend(leg.pos, bty="n", legend = classes, pch=plotchar, 
             lty=linetype, title="Gradient length",col=colors, 
             lwd=3, cex=1.25, pt.cex=2, horiz = F, y.intersp = 0.5,
             xjust=0.5)}
  }
}