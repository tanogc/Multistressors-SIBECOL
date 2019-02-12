# Appendix 1. R code and functions to apply the cookbook to analyse biological responses 
#             to multiple stressors 
#
# Analysing the impact of multiple stressors on aquatic biology: cookbook with applications in R 
# C. K. Feld, P. Segurado & C. Gutierrez-Canovas
#
# Script to perform the simulations to test cookbook response to artificial datasets
# with varying sample size and gradient length
#
# Code written by C. K. Feld, P. Segurado & C. Gutierrez-Canovas 
#
# Problems with the code, send email to GutierrezCanovasC@cardiff.ac.uk

# Install the required packages, if necessary

install.packages("MuMIn") # Multi-Model Inference
install.packages("randomForestSRC") # RF
install.packages("ggRandomForests") # RF
install.packages("gbm") # BRT
install.packages("dismo") # BRT

# Setting working directory
setwd("/my_folder/")

# Loading required functions
source("simul_functions.R")

set.seed(1234)

n <- 5000 # number of sites
ac <- 3 # accuracy SD units of error
runs <- 100 # number of random runs

BRT_mod=T
RF_mod=T

# Single stressor hierarchy s1 > s2 > s3 = s4
# Interaction hierarchy 1:3 > 1:2 > 2:4
# s1, s2, s3, s4, s1:s2, s1:s3, s2:s4
ses<-c(5, 3, 2, 2, 2, 3, 1)

# RF settings
ntree <- 2000 # tress for RF

# BRT settings
lr <- 0.005

# How to generate multi-stressor data
#sim.multi.str(n, ses, ac, mod.type="opposing",plot.int=T)->sim.set
#sim.multi.str(n, ses, ac, mod.type="antagonistic",plot.int=T)->sim.set
#sim.multi.str(n, ses, ac, mod.type="synergistic",plot.int=T)->sim.set
#sim.multi.str(n, ses, ac, mod.type="additive",plot.int=T)->sim.set
#sim.multi.str(n, ses, ac, mod.type="mixed",plot.int=T)->sim.set

#### Simulations

# For the simulations, we are going to focus on mixed interactions data
sim.multi.str(n, ses, ac, mod.type="mixed",plot.int=F)->sim.set

# Testing the effect of sites on whole gradient
sim_ntest(sim.set,n.test=25,ntree, lr, BRT_mod, RF_mod, runs)->sim_25s
sim_ntest(sim.set,n.test=50,ntree, lr, BRT_mod, RF_mod, runs)->sim_50s
sim_ntest(sim.set,n.test=75,ntree, lr, BRT_mod, RF_mod, runs)->sim_75s
sim_ntest(sim.set,n.test=100,ntree, lr, BRT_mod, RF_mod, runs)->sim_100s
sim_ntest(sim.set,n.test=150,ntree, lr, BRT_mod, RF_mod, runs)->sim_150s
sim_ntest(sim.set,n.test=300,ntree, lr, BRT_mod, RF_mod, runs)->sim_300s
sim_ntest(sim.set,n.test=500,ntree, lr, BRT_mod, RF_mod, runs)->sim_500s

# Testing the effect of sites on first quartile
sim_length(sim.set,n.test=25,stressor="s1",st_leng=0.25,st_level=0,ntree, lr, BRT_mod, RF_mod, runs)->sim_q1_25s
sim_length(sim.set,n.test=50,stressor="s1",st_leng=0.25,st_level=0,ntree, lr, BRT_mod, RF_mod, runs)->sim_q1_50s
sim_length(sim.set,n.test=75,stressor="s1",st_leng=0.25,st_level=0,ntree, lr, BRT_mod, RF_mod, runs)->sim_q1_75s
sim_length(sim.set,n.test=100,stressor="s1",st_leng=0.25,st_level=0,ntree, lr, BRT_mod, RF_mod, runs)->sim_q1_100s
sim_length(sim.set,n.test=150,stressor="s1",st_leng=0.25,st_level=0,ntree, lr, BRT_mod, RF_mod, runs)->sim_q1_150s
sim_length(sim.set,n.test=300,stressor="s1",st_leng=0.25,st_level=0,ntree, lr, BRT_mod, RF_mod, runs)->sim_q1_300s
sim_length(sim.set,n.test=500,stressor="s1",st_leng=0.25,st_level=0,ntree, lr, BRT_mod, RF_mod, runs)->sim_q1_500s

# Testing the effect of sites on first + second quartile
sim_length(sim.set,n.test=25,stressor="s1",st_leng=0.50,st_level=0,ntree, lr, BRT_mod, RF_mod, runs)->sim_q2_25s
sim_length(sim.set,n.test=50,stressor="s1",st_leng=0.50,st_level=0,ntree, lr, BRT_mod, RF_mod, runs)->sim_q2_50s
sim_length(sim.set,n.test=75,stressor="s1",st_leng=0.50,st_level=0,ntree, lr, BRT_mod, RF_mod, runs)->sim_q2_75s
sim_length(sim.set,n.test=100,stressor="s1",st_leng=0.50,st_level=0,ntree, lr, BRT_mod, RF_mod, runs)->sim_q2_100s
sim_length(sim.set,n.test=150,stressor="s1",st_leng=0.50,st_level=0,ntree, lr, BRT_mod, RF_mod, runs)->sim_q2_150s
sim_length(sim.set,n.test=300,stressor="s1",st_leng=0.50,st_level=0,ntree, lr, BRT_mod, RF_mod, runs)->sim_q2_300s
sim_length(sim.set,n.test=500,stressor="s1",st_leng=0.50,st_level=0,ntree, lr, BRT_mod, RF_mod, runs)->sim_q2_500s

# Testing the effect of sites on first + second + third quartile
sim_length(sim.set,n.test=25,stressor="s1",st_leng=0.75,st_level=0,ntree, lr, BRT_mod, RF_mod, runs)->sim_q3_25s
sim_length(sim.set,n.test=50,stressor="s1",st_leng=0.75,st_level=0,ntree, lr, BRT_mod, RF_mod, runs)->sim_q3_50s
sim_length(sim.set,n.test=75,stressor="s1",st_leng=0.75,st_level=0,ntree, lr, BRT_mod, RF_mod, runs)->sim_q3_75s
sim_length(sim.set,n.test=100,stressor="s1",st_leng=0.75,st_level=0,ntree, lr, BRT_mod, RF_mod, runs)->sim_q3_100s
sim_length(sim.set,n.test=150,stressor="s1",st_leng=0.75,st_level=0,ntree, lr, BRT_mod, RF_mod, runs)->sim_q3_150s
sim_length(sim.set,n.test=300,stressor="s1",st_leng=0.75,st_level=0,ntree, lr, BRT_mod, RF_mod, runs)->sim_q3_300s
sim_length(sim.set,n.test=500,stressor="s1",st_leng=0.75,st_level=0,ntree, lr, BRT_mod, RF_mod, runs)->sim_q3_500s

# Producing final results

res.sets<-list(sim_25s,sim_50s,sim_75s,sim_100s,sim_150s,sim_300s,sim_500s,
               sim_q1_25s,sim_q1_50s,sim_q1_75s,sim_q1_100s,sim_q1_150s,sim_q1_300s,sim_q1_500s,
               sim_q2_25s,sim_q2_50s,sim_q2_75s,sim_q2_100s,sim_q2_150s,sim_q2_300s,sim_q2_500s,
               sim_q3_25s,sim_q3_50s,sim_q3_75s,sim_q3_100s,sim_q3_150s,sim_q3_300s,sim_q3_500s)

len_sim<-c(rep(F,7),rep(T,21))
BRT_mods<-rep(T,length(res.sets))
RF_mods<-rep(T,length(res.sets))

sim_res(res.sets, len_sim, ses=c(-5,-3,-2,+2,-2,+3,+1), BRT_mods, RF_mods)->overall.res

# 9.5 x 9.5
sim_plots(overall.res$lm_brt_qual_sim_res, ses_test=F, var.names="r2", resp.name="r2",leg.pos="topleft")
sim_plots(overall.res$lm_rf_qual_sim_res, ses_test=F, var.names="r2", resp.name="r2",leg.pos="bottomright")
# 14 x 7.5
sim_plots(overall.res$lm_brt_qual_sim_res, ses_test=T, var.names=sim.set$var.names, resp.name="SES error percentage")
sim_plots(overall.res$lm_rf_qual_sim_res, ses_test=T, var.names=sim.set$var.names, resp.name="SES error percentage")

write.table(overall.res$brt_mean_rank,"brt_rank.txt",sep="\t")
write.table(overall.res$rf_mean_rank,"rf_rank.txt",sep="\t")

write.table(overall.res$brt_qual_sim_res,"brt_res.txt",sep="\t")
write.table(overall.res$rf_qual_sim_res,"rf_res.txt",sep="\t")

write.table(overall.res$top_int_brt,"int_brt_res.txt",sep="\t")
write.table(overall.res$top_int_rf,"int_rf_res.txt",sep="\t")

