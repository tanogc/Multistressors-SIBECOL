## Code to run the multi-stressor analysis

# Setting working directory
setwd("")

# Loading required libraries
# The required packages should be installed

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("variancePartition", version = "3.8")
install.packages("usdm")
install.packages("randomForestSRC")
install.packages("ggRandomForests")
install.packages("gbm")
install.packages("dismo")
install.packages("MuMIn")
install.packages("lattice")

library(usdm) # Collinearity
library(lattice) # dotplots
library(randomForestSRC) # RF
library(ggRandomForests) # RF
library(gbm) # BRT
library(dismo) # BRT
library(MuMIn) # Multi-model inference
library(variancePartition)

# Loading required functions
source("simul_functions.R")

# Simulated dataset

set.seed (1234) # sets a numerical starting point 
n <- 100 # number of sites
ac <- 3 # accuracy SD units of error

# Single stressor hierarchy s1 > s2 > s3 = s4
# Interaction hierarchy 1:3 > 1:2 > 2:4
# s1, s2, s3, s4, s1:s2, s1:s3, s2:s4
ses<-c(5, 3, 2, 2, 2, 3, 1)

# Simulating data
sim.multi.str(n, ses, ac, mod.type="additive",plot.int=T)->sim.set.add
sim.multi.str(n, ses, ac, mod.type="antagonistic",plot.int=T)->sim.set.ant
sim.multi.str(n, ses, ac, mod.type="synergistic",plot.int=T)->sim.set.syn
sim.multi.str(n, ses, ac, mod.type="opposing",plot.int=T)->sim.set.opo

sim.multi.str(n, ses, ac, mod.type="mixed",plot.int=F)->sim.set

sim.set$sim.dat->dat

###########################
# Variable transformation
###########################

q<-sapply(dat,class)=="numeric" | sapply(dat,class)=="integer"# selecting quantitative variables

par(mfrow=c(3,3))
for (i in which(q==T)) hist(dat[,i], main=names(dat)[i])
par(mfrow=c(1,1))

# Collinearity

# calculates pairwise Pearson correlation coefficients for all variables 
# of the object dat; note that the function 	is only applicable to numerical variables

round(as.dist(cor(dat[,-1])),2)

# Exploring collinearity using Variance Inflation Factor
vifstep (dat[,-1], th=7) # threshold set to VIF<7

###########################
# Exploratory analysis
###########################

cor.res <- round(cor(dat[,-1],dat$y),2)

par(mfrow=c(3,3))

for (i in which(q==T)) 
  if(i>1) {
    plot(dat$y~dat[,i], main=names(dat)[i], ylab="y",xlab=names(dat)[i])
    mod <- lm(dat$y ~ dat[, i])
    abline(mod, lwd=3, col="blue")
    cor1 <- bquote(italic(r) == .(round(cor.res[i-1], 3)))
    mtext(cor1, side = 3, adj = 0.975, line = -10.5, cex = 1)
  }

par(mfrow=c(1,1))


###########################
# Random Forsest (RF)
###########################

# sets a numerical starting point, which will be set randomly if not set by the user 
set.seed (1234)
 
# How many variables per node? (mtry)
floor(ncol(dat[,-1])/3)

# function to run RF my.rf # output of model details, e.g. goodness-of-fit, out-of-bag (OBB) error
my.rf <- rfsrc (y ~ ., mtry = 6, ntree = 2000, importance="permute", 	data = dat) 

my.rf

# For regression we set mtry as upper integer resulting from the number of predictors / 3

# plot to determine the optimum number of trees
plot (na.omit(gg_error (my.rf)) )

# gg_vimp() provides the predictor's importance
my.rf.vimp <- gg_vimp (my.rf)
my.rf.vimp
plot (my.rf.vimp) #plot variable importance

# Saving plot with variable importance

# Dotplot of the variable importance for the RF model
pdf(file="rf_varimp.pdf",onefile=T,width=10,height=10)

my.rf.vimp$vimp->vimp
my.rf.vimp$vars->names(vimp)

dotplot(sort(vimp),xlab = list("importance",cex = 1.75),par.settings = list(axis.text = list(cex = 2), 
                                                                            par.xlab.text = list(cex = 2), 
                                                                            par.ylab.text = list(cex = 2)),
        panel = function(x,y){ 
          panel.dotplot(x, y, pch=16, cex=2,col="black") 
          panel.abline(v=0, col="grey", lty="longdash", lwd=2) 
          
        } 
) 
dev.off()

# plots the marginal response for each predictor
my.rf.part.plot <- plot.variable(my.rf, partial=TRUE, sorted=F, show.plots = FALSE )
gg.part <- gg_partial(my.rf.part.plot)
plot(gg.part, xvar = names(dat[,-1]), panel=TRUE, se=T)

# now adding smooth.lines
plot.variable (my.rf, partial = TRUE, smooth.lines = TRUE) 

# Selecting the best subset of predictors
md.obj <- max.subtree (my.rf)
# extracts the names of the variables in object md.obj
md.obj$topvars 

# method = "vimp" 	computes the additive and joint effects of each pair 
# of stressors in 	the response variable
my.rf.interaction <- find.interaction (my.rf, xvar.names = md.obj$topvars, 	
                                       importance= "permute", method = "vimp", 
                                       nrep = 3) 

###########################
# Boosted Regression Trees (BRT)
###########################

# BRT with 2-way interactions and 67% of the data used to train the model
my.brt <- gbm.step (data=dat, gbm.x=c(2:ncol(dat)), gbm.y=1,  	
                    family = "gaussian", tree.complexity = 2, 
                    learning.rate = 0.01, 	bag.fraction = 0.67) 

# Explained deviance based on CV
my.brt$self.statistics$mean.null->null.dev
my.brt$cv.statistics$deviance.mean->resid.dev

1-resid.dev/null.dev

# bar plot of predictor's contributions to the deviance explained in the model
brt.imp <- summary (my.brt)
brt.imp

# Dotplot of the variable importance for the BRT model
pdf(file="brt_varimp.pdf",onefile=T,width=10,height=10)

brt.imp$rel.inf->vimp
brt.imp$var->names(vimp)

dotplot(sort(vimp),xlab = list("importance",cex = 1.75),par.settings = list(axis.text = list(cex = 2), 
                                                                            par.xlab.text = list(cex = 2), 
                                                                            par.ylab.text = list(cex = 2)),
        panel = function(x,y){ 
          panel.dotplot(x, y, pch=16, cex=2,col="black") 
          panel.abline(v=0, col="grey", lty="longdash", lwd=2) 
          
        } 
) 
dev.off()

# Partial dependence plot with a fitted line and a smoother line overlaid, 
# y-axis normalised
plot.my.brt <- gbm.plot (my.brt, smooth=TRUE, n.plots=9, write.title = FALSE, 	plot.layout=c(3, 3)) 

# Simplifying BRT model
my.brt.simp <- gbm.simplify (my.brt) 

int.my.brt <- gbm.interactions (my.brt) # creates an object with the results
int.my.brt$interactions # shows interactions
int.my.brt$rank.list # shows a list of ranked interactions

# z.range defines the range of the z-axis, theta y phi parametercs set perspective 
# i.e. between -20 and 20 in this example
gbm.perspec (my.brt, 1, 2, z.range=c(-20, 20), theta = 20, phi = 30) 

###########################
# Multi-model inference
###########################

# Collinearity

# calculates pairwise Pearson correlation coefficients for all variables 
# of the object dat; note that the function 	is only applicable to numerical variables

round(as.dist(cor(dat[,c("s1", "s2", "s3", "s6", "s9")])),2)

# Exploring collinearity using Variance Inflation Factor
vifstep (dat[,c("s1", "s2", "s3", "s6", "s9")], th=7) # threshold set to VIF<7

# Global model
mod <- lm (y ~ s1 + s2 + s3 + s6 + s9 + s1:s2 + s1:s3, data=dat)
summary(mod)

r.squaredGLMM(mod)

round(calcVarPart(mod),2)->mod.varPart
mod.varPart

# Dotplot of the variable importance for the generalised linear model
pdf(file="brt_varimp.pdf",onefile=T,width=10,height=10)

dotplot(sort(mod.varPart[-length(mod.varPart)]),xlab = list("importance",cex = 1.75),par.settings = list(axis.text = list(cex = 2), 
                                                                            par.xlab.text = list(cex = 2), 
                                                                            par.ylab.text = list(cex = 2)),
        panel = function(x,y){ 
          panel.dotplot(x, y, pch=16, cex=2,col="black") 
          panel.abline(v=0, col="grey", lty="longdash", lwd=2) 
          
        } 
) 
dev.off()

options(na.action = "na.fail") # necessary to run dredge()

# runs all possible models 	and ranks the output according to the AIC, R2 is printed 
mod_d <- dredge (mod, rank = "AIC", 
                      extra = c(R2=function(x) r.squaredGLMM(x)))
mod_d # prints the model ranking

mod_set <- get.models (mod_d, subset=delta<=2) # subset 	with delta AIC ≤2 
mod_set <- get.models (mod_d, subset = cumsum  	(mod_d$weight)<=.95) # subset of top models based on AICw≤0.95

mod_av <- model.avg (mod_set, revised.var = TRUE) # runs model 	averaging 
summary (mod_av) # summary statistics of the model averaging

d.glm.df<-data.frame(mod_d) # storing models
length(mod_set)->n.mod # number of models in the subset

# Storing the top models in d.res
d.res<-data.frame(d.glm.df[1:n.mod,])

weighted.mean(mod_d$R2.R2m[1:n.mod],mod_d$weight[1:n.mod])->r2m.w # weighted r2m
weighted.mean(mod_d$R2.R2c[1:n.mod],mod_d$weight[1:n.mod])->r2c.w # weighted r2c

r2m.w
r2c.w

# Exporting the top models to a text file
write.table(d.res,"d.res.txt",sep="\t",row.names=F)

# This function plot model assumptions (residual normality and homogeneity of variances)
# Plots can be stored in any folder. Default folder is the working directory

for (i in 1:n.mod) 
{
  mod_set[[i]]->mod
  pdf(file=paste("mod",i,"_assum",".pdf",sep=""),onefile=T,width=6.24,height=4.05)
  par(mfrow=c(1,2),cex=1, cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
  m.resid<-resid(mod)
  hist(m.resid,main="",xlab="residuals")
  mtext(paste("model",i),line=0.5,at=4,cex=2)
  plot(fitted(mod),m.resid,xlab="fitted values",ylab="residuals",main="")  
  dev.off()
}


# Plotting model results
plot(spp.ric~env.spp$cond, ylab="",xlab=expression(paste("Conductivity (",mu,S," ",cm^-1,")",sep="")), xaxt = "n", main="",col=c("black","black","black")[env.spp2$hydro], pch=c(19,2,3)[env.spp2$hydro])
axis(1,x.labs,round(exp(back_st(x.labs,c.m,c.sd)),0))

par(mfrow=c(2,2), cex=1.5, cex.lab=1.8, cex.axis=1.5, pch=19, mar=c(5,5,3,1))

# Creating the sequence of predictor values to represent the fitted values
s1.seq<-seq(min(dat$s1),max(dat$s1),length.out=1000)
s2.seq<-seq(min(dat$s2),max(dat$s2),length.out=1000)
s3.seq<-seq(min(dat$s3),max(dat$s3),length.out=1000)

lwd.set<-6
col.int<-rainbow(3)

##################################################
# Scatterplot showing the relationship between y and s1 x s2
##################################################

plot(dat$y~dat$s1, xlab="s1", ylab="y")
mtext("a)", line = 1.0, adj = -0.2, cex = 3, font = 2)

#plotting legend in the graph
legend("topright", legend = c("s2 low","s2 intermediate","s2 high"),
       bty = "n", col = col.int,lwd=6,cex=1.5)


# Fitted values, s2 at percentile 10th
lines(s1.seq,predict(mod_av,data.frame(s1=s1.seq,
                                    s2=as.numeric(quantile(dat$s2,0.10)),
                                    s3=mean(dat$s3),
                                    s6=mean(dat$s6),
                                    s9=mean(dat$s9)),
                                    re.form=NA,type="response"),
                                    lwd=lwd.set, col=col.int[1])

# Fitted values, s2 at percentile 50th
lines(s1.seq,predict(mod_av,data.frame(s1=s1.seq,
                                       s2=as.numeric(quantile(dat$s2,0.50)),
                                       s3=mean(dat$s3),
                                       s6=mean(dat$s6),
                                       s9=mean(dat$s9)),
                                       re.form=NA,type="response"),
                                       lwd=lwd.set, col=col.int[2])

# Fitted values, s2 at percentile 90th
lines(s1.seq,predict(mod_av,data.frame(s1=s1.seq,
                                       s2=as.numeric(quantile(dat$s2,0.90)),
                                       s3=mean(dat$s3),
                                       s6=mean(dat$s6),
                                       s9=mean(dat$s9)),
                                       re.form=NA,type="response"),
                                       lwd=lwd.set, col=col.int[3])

##################################################
# Scatterplot showing the relationship between y and s1 x s3
##################################################

plot(dat$y~dat$s1, xlab="s1", ylab="y")
mtext("b)", line = 1.0, adj = -0.2, cex = 3, font = 2)

#plotting legend in the graph
legend("topright", legend = c("s3 low","s3 intermediate","s3 high"),
       bty = "n", col = col.int,lwd=6,cex=1.5)


# Fitted values, s3 at percentile 10th
lines(s1.seq,predict(mod_av,data.frame(s1=s1.seq,
                                       s2=mean(dat$s2),
                                       s3=as.numeric(quantile(dat$s3,0.10)),
                                       s6=mean(dat$s6),
                                       s9=mean(dat$s9)),
                                       re.form=NA,type="response"),
                                       lwd=lwd.set, col=col.int[1])

# Fitted values, s3 at percentile 50th
lines(s1.seq,predict(mod_av,data.frame(s1=s1.seq,
                                       s2=mean(dat$s2),
                                       s3=as.numeric(quantile(dat$s3,0.50)),
                                       s6=mean(dat$s6),
                                       s9=mean(dat$s9)),
                                       re.form=NA,type="response"),
                                       lwd=lwd.set, col=col.int[2])

# Fitted values, s3 at percentile 90th
lines(s1.seq,predict(mod_av,data.frame(s1=s1.seq,
                                       s2=mean(dat$s2),
                                       s3=as.numeric(quantile(dat$s3,0.90)),
                                       s6=mean(dat$s6),
                                       s9=mean(dat$s9)),
                                       re.form=NA,type="response"),
                                       lwd=lwd.set, col=col.int[3])
