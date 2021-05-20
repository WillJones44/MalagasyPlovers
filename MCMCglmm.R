# http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

library(MCMCglmm)
library(ape)
library(readxl)
library(phytools)
library(tidyverse)
library(geiger)
library(Rmisc)

## get the tree

plover <- read_excel("C:/Users/willj/Dropbox/Plovers/plover manuscripts/Demography/Rfiles/SizeSurvival.xlsx", sheet = "datos")

tree<-read.newick("C:/Users/willj/Dropbox/Plovers/plover manuscripts/Demography/Rfiles/phylogeny.nwk")
data<-column_to_rownames(plover, var = "Taxa")

plot(tree)


obj<-name.check(tree,data)
obj

shorttree<-drop.tip(tree, obj$tree_not_data)
name.check(shorttree,data)
bm<-corBrownian(1, shorttree)
bm

plot(shorttree)

## load the variables

library(MCMCglmm)

mcmc <- read_excel("C:/Users/willj/Dropbox/Plovers/plover manuscripts/Demography/Rfiles/SizeSurvival.xlsx", sheet = "mcmcglmm")

mcmc<-column_to_rownames(mcmc, var = "Pop")


phylo<-shorttree
data<-mcmc

inv.phylo<-inverseA(phylo,nodes="TIPS",scale=TRUE)

prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
prior2<-list(G=list(G1=list(V=1,nu=100,alfa.mu=0,alfa.V=1)),R=list(V=1,fix=1))


# All effects in one model

model_simple<-MCMCglmm(Survival~Clutch*EquatorDist+LogMass,random=~phylo, prior = prior, family="gaussian",ginverse=list(phylo=inv.phylo$Ainv), data=data,nitt=5000000,burnin=10000,thin=5000)

summary(model_simple)

plot(model_simple)

lambda <- model_simple$VCV[,'phylo']/
  (model_simple$VCV[,'phylo']+model_simple$VCV[,'units'])
mean(lambda)
HPDinterval(lambda)
posterior.mode(lambda)






















##### OLD!#####

#Distance from equator

model_Dist<-MCMCglmm(Survival~EquatorDist,random=~phylo, prior=prior, family="gaussian",ginverse=list(phylo=inv.phylo$Ainv), data=data,nitt=500000,burnin=1000,thin=500)

summary(model_Dist)


plot(model_Dist)

Dlambda <- model_Dist$VCV[,'phylo']/
  (model_Dist$VCV[,'phylo']+model_Dist$VCV[,'units'])
mean(Dlambda)
HPDinterval(Dlambda)

#Clutch

model_Clutch<-MCMCglmm(Survival~Clutch,random=~phylo, family="gaussian",ginverse=list(phylo=inv.phylo$Ainv), data=data,nitt=500000,burnin=1000,thin=500)

summary(model_Clutch)

plot(model_Clutch)

Clambda <- model_Clutch$VCV[,'phylo']/
  (model_Clutch$VCV[,'phylo']+model_Clutch$VCV[,'units'])
mean(Clambda)
HPDinterval(Clambda)

#Mass

model_Mass<-MCMCglmm(Survival~LogMass,random=~phylo, family="gaussian",ginverse=list(phylo=inv.phylo$Ainv), data=data,nitt=500000,burnin=1000,thin=500)

summary(model_Mass)

plot(model_Mass)

Mlambda <- model_Mass$VCV[,'phylo']/
  (model_Mass$VCV[,'phylo']+model_Mass$VCV[,'units'])
mean(Mlambda)
HPDinterval(Mlambda)




