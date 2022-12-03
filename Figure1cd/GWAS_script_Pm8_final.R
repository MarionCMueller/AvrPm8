##########################################################
### Script to run GWAS Analysis for AvrPm8 ###############
##########################################################

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("https://zzlab.net/GAPIT/previous/gapit_functions_20220411.txt")

myG<-read.table("SNPcall_79_onISR7_hapmap.txt",h=F)
myY<-read.table("Phenotypes_Pm8.txt",h=T)
dir.create("GWAS_79isolates")
setwd("GWAS_79isolates")
GAPIT1<-GAPIT(G=myG,Y=myY,PCA.total =3, Model.selection = TRUE,model="GLM",kinship.algorithm = "VanRanden")
