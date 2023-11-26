library(poolSeq)
library(ACER)

#I need to do this locally one chromosome at a time. This takes forever but I just want to be done. 

#read in the file
#these first few will stay the same. 
reps <- c(1:6)
#Always one generation between the pools
gen <- rep(1,6)

pops <-c('81High', '81Low', 
         '96High', '96Low',
         '418High', '418Low')

#Creating the vars for the CMH test 
rep<-c(1,1,2,2,3,3) #Number of replicates
Ne<-c(2538, 2857, 2415) #Using the population size I phenotypes 
tp<-rep(rep(c(0,1)),3) #Generations of evolution for each sample. Set to one because there is only one gen between pools. 
ps<-rep(50, 6) #Pool size


####
syncX <- read.sync(file="../Data/allShapeCoreGenome_X.sync", 
                  gen=gen, repl=reps, 
                  polarization = "minor", 
                  keepOnlyBiallelic = TRUE)

af.matX <- matrix(NA, nrow = nrow(syncX@alleles), ncol = 6)
colnames(af.matX) <- pops

for (i in 1:ncol(af.matX)){
  tempdat <- af(syncX, repl = i, gen = 1)
  af.matX[,i] <- as.matrix(tempdat)
}

af.matX <- na.omit(af.matX)
head(af.matX)
dim(af.matX)

dim(syncX@alleles)

#now to make a coverage one. 

cov.matX <- matrix(NA, nrow = nrow(syncX@alleles), ncol = 6)
# cov.mat[,1:2] <- sync@alleles[1,]
colnames(cov.matX) <- pops

for (i in 1:ncol(cov.matX)){
  tempdat <- coverage(syncX, repl = i, gen = 1)
  cov.matX[,i] <- as.matrix(tempdat)
}

crapX <- data.frame(cov.matX, syncX@alleles[,1:2])
crapX[crapX==0] <- NA
crap2X <- na.omit(crapX)
locationX <- crap2X[,7:8]

cov.matX[cov.matX==0] <- NA
cov.matX <- na.omit(cov.matX)

dim(cov.matX)

head(cov.matX)

pvalX <- adapted.cmh.test(freq=af.matX, coverage=cov.matX, 
                         Ne=Ne, gen=tp, repl=rep, poolSize=ps)

#these are all 1. So everything goes away when we account for drift?
padjX <- p.adjust(pvalX, "fdr")
#This doesn't line up with the locations from the original file because we droped those NA lines. 
afdatX <- cbind(locationX, af.matX, pvalX, padjX)

write.csv(afdatX, "../Data/allShapeMales_cmh_acer_X.csv")

rm(list=ls(pattern="X"))
#####################


sync2L <- read.sync(file="../Data/allShapeCoreGenome_2L.sync", 
                    gen=gen, repl=reps, 
                    polarization = "minor", 
                    keepOnlyBiallelic = TRUE)

af.mat2L <- matrix(NA, nrow = nrow(sync2L@alleles), ncol = 6)
colnames(af.mat2L) <- pops

for (i in 1:ncol(af.mat2L)){
  tempdat <- af(sync2L, repl = i, gen = 1)
  af.mat2L[,i] <- as.matrix(tempdat)
}

af.mat2L <- na.omit(af.mat2L)
head(af.mat2L)
dim(af.mat2L)

dim(sync2L@alleles)

#now to make a coverage one. 

cov.mat2L <- matrix(NA, nrow = nrow(sync2L@alleles), ncol = 6)
# cov.mat[,1:2] <- sync@alleles[1,]
colnames(cov.mat2L) <- pops

for (i in 1:ncol(cov.mat2L)){
  tempdat <- coverage(sync2L, repl = i, gen = 1)
  cov.mat2L[,i] <- as.matrix(tempdat)
}

crap2L <- data.frame(cov.mat2L, sync2L@alleles[,1:2])
crap2L[crap2L==0] <- NA
crap22L <- na.omit(crap2L)
location2L <- crap22L[,7:8]

cov.mat2L[cov.mat2L==0] <- NA
cov.mat2L <- na.omit(cov.mat2L)

dim(cov.mat2L)

head(cov.mat2L)

pval2L <- adapted.cmh.test(freq=af.mat2L, coverage=cov.mat2L, 
                           Ne=Ne, gen=tp, repl=rep, poolSize=ps)

#these are all 1. So everything goes away when we account for drift?
padj2L <- p.adjust(pval2L, "fdr")
#This doesn't line up with the locations from the original file because we droped those NA lines. 
afdat2L <- cbind(location2L, af.mat2L, pval2L, padj2L)

write.csv(afdat2L, "../Data/allShapeMales_cmh_acer_2L.csv")

rm(list=ls(pattern="2L"))


##########

sync2R <- read.sync(file="../Data/allShapeCoreGenome_2R.sync", 
                    gen=gen, repl=reps, 
                    polarization = "minor", 
                    keepOnlyBiallelic = TRUE)

af.mat2R <- matrix(NA, nrow = nrow(sync2R@alleles), ncol = 6)
colnames(af.mat2R) <- pops

for (i in 1:ncol(af.mat2R)){
  tempdat <- af(sync2R, repl = i, gen = 1)
  af.mat2R[,i] <- as.matrix(tempdat)
}

af.mat2R <- na.omit(af.mat2R)
head(af.mat2R)
dim(af.mat2R)

dim(sync2R@alleles)

#now to make a coverage one. 

cov.mat2R <- matrix(NA, nrow = nrow(sync2R@alleles), ncol = 6)
# cov.mat[,1:2] <- sync@alleles[1,]
colnames(cov.mat2R) <- pops

for (i in 1:ncol(cov.mat2R)){
  tempdat <- coverage(sync2R, repl = i, gen = 1)
  cov.mat2R[,i] <- as.matrix(tempdat)
}

crap2R <- data.frame(cov.mat2R, sync2R@alleles[,1:2])
crap2R[crap2R==0] <- NA
crap22R <- na.omit(crap2R)
location2R <- crap22R[,7:8]

cov.mat2R[cov.mat2R==0] <- NA
cov.mat2R <- na.omit(cov.mat2R)

dim(cov.mat2R)

head(cov.mat2R)

pval2R <- adapted.cmh.test(freq=af.mat2R, coverage=cov.mat2R, 
                           Ne=Ne, gen=tp, repl=rep, poolSize=ps)

#these are all 1. So everything goes away when we account for drift?
padj2R <- p.adjust(pval2R, "fdr")
#This doesn't line up with the locations from the original file because we droped those NA lines. 
afdat2R <- cbind(location2R, af.mat2R, pval2R, padj2R)

write.csv(afdat2R, "../Data/allShapeMales_cmh_acer_2R.csv")

rm(list=ls(pattern="2R"))


################

sync3L <- read.sync(file="../Data/allShapeCoreGenome_3L.sync", 
                    gen=gen, repl=reps, 
                    polarization = "minor", 
                    keepOnlyBiallelic = TRUE)

af.mat3L <- matrix(NA, nrow = nrow(sync3L@alleles), ncol = 6)
colnames(af.mat3L) <- pops

for (i in 1:ncol(af.mat3L)){
  tempdat <- af(sync3L, repl = i, gen = 1)
  af.mat3L[,i] <- as.matrix(tempdat)
}

af.mat3L <- na.omit(af.mat3L)
head(af.mat3L)
dim(af.mat3L)

dim(sync3L@alleles)

#now to make a coverage one. 

cov.mat3L <- matrix(NA, nrow = nrow(sync3L@alleles), ncol = 6)
# cov.mat[,1:2] <- sync@alleles[1,]
colnames(cov.mat3L) <- pops

for (i in 1:ncol(cov.mat3L)){
  tempdat <- coverage(sync3L, repl = i, gen = 1)
  cov.mat3L[,i] <- as.matrix(tempdat)
}

crap3L <- data.frame(cov.mat3L, sync3L@alleles[,1:2])
crap3L[crap3L==0] <- NA
crap23L <- na.omit(crap3L)
location3L <- crap23L[,7:8]

cov.mat3L[cov.mat3L==0] <- NA
cov.mat3L <- na.omit(cov.mat3L)

dim(cov.mat3L)

head(cov.mat3L)

pval3L <- adapted.cmh.test(freq=af.mat3L, coverage=cov.mat3L, 
                           Ne=Ne, gen=tp, repl=rep, poolSize=ps)

#these are all 1. So everything goes away when we account for drift?
padj3L <- p.adjust(pval3L, "fdr")
#This doesn't line up with the locations from the original file because we droped those NA lines. 
afdat3L <- cbind(location3L, af.mat3L, pval3L, padj3L)

write.csv(afdat3L, "../Data/allShapeMales_cmh_acer_3L.csv")

rm(list=ls(pattern="3L"))

##############


sync3R <- read.sync(file="../Data/allShapeCoreGenome_3R.sync", 
                    gen=gen, repl=reps, 
                    polarization = "minor", 
                    keepOnlyBiallelic = TRUE)

af.mat3R <- matrix(NA, nrow = nrow(sync3R@alleles), ncol = 6)
colnames(af.mat3R) <- pops

for (i in 1:ncol(af.mat3R)){
  tempdat <- af(sync3R, repl = i, gen = 1)
  af.mat3R[,i] <- as.matrix(tempdat)
}

af.mat3R <- na.omit(af.mat3R)
head(af.mat3R)
dim(af.mat3R)

dim(sync3R@alleles)

#now to make a coverage one. 

cov.mat3R <- matrix(NA, nrow = nrow(sync3R@alleles), ncol = 6)
# cov.mat[,1:2] <- sync@alleles[1,]
colnames(cov.mat3R) <- pops

for (i in 1:ncol(cov.mat3R)){
  tempdat <- coverage(sync3R, repl = i, gen = 1)
  cov.mat3R[,i] <- as.matrix(tempdat)
}

crap3R <- data.frame(cov.mat3R, sync3R@alleles[,1:2])
crap3R[crap3R==0] <- NA
crap23R <- na.omit(crap3R)
location3R <- crap23R[,7:8]

cov.mat3R[cov.mat3R==0] <- NA
cov.mat3R <- na.omit(cov.mat3R)

dim(cov.mat3R)

head(cov.mat3R)

pval3R <- adapted.cmh.test(freq=af.mat3R, coverage=cov.mat3R, 
                           Ne=Ne, gen=tp, repl=rep, poolSize=ps)

#these are all 1. So everything goes away when we account for drift?
padj3R <- p.adjust(pval3R, "fdr")
#This doesn't line up with the locations from the original file because we droped those NA lines. 
afdat3R <- cbind(location3R, af.mat3R, pval3R, padj3R)

write.csv(afdat3R, "../Data/allShapeMales_cmh_acer_3R.csv")

rm(list=ls(pattern="3R"))

##################
sync4 <- read.sync(file="../Data/allShapeCoreGenome_4.sync", 
                   gen=gen, repl=reps, 
                   polarization = "minor", 
                   keepOnlyBiallelic = TRUE)

af.mat4 <- matrix(NA, nrow = nrow(sync4@alleles), ncol = 6)
colnames(af.mat4) <- pops

for (i in 1:ncol(af.mat4)){
  tempdat <- af(sync4, repl = i, gen = 1)
  af.mat4[,i] <- as.matrix(tempdat)
}

af.mat4 <- na.omit(af.mat4)
head(af.mat4)
dim(af.mat4)

dim(sync4@alleles)

#now to make a coverage one. 

cov.mat4 <- matrix(NA, nrow = nrow(sync4@alleles), ncol = 6)
# cov.mat[,1:2] <- sync@alleles[1,]
colnames(cov.mat4) <- pops

for (i in 1:ncol(cov.mat4)){
  tempdat <- coverage(sync4, repl = i, gen = 1)
  cov.mat4[,i] <- as.matrix(tempdat)
}

crap4 <- data.frame(cov.mat4, sync4@alleles[,1:2])
crap4[crap4==0] <- NA
crap24 <- na.omit(crap4)
location4 <- crap24[,7:8]

cov.mat4[cov.mat4==0] <- NA
cov.mat4 <- na.omit(cov.mat4)

dim(cov.mat4)

head(cov.mat4)

pval4 <- adapted.cmh.test(freq=af.mat4, coverage=cov.mat4, 
                          Ne=Ne, gen=tp, repl=rep, poolSize=ps)

#these are all 1. So everything goes away when we account for drift?
padj4 <- p.adjust(pval4, "fdr")
#This doesn't line up with the locations from the original file because we droped those NA lines. 
afdat4 <- cbind(location4, af.mat4, pval4, padj4)

write.csv(afdat4, "../Data/allShapeMales_cmh_acer_4.csv")

rm(list=ls(pattern="4"))

