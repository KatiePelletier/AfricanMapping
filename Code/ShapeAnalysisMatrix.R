#Here I am taking a look at the directions of variation in (two, for now) populations, because I want to compare within and between. Also since there are problems with the genetics I need something to make this a paper. 

library(geomorph)
library(emmeans)
library(tidyverse)
library(effects)
library(lme4)
library(cowplot)
library(magick)
library(boot)


projFunction <- function(x, y) {
  scalarProj <- (x %*% y) / norm(y, type = "2")
  return(scalarProj)
}

load("../Data/allMalesForThesisAnalysis.rda")


ha.m.lm <- colMeans(pure.male[pure.male$population == "HA",1:30])
la.m.lm <- colMeans(pure.male[pure.male$population == "LA",1:30])

male.diff.lm <- ha.m.lm - la.m.lm

#and one thats size corrected. 
pure.male2 <- data.frame(pure.male, lm(as.matrix(pure.male[,1:30]) 
                                       ~ logCS, data = pure.male)$residuals)

ha.m.noAllo <- colMeans(pure.male2[pure.male2$population == "HA",37:66])
la.m.noAllo <- colMeans(pure.male2[pure.male2$population == "LA",37:66])

male.diff.noAllo <- ha.m.noAllo - la.m.noAllo

###################################################################
#a PCA of both crosses together. 

#first wihtout size included. 
pca.forCross <- prcomp(cross.male[,1:30])
summary(pca.forCross)

cross.pca <- data.frame(cross.male, pca.forCross$x)

ggplot(cross.pca, aes(x = PC1, y = PC2, col = line)) + 
  geom_point(alpha = 0.4)

#unexpected  
ggplot(cross.pca, aes(x = CS, y = PC1, col = line)) + 
  geom_point(alpha = 0.4)

png("../Figures/zi192CrossPCA12.png")
ggplot(cross.pca, aes(x = PC1, y = PC2, col = line)) + 
  geom_point(alpha = 0.4) + 
  theme_classic()
dev.off()


#Now I want to take out size. I am going to model it here because this at least seems to work. 

cross.male2 <- data.frame(cross.male, lm(as.matrix(cross.male[,1:30]) 
                                         ~ logCS, data = cross.male)$residuals)

pca.cross.allo <- prcomp(cross.male2[,37:66])
summary(pca.cross.allo)

cross.pca.noAllo <- data.frame(cross.male2, pca.cross.allo$x)

#a lot of the diffrences were just size related. 
ggplot(cross.pca.noAllo, aes(x = PC1, y = PC2, col = line)) + 
  geom_point(alpha = 0.3)

#one is just bigger than the other. neat. 
ggplot(cross.pca.noAllo, aes(x = logCS, y = PC1, col = line)) + 
  geom_point(alpha = 0.3)

png("../Figures/comparingShapesPCA.png")
ggplot(cross.pca.noAllo, aes(x = PC1, y = PC2, col = line)) + 
  geom_point(alpha = 0.3) + 
  theme_classic()
dev.off()

#a PCA of both crosses together. 

#first wihtout size included. 
pca.forCross <- prcomp(cross.male[,1:30])
summary(pca.forCross)

cross.pca <- data.frame(cross.male, pca.forCross$x)

ggplot(cross.pca, aes(x = PC1, y = PC2, col = line)) + 
  geom_point(alpha = 0.4)

#unexpected  
ggplot(cross.pca, aes(x = CS, y = PC1, col = line)) + 
  geom_point(alpha = 0.4)

png("../Figures/zi192CrossPCA12.png")
ggplot(cross.pca, aes(x = PC1, y = PC2, col = line)) + 
  geom_point(alpha = 0.4) + 
  theme_classic()
dev.off()


#Now I want to take out size. I am going to model it here because this at least seems to work. 

cross.male2 <- data.frame(cross.male, lm(as.matrix(cross.male[,1:30]) 
                                         ~ logCS, data = cross.male)$residuals)

pca.cross.allo <- prcomp(cross.male2[,37:66])
summary(pca.cross.allo)

cross.pca.noAllo <- data.frame(cross.male2, pca.cross.allo$x)

#a lot of the diffrences were just size related. 
ggplot(cross.pca.noAllo, aes(x = PC1, y = PC2, col = line)) + 
  geom_point(alpha = 0.3)

#one is just bigger than the other. neat. 
ggplot(cross.pca.noAllo, aes(x = logCS, y = PC1, col = line)) + 
  geom_point(alpha = 0.3)

png("../Figures/comparingShapesPCA.png")
ggplot(cross.pca.noAllo, aes(x = PC1, y = PC2, col = line)) + 
  geom_point(alpha = 0.3) + 
  theme_classic()
dev.off()



##################Making a plot comparing the shapes in PC space with and without allometry, for a supplement#####################

lm.sizePC1temp <- ggplot(cross.pca, aes(x = logCS, y = PC1, col = line)) + 
  geom_point(alpha = 0.5) + 
  scale_colour_manual(values=c('black', 'grey57', 'grey86')) +
  theme_classic() + 
  guides(col=guide_legend("Genotype"))

lm.sizePC1 <- ggplot(cross.pca, aes(x = logCS, y = PC1, col = line)) + 
  geom_point(alpha = 0.5) + 
  scale_colour_manual(values=c('black', 'grey57', 'grey86')) +
  theme_classic() +
  theme(legend.position = "none")

lm.PC1PC2 <- ggplot(cross.pca, aes(x = PC1, y = PC2, col = line)) + 
  geom_point(alpha = 0.5) + 
  scale_colour_manual(values=c('black', 'grey57', 'grey86')) +
  theme_classic() +
  theme(legend.position = "none")

legend <- get_legend(
  # create some space to the left of the legend
  lm.sizePC1temp + theme(legend.box.margin = margin(0, 0, 0, 12))
)

allo.sizePC1 <- ggplot(cross.pca.noAllo, aes(x = logCS, y = PC1, col = line)) + 
  geom_point(alpha = 0.5) + 
  scale_colour_manual(values=c('black', 'grey57', 'grey86')) +
  theme_classic() +
  theme(legend.position = "none")

allo.PC1PC2 <- ggplot(cross.pca.noAllo, aes(x = PC1, y = PC2, col = line)) + 
  geom_point(alpha = 0.5) + 
  scale_colour_manual(values=c('black', 'grey57', 'grey86')) +
  theme_classic() +
  theme(legend.position = "none")


pcas_plot <- plot_grid(lm.sizePC1, lm.PC1PC2, allo.sizePC1, allo.PC1PC2, 
                       nrow = 2, ncol = 2, 
                       labels = c("A", "B", "C", "D"))

final.pcs.plot <- plot_grid(pcas_plot, legend, rel_widths = c(3, 0.4))


png("../Figures/zi192Pools_pca_withAndWithoutSize.png", width =2000, height = 1000, units = "px",res = 200)

final.pcs.plot

dev.off()

#################################################

#Now I want to put these into geomorph to use pairwise to compare things using pairwise
cord <- as.matrix(all.male[,1:30])
shape <- arrayspecs(cord, 15, 2)

gdf.cross.males <- geomorph.data.frame(shape = shape,
                                       CS = all.male$CS, 
                                       line = all.male$line, 
                                       logCS = all.male$logCS, 
                                       logCS_c = all.male$logCS_c)


#first I need a model and a null. this model is actually pretty simple. 

cross.males.fullMod <- procD.lm(shape ~ logCS*line, 
                                data = gdf.cross.males)


#setting up the nulls. 


# cross.males.lineNull <- procD.lm(shape ~ line, 
#                                 data = gdf.cross.males)

cross.males.morphoNull <- procD.lm(shape ~ logCS + line:logCS, 
                                   data = gdf.cross.males)

#Cs:line is actually the allo null, as pointed out by ID
# cross.males.alloNull <- procD.lm(shape ~ logCS + line, 
#                                  data = gdf.cross.males)



summary(cross.males.fullMod)

anova(cross.males.fullMod)


morpho.pair <- pairwise(cross.males.fullMod, cross.males.morphoNull, 
                        groups = gdf.cross.males$line)


#diffrent mean shapes. Makes sense. 
summary(morpho.pair)
#compare with the parents distances. 
# Pairwise distances between means, plus statistics
#                               d  UCL (95%)           Z Pr > d
# EF43:Zi192Ef43      0.01672716 0.01193464  3.80809502  0.001
# EF81:Zi192Ef81      0.02296323 0.03311468 -0.17270440  0.574
# EF96:Zi192Ef96      0.02649739 0.03456019  0.41547123  0.355
# ZI192:Zi192Ef43     0.02504553 0.01944099  3.92231252  0.001
# ZI192:Zi192Ef81     0.01742393 0.01333442  3.71573742  0.001
# ZI192:Zi192Ef96     0.01997683 0.01923898  1.95224678  0.027
# Zi192Ef43:Zi192Ef81 0.01084858 0.01160107 -0.41631695  0.671
# Zi192Ef43:Zi192Ef96 0.01455783 0.01523553 -0.31750958  0.621
# Zi192Ef81:Zi192Ef96 0.01301968 0.01346018  0.06716516  0.466

#Similar direction 
#summary(morpho.pair, test.type = "VC")

#I don't really know how to consider this one. Is this a lot or a little? maybe diffrent mags 
summary(morpho.pair, test.type = "DL")

#again, not sure if this is meaningful. Will need more data for this. 
summary(morpho.pair, test.type = "var")
#compare to the parents again to compare genetic to enviormental. 

###########################################################
#Ian recomended these:
#evolqg and vcvComp to get functions to compare matracies. 
#vcvComp has an actual vingette so I'm using that first 

library(vcvComp)

#What I think I want to ask about is the relative variance in the cross vs the parents because I think that is an interesting question for what I want. 

#this is following from the example in the vcvComp vignette. 

#step 1 is a PCA, I might have this elsewhere but honestly, I want to just keep this analysis together. 
#I am doing this with size left in. I don't think that I want this but I want to talk about the best why to do this today. 
males.lm.size.resid <- lm(as.matrix(all.male.zi192[,1:30])~ logCS, 
                          data = all.male.zi192)$residuals


all.male.zi192.pcs <- prcomp(males.lm.size.resid)$x[,1:26]

#Now to compute cov matricies 
geno.cov <- cov.group(all.male.zi192.pcs, groups = all.male.zi192$line)

#comparing variance between 
eigen.phen <- mat.sq.dist(geno.cov , dist. = "Riemannian")  # Riemannian distances
prcoa <- pr.coord(eigen.phen)  # ordination
#That is a lot of between group variance. 
prcoa$Variance 

plot(prcoa$Variance$eigenvalues)

shit <- data.frame(prcoa$PCoords)
shit$name <- rownames(shit)

#After I did this and accounted for size, this is actually sort of cool. 
ggplot(shit, aes(x = PCo1, y = PCo2, label = name)) + 
  geom_point() + 
  geom_text(vjust = 1)

#need numbers 
with(all.male.zi192, table(line))

#comparing each cross to the eth parent:
#3.915323e-34 but grain of salt with this number because of the problems with sample sizes 
prop.vcv.test(n = c(1185,  53), geno.cov[,,"Zi192Ef43"], geno.cov[,,"EF43"])  # ML test

#2.867065e-43
prop.vcv.test(n = c( 1288 ,  45 ), geno.cov[,,"Zi192Ef81"], geno.cov[,,"EF81"])  

#8.896229e-30
prop.vcv.test(n = c(880,   56 ), geno.cov[,,"Zi192Ef96"], geno.cov[,,"EF96"]) 

#now across all six 
#swap this?
relEigen.ef43 <- relative.eigen(geno.cov[,,"Zi192Ef43"], geno.cov[,,"EF43"])
relEigen.ef43$relValues  # relative eigenvalues


#using this base code, but I can make pretty later.
#This is totally unsurprising because there is so many more wings in the cross relative to the parent. 

plot(relEigen.ef43$relValues[1:relEigen.ef43$q], 
     log = "y",  las = 1, col = "blue", type = "b", 
     main = "zi192ef43 relative to ef43", cex = 0.8, 
     cex.main = 1, cex.axis = 0.8, cex.sub = 0.7, 
     sub = paste("Relative generalized variance =", relEigen.ef43$relGV), 
     xlab = NA, ylab = "Relative eigenvalues")
abline(h = 1)

#Idea: could I sub sample the cross to get about an equal number of observations to the parent, measure the relative variance and then bootstrap this? This would give me more meaningful variance estimates I think. Or I need to transform this is someway that handles the sample size. Ask ID

#What if I compare the F20s to eachother 
#These are very small numbers and suprising to me. 

#1.045531e-198
prop.vcv.test(n = c(1185,  1288), geno.cov[,,"Zi192Ef43"], geno.cov[,,"Zi192Ef81"])  # ML test

#2.162189e-204
prop.vcv.test(n = c(1288,  880), geno.cov[,,"Zi192Ef81"], geno.cov[,,"Zi192Ef96"])

#1.431902e-280
prop.vcv.test(n = c(1185,  880), geno.cov[,,"Zi192Ef43"], geno.cov[,,"Zi192Ef96"]) 

#now across all six 
#swap this?
relEigen.cross1 <- relative.eigen(geno.cov[,,"Zi192Ef43"], geno.cov[,,"Zi192Ef81"])
relEigen.cross1$relValues  # relative eigenvalues

plot(relEigen.cross1$relValues[1:relEigen.cross1$q], 
     log = "y",  las = 1, col = "blue", type = "b", 
     cex.main = 1, cex.axis = 0.8, cex.sub = 0.7, 
     sub = paste("Relative generalized variance =", relEigen.ef43$relGV), 
     xlab = NA, ylab = "Relative eigenvalues")
abline(h = 1)

relEigen.cross2 <- relative.eigen(geno.cov[,,"Zi192Ef81"], geno.cov[,,"Zi192Ef96"])
relEigen.cross2$relValues  # relative eigenvalues

plot(relEigen.cross2$relValues[1:relEigen.cross2$q], 
     log = "y",  las = 1, col = "blue", type = "b", 
     cex.main = 1, cex.axis = 0.8, cex.sub = 0.7, 
     sub = paste("Relative generalized variance =", relEigen.ef43$relGV), 
     xlab = NA, ylab = "Relative eigenvalues")
abline(h = 1)


relEigen.cross3 <- relative.eigen(geno.cov[,,"Zi192Ef43"], geno.cov[,,"Zi192Ef96"])
relEigen.cross3$relValues  # relative eigenvalues

relEigen.parent1 <- relative.eigen(geno.cov[,,"ZI192"], geno.cov[,,"EF96"])
relEigen.parent1$relValues  # relative eigenvalues

relEigen.parent2 <- relative.eigen(geno.cov[,,"ZI192"], geno.cov[,,"EF81"])
relEigen.parent2$relValues  # relative eigenvalues

relEigen.parent3 <- relative.eigen(geno.cov[,,"ZI192"], geno.cov[,,"EF43"])
relEigen.parent3$relValues  # relative eigenvalues

relEigen.parent4 <-  relative.eigen(geno.cov[,,"EF43"], geno.cov[,,"EF81"])
relEigen.parent5 <-  relative.eigen(geno.cov[,,"EF43"], geno.cov[,,"EF96"])
relEigen.parent6 <-  relative.eigen(geno.cov[,,"EF81"], geno.cov[,,"EF96"])

relEigen.parCross1 <-  relative.eigen(geno.cov[,,"Zi192Ef43"], geno.cov[,,"EF43"])

relEigen.parCross2 <-  relative.eigen(geno.cov[,,"Zi192Ef81"], geno.cov[,,"EF43"])
relEigen.parCross3 <-  relative.eigen(geno.cov[,,"Zi192Ef81"], geno.cov[,,"EF43"])


plot(relEigen.cross3$relValues[1:relEigen.cross3$q], 
     log = "y",  las = 1, col = "blue", type = "b", 
     cex.main = 1, cex.axis = 0.8, cex.sub = 0.7, 
     sub = paste("Relative generalized variance =", relEigen.ef43$relGV), 
     xlab = NA, ylab = "Relative eigenvalues")
abline(h = 1)

eigens.comp <- data.frame(relEigen.cross1$relValues, relEigen.cross2$relValues, relEigen.cross3$relValues, relEigen.parent1$relValues, relEigen.parent2$relValues, relEigen.parent3$relValues, relEigen.parent4$relValues, relEigen.parent5$relValues, relEigen.parent6$relValues, relEigen.parCross1$relValues, relEigen.parCross2$relValues, relEigen.parCross3$relValues)

dim(eigens.comp)

names(eigens.comp) <- c("comp.cross43.81", "comp.cross81.96", "comp.cross43.96", "comp.192-43", "comp.192-81", "comp.192-96", "comp.43-81", "comp.43-96", "comp.81-96", "comp.cross-pure43", "comp.cross-pure81", "comp.cross-pure96")
eigens.comp$index <- seq(1, nrow(eigens.comp), by = 1)
plot.crap <- pivot_longer(eigens.comp,  cols = starts_with("comp"), names_to = "ID", values_to = "eigen")

#I think this is telling me that the variation for the parents is distributed really diffrently between the zambian and the ethiopian parent. and that the crosses have less (and similar) variation pattrens. 
ggplot(plot.crap, aes(x = index, y = eigen, col = ID)) + 
  geom_point() + 
  geom_hline(yintercept = 1)


#More than what I did above, I think I want to compare between and within genotype changes overall? Again this probably will be meaninless. 

B <- cov.B(all.male.zi192.pcs, groups = all.male.zi192$line) #between
W <- cov.W(all.male.zi192.pcs, groups = all.male.zi192$line) #within

#number of groups and indiv. 
#0 lol. This was never going to work 
prop.vcv.test(n = c(7, 3567), B, W)

#getting the magnitude. Going ahead blindly is the best way to suceed. 
Bsc <- B / scaling.BW(B, W)  # scale B to W

# Create an array of group covariance matrices, B and W
S.bw <- array(c(geno.cov, W, Bsc), 
              dim = c(dim(geno.cov)[[1]], 
                      dim(geno.cov)[[2]], 
                      dim(geno.cov)[[3]] + 2))
dimnames(S.bw) <- list(dimnames(geno.cov)[[1]], 
                       dimnames(geno.cov)[[2]], 
                       c(dimnames(geno.cov)[[3]], "W", "B"))
# Ordination
eigen.phen.bw <- mat.sq.dist(S.bw, dist. = "Riemannian")
prcoa.bw <- pr.coord(eigen.phen.bw) #this doesnt work because I think there are 0s

# Relative PCA of B with respect to W
relEigenBW <- relative.eigen(B, W)
#Lol this is a terrible idea. Need to fix things.
#not suprising that there is much more variation between than within genotypes. 
relEigenBW$relValues  # relative eigenvalues

#other things to look at will get the distances between cov matrices. 
#I think this will be equivlent to what I've done before. 
#euclidean.dist 

#I also looked at evolQG and couldn't find a function to do the coef of variation here. Talk more to ID asap. 
#ID gave me some code he wrote for Brandon that uses evolqg. 

library(evolqg)

#this just makes the VcoV matrix for each line. 
#I can't find where these levels are set. So I will just rest now.
all.male2 <- all.male

all.male2$line <- factor(all.male2$line, levels = c("ZI192", "ZI251", "ZI418", "Zi192Ef43", "Zi192Ef81", "Zi192Ef96", "EF43", "EF81", "EF96"))

levels(all.male$line)
levels(all.male2$line)

VCV_mats <- with(all.male, 
                 by(all.male[, 1:30], 
                    line,
                    function(x) cov(x)))

#I want to size correct these estimates. Beccause fitting within line, only need size here. 
allo.mod <- lm(as.matrix(all.male[,1:30]) ~ all.male$logCS_c)
one.mat <- CalculateMatrix(allo.mod)

#this is fitting with logCS_c
VCV_mats2 <- with(all.male2, 
                  by(all.male2[, c(1:30, 36)], 
                     line,
                     function(x) CalculateMatrix(lm(as.matrix(x[,1:30]) ~ 
                                                      x[,31]))))

#Smaller. which makes sense to me. Becuase size is a major componet. 
CalcEigenVar(VCV_mats2$EF43, sd = TRUE, rel = F) # uses a correction factor
CalcEigenVar(VCV_mats$EF43, sd = TRUE, rel = F) # uses a correction factor
MeanMatrixStatistics(VCV_mats2$EF43)



#this says I have singular matricies. Is this because I don't have 30 avalible dimentions?
MMS_summary <- lapply(VCV_mats2, MeanMatrixStatistics)

mat_stats <- data.frame(matrix(unlist(MMS_summary),
                               nrow = 9, ncol = 10, byrow = T))


#I want to bootstrap these estimates. There is a function to do this within evolQG so that is really nice of them. False. This is simply a repetability measure within the matrix? They are clearly getting the bootstrapped matrices but not spitting them out. Easier to just write it myself. Probably worthwile writing for the distances? So I can plot CIs. 
PCAsimilarity(VCV_mats2$EF43, VCV_mats2$EF81)


#############

colnames(mat_stats) <- c("MeanSquaredCorrelation", "PC1", "ICV", 
                         "EigenSd", "Respondability", "Evolvability", 
                         "Conditional.Evolvability", "autonomy",
                         "flexibility", "constraints")

#this matches. 

#now to add poulation 
names(VCV_mats2)

#check this?
mat_stats$population <- c( "LA", "LA", "LA", "cross", "cross", "cross", "HA", "HA", "HA")
mat_stats$line <- c("zi192", "zi251", "zi418", "Zi192Ef43", "Zi192Ef81", "Zi192Ef96", "ef43", "ef81", "ef96")



#Well. This is fun. 
#Used the beeswarm because it was in ID code. Probably don't need it with three points. 
#There is a lot of variation in the high and low and not in the cross, like I saw for allometry. 

#Using ICV here because there are way more individuals in the cross than the parents (so this is scaled by the mean eigenvalue)
ggplot(mat_stats, aes(y =  ICV, x = population )) +
  geom_jitter() +
  labs(x = "", y = "Phenotypic Integration") 

ggplot(mat_stats, aes(y =  PC1, x = population )) +
  geom_jitter() +
  labs(x = "", y = "Eccentricity")


#Going to fit a high vs low model. Not a ton of data but not a ton of predictors. 
PI.mod <- lm(ICV ~ population, data = mat_stats)

summary(PI.mod)

#too much variance in here to see anything
anova(PI.mod)

#just because. It won't show anything. 
ecc.mod <- lm(PC1 ~ population, data = mat_stats)

summary(ecc.mod)

#too much variance in here to see anything
anova(ecc.mod)

#options? This says I can use BootstrapStat to get CIs. Going to need to wrap it by population but that is ok

#this is going to be faster and easier than trying to get the boot functions to work. 
crap <- filter(all.male, line == "ZI251")#first filter data by pop 
crap.resid <- lm(as.matrix(crap[,1:30]) ~ crap$logCS_c)$residuals#model fit resuduals. 
crap.cov <- cov(crap.resid)

CalcEigenVar(crap.cov, sd = TRUE, rel = T)

MeanMatrixStatistics(crap.cov)[2]
MeanMatrixStatistics(crap.cov)[3]
sum(diag(crap.cov))

#I am a dummy and should also probably compute the size of the matirx with this. BUT I do have the problem of not wanting to do start this again. so we will just let it be for now. 
MatStats <- function(data, x) {
  crap <- filter(data, line == x) 
  crap.resid <- lm(as.matrix(crap[,1:30]) ~ crap$logCS_c)$residuals
  crap.sample <- crap.resid[sample(x = nrow(crap.resid), 
                                   size = nrow(crap.resid), replace = T), ]
  crap.cov <- cov(crap.sample)
  ecc <- MeanMatrixStatistics(crap.cov)[2]
  ICV <- MeanMatrixStatistics(crap.cov)[3]
  mat.size <- sum(diag(crap.cov))
  
  return(c(ecc, ICV, mat.size))
} 

#The function works
#Suppressing warnings for stats I don't care about. 
testing <- suppressWarnings(MatStats(all.male, "EF43"))

#working like I think it should. 
suppressWarnings(MatStats(all.male, "EF43"))
suppressWarnings(MatStats(all.male, "EF43"))
suppressWarnings(MatStats(all.male, "EF43"))


#I want to turn this because I like it better in my head. 
mat.stats.ef43 <- t(replicate(1000, suppressWarnings(MatStats(all.male, "EF43"))))
mat.stats.ef81 <- t(replicate(1000, suppressWarnings(MatStats(all.male, "EF81"))))
mat.stats.ef96 <- t(replicate(1000, suppressWarnings(MatStats(all.male, "EF96"))))
mat.stats.zi192 <- t(replicate(1000, suppressWarnings(MatStats(all.male, "ZI192"))))
mat.stats.zi418 <- t(replicate(1000, suppressWarnings(MatStats(all.male, "ZI418"))))

mat.stats.zi251 <- t(replicate(1000, suppressWarnings(MatStats(all.male, "ZI251"))))
mat.stats.Zi192Ef43 <- t(replicate(1000, suppressWarnings(MatStats(all.male, "Zi192Ef43"))))
mat.stats.Zi192Ef81 <- t(replicate(1000, suppressWarnings(MatStats(all.male, "Zi192Ef81"))))
mat.stats.Zi192Ef96 <- t(replicate(1000, suppressWarnings(MatStats(all.male, "Zi192Ef96"))))

boot.mat.names <- c("ef43", "ef81", "ef96", "zi192", "zi251", "zi418", "Zi192Ef96", "Zi192Ef81", "Zi192Ef43")

boot.pc1prop <- data.frame(mat.stats.ef43[,1], mat.stats.ef81[,1], mat.stats.ef96[,1], 
                           mat.stats.zi192[,1], mat.stats.zi251[,1], mat.stats.zi418[,1], 
                           mat.stats.Zi192Ef96[,1], mat.stats.Zi192Ef81[,1], 
                           mat.stats.Zi192Ef43[,1])

names(boot.pc1prop) <- boot.mat.names

class.lines <- data.frame(c("zi192", "zi251", "zi418",  "Zi192Ef96", "Zi192Ef81", "Zi192Ef43", "ef43", "ef81", "ef96"), c("LA", "LA", "LA",  "cross", "cross", "cross", "HA", "HA", "HA"))

names(class.lines) <- c("line", "population")

boot.pc1prop.long <- pivot_longer(boot.pc1prop, everything(), names_to = "line", values_to = "PC1.prop")

boot.pc1prop.long <- left_join(boot.pc1prop.long, class.lines, by = "line")


boot.pc1prop.long$line <- factor(boot.pc1prop.long$line , levels= c("zi192", "zi251", "zi418",  "Zi192Ef96", "Zi192Ef81", "Zi192Ef43", "ef43", "ef81", "ef96"))


PC1prop.plot <- ggplot(boot.pc1prop.long, aes(x = line, y = PC1.prop, fill = population)) +
  stat_eye(alpha = 0.6, 
           p_limits = c(0.05, 0.95), 
           show_interval = FALSE,show_point = FALSE) +
  geom_point(data = mat_stats, aes(x = line, y = PC1)) + 
  theme_classic() + 
  ylab("Eccentricity") + 
  xlab("") + 
  theme(text = element_text(size = 20))


boot.ICV<- data.frame(mat.stats.ef43[,2], mat.stats.ef81[,2], mat.stats.ef96[,2], 
                      mat.stats.zi192[,2], mat.stats.zi251[,2], mat.stats.zi418[,2], 
                      mat.stats.Zi192Ef96[,2], mat.stats.Zi192Ef81[,2], 
                      mat.stats.Zi192Ef43[,2])

names(boot.ICV) <- boot.mat.names

boot.ICV.long <- pivot_longer(boot.ICV, everything(), names_to = "line", values_to = "ICV")
boot.ICV.long <- left_join(boot.ICV.long, class.lines, by = "line")


boot.ICV.long$line <- factor(boot.ICV.long$line , levels=c("zi192", "zi251", "zi418",  "Zi192Ef96", "Zi192Ef81", "Zi192Ef43", "ef43", "ef81", "ef96"))

#gives exactly the same answer as the PC1 proportionality. 
ICV.plot <- ggplot(boot.ICV.long, aes(x = line, y = ICV, fill = population)) + 
  stat_eye(alpha = 0.6, 
           p_limits = c(0.05, 0.95), 
           show_interval = FALSE,show_point = FALSE) + 
  geom_point(data = mat_stats, aes(x = line, y = ICV)) + 
  theme_classic() + 
  ylab("Integration") + 
  xlab("") + 
  theme(text = element_text(size = 20))


boot.size <- data.frame(mat.stats.ef43[,3], mat.stats.ef81[,3], mat.stats.ef96[,3], 
                        mat.stats.zi192[,3], mat.stats.zi251[,3], mat.stats.zi418[,3], 
                        mat.stats.Zi192Ef96[,3], mat.stats.Zi192Ef81[,3], 
                        mat.stats.Zi192Ef43[,3])

names(boot.size) <- boot.mat.names

boot.size.long <- pivot_longer(boot.size, everything(), names_to = "line", values_to = "size")

# class.lines <- data.frame(c("zi192", "zi251", "zi418",  "Zi192Ef96", "Zi192Ef81", "Zi192Ef43", "ef43", "ef81", "ef96"), c("LA", "LA", "LA",  "cross", "cross", "cross", "HA", "HA", "HA"))
# 
# names(class.lines) <- c("line", "population")
boot.size.long <- left_join(boot.size.long, class.lines, by = "line")

boot.size.long$line <- factor(boot.size.long$line , levels=c("zi192", "zi251", "zi418",  "Zi192Ef96", "Zi192Ef81", "Zi192Ef43", "ef43", "ef81", "ef96"))


#Obviously not the best measure when sample sizes are really uneven 
matSize.plot <- ggplot(boot.size.long, aes(x = line, y = size, fill = population)) + 
  stat_eye(alpha = 0.6, 
            p_limits = c(0.05, 0.95), 
            show_interval = FALSE,show_point = FALSE) + 
  #geom_point(data = mat_stats, aes(x = line, y = PC1)) + 
  theme_classic() + 
  ylab("Matrix Size") + 
  xlab("") + 
  theme(text = element_text(size = 20))

#I also want total variance but forgot to include above. 

#measures the observed. 

TotVar_all <- data.frame(colSums(sapply(VCV_mats2, function(x) svd(x)$d)))
names(TotVar_all) <- "totalVaricance"

TotVar_all$line <- c("zi192", "zi251", "zi418",  "Zi192Ef96", "Zi192Ef81", "Zi192Ef43", "ef43", "ef81", "ef96")
TotVar_all$population <- c("LA", "LA", "LA", "cross", "cross","cross", "HA","HA","HA")

sum(svd(VCV_mats2[[1]])$d)

TotVarBoot <- function(data, x) {
  crap <- filter(data, line == x) 
  crap.resid <- lm(as.matrix(crap[,1:30]) ~ crap$logCS_c)$residuals
  crap.sample <- crap.resid[sample(x = nrow(crap.resid), 
                                   size = nrow(crap.resid), replace = T), ]
  crap.cov <- cov(crap.sample)
  totVar <- sum(svd(crap.cov)$d)
  
  return(totVar)
} 

ef43.totVar <- replicate(1000, TotVarBoot(all.male, "EF43"))
ef81.totVar <- replicate(1000, TotVarBoot(all.male, "EF81"))
ef96.totVar <- replicate(1000, TotVarBoot(all.male, "EF96"))
zi192.totVar <- replicate(1000, TotVarBoot(all.male, "ZI192"))
zi251.totVar <- replicate(1000, TotVarBoot(all.male, "ZI251"))
zi418.totVar <- replicate(1000, TotVarBoot(all.male, "ZI418"))
Zi192Ef43.totVar <- replicate(1000, TotVarBoot(all.male, "Zi192Ef43"))
Zi192Ef81.totVar <- replicate(1000, TotVarBoot(all.male, "Zi192Ef81"))
Zi192Ef96.totVar <- replicate(1000, TotVarBoot(all.male, "Zi192Ef96"))


boot.totVar <- data.frame(ef43.totVar, ef81.totVar, ef96.totVar, zi192.totVar, 
                          zi251.totVar, zi418.totVar,Zi192Ef43.totVar, 
                          Zi192Ef81.totVar, Zi192Ef96.totVar)

names(boot.totVar) <- boot.mat.names

class.lines <- data.frame(c("zi192", "zi251", "zi418",  "Zi192Ef96", "Zi192Ef81", "Zi192Ef43", "ef43", "ef81", "ef96"), c("LA", "LA", "LA",  "cross", "cross", "cross", "HA", "HA", "HA"))

names(class.lines) <- c("line", "population")

boot.totVar.long <- pivot_longer(boot.totVar, everything(), names_to = "line", values_to = "totalVaricance")

boot.totVar.long <- left_join(boot.totVar.long, class.lines, by = "line")


boot.totVar.long$line <- factor(boot.pc1prop.long$line , levels= c("zi192", "zi251", "zi418",  "Zi192Ef96", "Zi192Ef81", "Zi192Ef43", "ef43", "ef81", "ef96"))


totVar.plot <- ggplot(boot.totVar.long, aes(x = line, y = totalVaricance, fill = population)) +
  stat_eye(alpha = 0.6, 
           p_limits = c(0.05, 0.95), 
           show_interval = FALSE,show_point = FALSE) +  
  geom_point(data = TotVar_all, aes(x = line, y = totalVaricance)) + 
  theme_classic() + 
  ylab("Total Variance") + 
  xlab("") + 
  theme(text = element_text(size = 20))




#I also was to compute the pairwise distances and pairwise comparisons between matricies. I think that then I want to fit a model for the relationship between estimated relatedness (clearly not with an actual measure) and the pairwise statistic?

#Maybe I should have saved the matrices above because then I could and compare them here.I think that more likley, it will be better to just do it here (this might be super slow) 

#I still want to fit CS within line to not have to worry about shit 
MatPairwiseStats <- function(data, x, y) {
  crap <- filter(data, line == x) 
  crap.resid <- lm(as.matrix(crap[,1:30]) ~ crap$logCS_c)$residuals
  crap.sample <- crap.resid[sample(x = nrow(crap.resid), 
                                   size = nrow(crap.resid), replace = T), ]
  crap.cov <- cov(crap.sample)
  
  crap2 <- filter(data, line ==  y) 
  crap2.resid <- lm(as.matrix(crap2[,1:30]) ~ crap2$logCS_c)$residuals
  crap2.sample <- crap2.resid[sample(x = nrow(crap2.resid), 
                                     size = nrow(crap2.resid), replace = T), ]
  crap2.cov <- cov(crap2.sample)
  
  mat.cor <- PCAsimilarity(crap.cov, crap2.cov)
  #I don't know why this isn't working right now. 
  #mat.dist <- MatrixDistance(crap.cov, crap2.cov)
  
  return(mat.cor)
} 

#seems to work 
MatPairwiseStats(all.male, "EF43", "EF81")
test <- MatPairwiseStats(all.male, "EF43", "EF81")

class(test)

#this is suprisingly quick. 
ef43.ef81MatCor <- replicate(1000, MatPairwiseStats(all.male, "EF43", "EF81"))
ef43.ef96MatCor <- replicate(1000, MatPairwiseStats(all.male, "EF43", "EF96"))
ef43.zi192MatCor <- replicate(1000, MatPairwiseStats(all.male, "EF43", "ZI192"))
ef43.zi418MatCor <- replicate(1000, MatPairwiseStats(all.male, "EF43", "ZI418"))
ef43.zi251MatCor <- replicate(1000, MatPairwiseStats(all.male, "EF43", "ZI251"))
ef43.zi192ef43MatCor <- replicate(1000, MatPairwiseStats(all.male, "EF43", "Zi192Ef43"))
ef43.zi192ef81MatCor <- replicate(1000, MatPairwiseStats(all.male, "EF43", "Zi192Ef81"))
ef43.zi192ef96MatCor <- replicate(1000, MatPairwiseStats(all.male, "EF43", "Zi192Ef96"))

ef81.ef96MatCor <- replicate(1000, MatPairwiseStats(all.male, "EF81", "EF96"))
ef81.zi192MatCor <- replicate(1000, MatPairwiseStats(all.male, "EF81", "ZI192"))
ef81.zi418MatCor <- replicate(1000, MatPairwiseStats(all.male, "EF81", "ZI418"))
ef81.zi251MatCor <- replicate(1000, MatPairwiseStats(all.male, "EF81", "ZI251"))
ef81.zi192ef43MatCor <- replicate(1000, MatPairwiseStats(all.male, "EF81", "Zi192Ef43"))
ef81.zi192ef81MatCor <- replicate(1000, MatPairwiseStats(all.male, "EF81", "Zi192Ef81"))
ef81.zi192ef96MatCor <- replicate(1000, MatPairwiseStats(all.male, "EF81", "Zi192Ef96"))

ef96.zi192MatCor <- replicate(1000, MatPairwiseStats(all.male, "EF96", "ZI192"))
ef96.zi418MatCor <- replicate(1000, MatPairwiseStats(all.male, "EF96", "ZI418"))
ef96.zi251MatCor <- replicate(1000, MatPairwiseStats(all.male, "EF96", "ZI251"))
ef96.zi192ef43MatCor <- replicate(1000, MatPairwiseStats(all.male, "EF96", "Zi192Ef43"))
ef96.zi192ef81MatCor <- replicate(1000, MatPairwiseStats(all.male, "EF96", "Zi192Ef81"))
ef96.zi192ef96MatCor <- replicate(1000, MatPairwiseStats(all.male, "EF96", "Zi192Ef96"))

zi192.zi418MatCor <- replicate(1000, MatPairwiseStats(all.male, "ZI192", "ZI418"))
zi192.zi251MatCor <- replicate(1000, MatPairwiseStats(all.male, "ZI192", "ZI251"))
zi192.zi192ef43MatCor <- replicate(1000, MatPairwiseStats(all.male, "ZI192", "Zi192Ef43"))
zi192.zi192ef81MatCor <- replicate(1000, MatPairwiseStats(all.male, "ZI192", "Zi192Ef81"))
zi192.zi192ef96MatCor <- replicate(1000, MatPairwiseStats(all.male, "ZI192", "Zi192Ef96"))

zi418.zi251MatCor <- replicate(1000, MatPairwiseStats(all.male, "ZI418", "ZI251"))
zi418.zi192ef43MatCor <- replicate(1000, MatPairwiseStats(all.male, "ZI418", "Zi192Ef43"))
zi418.zi192ef81MatCor <- replicate(1000, MatPairwiseStats(all.male, "ZI418", "Zi192Ef81"))
zi418.zi192ef96MatCor <- replicate(1000, MatPairwiseStats(all.male, "ZI418", "Zi192Ef96"))

zi251.zi192ef43MatCor <- replicate(1000, MatPairwiseStats(all.male, "ZI251", "Zi192Ef43"))
zi251.zi192ef81MatCor <- replicate(1000, MatPairwiseStats(all.male, "ZI251", "Zi192Ef81"))
zi251.zi192ef96MatCor <- replicate(1000, MatPairwiseStats(all.male, "ZI251", "Zi192Ef96"))

zi192ef43.zi192ef81MatCor <- replicate(1000, MatPairwiseStats(all.male, "Zi192Ef43", "Zi192Ef81"))
zi192ef43.zi192ef96MatCor <- replicate(1000, MatPairwiseStats(all.male, "Zi192Ef43", "Zi192Ef96"))

zi192ef81.zi192ef96MatCor <- replicate(1000, MatPairwiseStats(all.male, "Zi192Ef81", "Zi192Ef96"))

comparisons.pairwise <- c("ef43.ef81", "ef43.ef96", "ef43.zi192", "ef43.zi418", "ef43.zi251", "ef43.zi192ef43", "ef43.zi192ef81", "ef43.zi192ef96", "ef81.ef96", "ef81.zi192", "ef81.zi418", "ef81.zi251", "ef81.zi192ef43", "ef81.zi192ef81", "ef81.zi192ef96", "ef96.zi192", "ef96.zi418", "ef96.zi251", "ef96.zi192ef43", "ef96.zi192ef81", "ef96.zi192ef96", "zi192.zi418", "zi192.zi251", "zi192.zi192ef43", "zi192.zi192ef81", "zi192.zi192ef96", "zi418.zi251", "zi418.zi192ef43", "zi418.zi192ef81", "zi418.zi192ef96", "zi251.zi192ef43", "zi251.zi192ef81", "zi251.zi192ef96", "zi192ef43.zi192ef81", "zi192ef43.zi192ef96", "zi192ef81.zi192ef96")

#making a list of these vectors so its easy to use sapply
comparisons.pairwise.vecs <- list(ef43.ef81MatCor, ef43.ef96MatCor, ef43.zi192MatCor, ef43.zi418MatCor, ef43.zi251MatCor, ef43.zi192ef43MatCor, ef43.zi192ef81MatCor, ef43.zi192ef96MatCor, ef81.ef96MatCor, ef81.zi192MatCor, ef81.zi418MatCor, ef81.zi251MatCor, ef81.zi192ef43MatCor, ef81.zi192ef81MatCor, ef81.zi192ef96MatCor, ef96.zi192MatCor, ef96.zi418MatCor, ef96.zi251MatCor, ef96.zi192ef43MatCor, ef96.zi192ef81MatCor, ef96.zi192ef96MatCor, zi192.zi418MatCor, zi192.zi251MatCor, zi192.zi192ef43MatCor, zi192.zi192ef81MatCor, zi192.zi192ef96MatCor, zi418.zi251MatCor, zi418.zi192ef43MatCor, zi418.zi192ef81MatCor, zi418.zi192ef96MatCor, zi251.zi192ef43MatCor, zi251.zi192ef81MatCor, zi251.zi192ef96MatCor, zi192ef43.zi192ef81MatCor, zi192ef43.zi192ef96MatCor, zi192ef81.zi192ef96MatCor)

pairwiseMat.cor <- matrix(NA, ncol = 5, nrow = length(comparisons.pairwise))

pairwiseMat.cor[,1] <- comparisons.pairwise

#Now I want to get the CI and mean for each of these for plotting. 
quantile(comparisons.pairwise.vecs[[1]], probs = c(0.05, 0.95))

#Dont know about these quotes. 
pairwiseMat.cor[,2] <- sapply(comparisons.pairwise.vecs, mean)
pairwiseMat.cor[,3] <- sapply(comparisons.pairwise.vecs, 
                              function(x) quantile(x, probs = c(0.05, 0.95))[1])
pairwiseMat.cor[,4] <- sapply(comparisons.pairwise.vecs, 
                              function(x) quantile(x, probs = c(0.05, 0.95))[2])

pairwiseMat.cor <- data.frame(pairwiseMat.cor)
pairwiseMat.cor[2:4] <- lapply(pairwiseMat.cor[2:4], as.numeric)

#Now I want to make some measure of how similar the two lines are genetically 
#This is a first aproximation because I don't actually have all the genomes I would need for this. Ideally, I would make a G matrix and measure the distances. 
#coded as: 0.1 = shared population, 1 = shared parent, 0 = no relationship. This overestimates the shared genetics within populations for sure. 

pairwiseMat.cor[,1]

pairwiseMat.cor[,5] <- c(0.2, 0.2, 0.2, 0, 0, 0.5, 0.1, 0.1, 0.2, 0, 0, 0, 0.1, 0.5, 0.1,0, 0, 0, 0.1, 0.1, 0.5, 0.2, 0.2, 0.5, 0.5, 0.5, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.6, 0.6, 0.6)

names(pairwiseMat.cor) <- c("comparison", "mean", "lower.lim", "upper.lim", "relatedness")

#This plot gets a big thumbs down from me. BUT the crosses to eachother are more similar to eachother and their parents than any of the parents are to eachother. Maybe group into three cats? Or do the model based on: crossxcross, crossxparent, crossxnonparent, HAxHA, LAxLA, HAxLA?

#Also need a much better way to plot. the answer might just be with the comparisons and the model and plotting those estimates because this is a mess. 
ggplot(pairwiseMat.cor, aes(x = relatedness, y = mean)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = lower.lim, ymax = upper.lim)) 


#Now I want to do this with the observations. I calculated this like 500 lines of code aga

VCV_mats2

#first argument: 	
#Single covariance matrix or list of covariance matrices. If cov.x is a single matrix, it is compared to cov.y. If cov.x is a list and no cov.y is supplied, all matrices are compared to each other. If cov.x is a list and cov.y is supplied, all matrices in cov.x are compared to cov.y.

#this doesn't work. 
PCAsimilarity(cov.x = VCV_mats2[1], cov.y = VCV_mats2[2])

#this works. So it all should work? this looks exactly like the set up in the example. 
PCAsimilarity(cov.x = VCV_mats2[[1]], cov.y = VCV_mats2[[9]])

#this does. 
PCAsimilarity(cov.x = VCV_mats2$EF43, cov.y = VCV_mats2$EF81)

#this is a list of cov matricies 
#if I do this manually, they all work? Its not how they are named because I checked that. 
PCAsimilarity(cov.x = VCV_mats2)

#Do want to get the vector of repetability for correction. Maybe that will help?
rep.vec <- rep(NA, 9)

for (i in 1:9) {
  dat <- filter(all.male, line == names(VCV_mats2)[i])
  resid <- lm(as.matrix(dat[,1:30]) ~ dat$logCS_c)$residuals
  rep.val <- BootstrapRep(resid, PCAsimilarity, correlation = FALSE)
  rep.vec[i] <- rep.val
}

PCAsimilarity(cov.x = VCV_mats2, ret.dim = 26)

#this is always right and what I expect. 
dim(eigen(VCV_mats2[[9]])$vectors)

#Adding this measure doesn't change anything at all (probably because my values are so high?)
PCAsimilarity(cov.x = VCV_mats2[[1]], cov.y = VCV_mats2[[9]],  ret.dim = 26, repeat.vector = rep.vec[c(1,9)])

#what if my list is shorter?
#this works 
test <- list(VCV_mats2[[1]], VCV_mats2[[2]], VCV_mats2[[3]])
PCAsimilarity(test)

#I guess my list was not listy enough? fuck this... 
whyDoesThisWork <- list(VCV_mats2[[1]], VCV_mats2[[2]], VCV_mats2[[3]], VCV_mats2[[4]], VCV_mats2[[5]], VCV_mats2[[6]], VCV_mats2[[7]], VCV_mats2[[8]], VCV_mats2[[9]])

#uncorrected below the diag and corrected above. 
vec.corr.toSave <- PCAsimilarity(whyDoesThisWork, ret.dim = 26, repeat.vector = rep.vec)

write.csv(vec.corr.toSave, "../Data/MatrixCorr_PCAsim.csv", quote = FALSE)

#what the fuck. That is a list? Maybe its something else pretending to be a list? Who the fuck knows. 
identical(whyDoesThisWork, VCV_mats2)
setdiff(whyDoesThisWork, VCV_mats2)


#The plan is to measure pairwise corralations between matricies and then make an ordered factor for the type of comparison and fit that as an effect. This is a very imperfect system because it has to be overestimating relatedness. but I at least don't have to do this weird number thing. 


mat.corr.comparisons <- read.csv("../Data/ComparisonTable.csv")
mat.corr.comparisons$comparison <- paste(mat.corr.comparisons$Gentotype1, 
                                         mat.corr.comparisons$Genotype2, sep = ".")

with(mat.corr.comparisons, table(RelatedFactor))

str(mat.corr.comparisons)
#I want to also order my 'relateness' factor before I fit my model.
#Also decided this was a better name. 
mat.corr.comparisons$RelatedFactor <- sub("CrossToOther", "CrossToNonparent", mat.corr.comparisons$RelatedFactor)

mat.corr.comparisons$RelatedFactor <- factor(mat.corr.comparisons$RelatedFactor, 
                                             levels = c("BetweenPop", "WithinPop",
                                                        "CrossToNonparent",
                                                        "CrossToParent",
                                                        "CrossToCross" ))

#Now a quick model. 
mat.corr.mod <- lm(MatCor.Corrected ~ RelatedFactor, 
                   data = mat.corr.comparisons)


#this makes a lot of sense with what I saw before. Little change within or between pop in terms of how similar the matrix is but crosses to the parents and crosses to each other have a much more similar variance
summary(mat.corr.mod)

anova(mat.corr.mod)

plot(emmeans(mat.corr.mod, ~ RelatedFactor))
emmeans(mat.corr.mod, ~ RelatedFactor)

mat.corr.plotshit <- data.frame(emmeans(mat.corr.mod, ~ RelatedFactor))

matrix.corr.plot <- ggplot(mat.corr.plotshit, aes(x = RelatedFactor, y = emmean)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL)) + 
  theme_classic() + 
  ylab("Covariance Matrix Corralation") + 
  xlab("") + 
  ylim(c(0.5, 1.1)) + 
  theme(text = element_text(size = 20))


#########################################)
#I also want to do the projections and then model the shape score. 
all.male.zi192.2 <-data.frame(all.male.zi192, 
                              lm(as.matrix(all.male.zi192[,1:30])~ logCS, 
                                 data = all.male.zi192)$resid)

all.male.zi192.2$shapeScore <- projFunction(as.matrix(all.male.zi192.2[,37:66]), 
                                            male.diff.noAllo)



#I'm still including size here because I want to also ask about that relationship and it matters. 
#but after looking at the anova, maybe thats the wrong choice? and should just have the interaction term?
shapeScoreMod <- lm(shapeScore ~ logCS*line, data = all.male.zi192.2)

summary(shapeScoreMod)
#why is that one?
anova(shapeScoreMod)

#so even this is different 
plot(emmeans(shapeScoreMod, ~ line))

#but this is the more interesting thing to look at 
#this will be better with more. But I think that the conclusion here is that the relationship is really different between crosses. 

#png("../Figures/AfricanShapeScoreEffects.png")
plot(allEffects(shapeScoreMod))
#dev.off()

crap <- allEffects(shapeScoreMod)

################################
###################################################################
#adding the parents in. 

#first wihtout size included. 
pca.forAll <- prcomp(all.male[,1:30])
summary(pca.forAll)

All.pca <- data.frame(all.male, pca.forAll$x)

ggplot(All.pca, aes(x = PC1, y = PC2, col = line)) + 
  geom_point(alpha = 0.4)

#unexpected  
ggplot(All.pca, aes(x = CS, y = PC1, col = line)) + 
  geom_point(alpha = 0.4)

png("../Figures/zi192AllPCA12.png")
ggplot(All.pca, aes(x = PC1, y = PC2, col = line)) + 
  geom_point(alpha = 0.1) + 
  theme_classic()
dev.off()

#Now I want to take out size. I am going to model it here because this at least seems to work. 

All.male2 <- data.frame(all.male, lm(as.matrix(all.male[,1:30]) 
                                     ~ logCS, data = all.male)$residuals)

pca.All.allo <- prcomp(All.male2[,37:66])
summary(pca.All.allo)

All.pca.noAllo <- data.frame(All.male2, pca.All.allo$x)

#a lot of the diffrences were just size related. 
ggplot(All.pca.noAllo, aes(x = PC1, y = PC2, col = line)) + 
  geom_point(alpha = 0.3)

#fine enough.  
ggplot(All.pca.noAllo, aes(x = logCS, y = PC1, col = line)) + 
  geom_point(alpha = 0.3)

# png("../Figures/comparingShapes_All_PCA.png")
# ggplot(All.pca.noAllo, aes(x = PC1, y = PC2, col = line)) + 
#   geom_point(alpha = 0.3) + 
#   theme_classic()
# dev.off()


#########################################
#I also want to do the projections and then model the shape score. 

All.male2$shapeScore <- projFunction(as.matrix(All.male2[,37:66]), 
                                     male.diff.noAllo)


#I'm still including size here because I want to also ask about that relationship and it matters. 
#but after looking at the anova, maybe thats the wrong choice? and should just have the interaction term?
All.male2$population <- relevel(as.factor(All.male2$population), ref = "LA")

#random effect of line because I am looking for population level relationships. Maybe not a great idea since there are so many problems within line
shapeScoreMod <- lmer(shapeScore ~ logCS_c*population + (1|population:line), data = All.male2)

summary(shapeScoreMod)

Anova(shapeScoreMod)

#Well this is encouraging.  
plot(emmeans(shapeScoreMod, ~ population))

population.shape.est <- data.frame(emmeans(shapeScoreMod, ~ population))


shape.score.plot <- ggplot(population.shape.est, aes(x = population, y = emmean)) + 
  geom_jitter(data = All.male2, aes(x = population, y = shapeScore), alpha = 0.1, col = "red") + 
  geom_point() + 
  geom_errorbar(aes(ymin =  asymp.LCL, ymax = asymp.UCL)) + 
  theme_classic() + 
  ylab("Shape Score") + 
  xlab("") + 
  theme(text = element_text(size = 20))

#Not all that interesting to me. 
plot(allEffects(shapeScoreMod))


shit <- data.frame(allEffects(shapeScoreMod))


#Making the second figure. 
#To choose from. 
shape.score.plot
matrix.corr.plot
totVar.plot
ICV.plot
matSize.plot
PC1prop.plot


#moving the legend inside the total variance plot and remove it from the eccentricity plot. 
ICV.plot2 <- ICV.plot + theme(legend.position = "none")
totVar.plot2 <- totVar.plot + theme(legend.position = c(.95, .95), legend.justification = c("right", "top")) + guides(fill=guide_legend(title="Genotype"))


#putting it all together


png("../Figures/shapeMatrixComparisonMainFigs.png", width =6000, height = 3000, units = "px",res = 300)
plot_grid(ICV.plot2, totVar.plot2,  matrix.corr.plot, shape.score.plot, nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))
dev.off()


png("../Figures/eccPlotMatrixCompSupp.png", width =4000, height =2000, units = "px",res = 300)
PC1prop.plot + guides(fill=guide_legend(title="Genotype"))
dev.off()







