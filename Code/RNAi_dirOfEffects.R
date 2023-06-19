#looking at the direction of effects with the RNAi experiment and the selection vector for the genes I liked in the 3R region. 

library(geomorph)
library(tidyverse)
library(abind) #easily bind matrix together in an array. 
library(cowplot)
library(magick)


source('~/Dropbox/DworkinLabSharedMaterial/scripts/WRP_FUNCTIONS.R', chdir = TRUE)

source('~/Dropbox/DworkinLabSharedMaterial/scripts/WINGPLOTSOURCE.R', chdir = TRUE)

wing15.links <-c(1, 7, 
                 7, 12, 
                 12, 13, 
                 13, 14,
                 14, 15,
                 15, 11, 
                 11, 10, 
                 10, 9,
                 9, 8, 
                 8, 6, 
                 6, 12, 
                 6, 2, 
                 8, 13, 
                 10, 14, 
                 9, 3, 
                 4, 5,
                 5, 3, 
                 5, 11)

wing15.links <- matrix(wing15.links, ncol = 2, byrow = TRUE)


dat <- read.delim("../Data/rnaiAllGenes.txt", header = F)
names(dat) <- c("ID", "x1", "y1", "x2", "y2", "x3", "y3","x4", "y4", "x5", "y5","x6", "y6", "x7", "y7","x8", "y8", "x9", "y9", "x10", "y10","x11", "y11","x12", "y12","x13", "y13", "x14", "y14", "x15", "y15")

clean.dat <- (dat %>%
                separate(ID, into = c("perp", "driver", "rnai", "rep", "sex", "fly")))

clean.dat$sex <- tolower(clean.dat$sex)

clean.dat$cross <- paste(clean.dat$driver, clean.dat$rnai, sep = "_")

str(clean.dat)


######I want to superimpose this in a common shape space with other stuff I have, to do this I have a mean shape that I have printed and decided that I will allign to. 
#this comes from some other african stuff I did with the tempature version of this experiment

common.mean <- read.csv("../Data/afMap_commonMean.csv")


#So I need to put my data into geomorph and then add this common mean as the in the array so that everything will allign in that shape space. 
cord <- as.matrix(clean.dat[,7:36])
shape <- arrayspecs(cord, 15, 2)

shape2 <- abind(common.mean, shape, along=3)

gdf <- geomorph.data.frame(landmarks = shape2,
                           cross = c(NA, clean.dat$cross),
                           rep = c(NA, clean.dat$rep),
                           sex = c(NA, clean.dat$sex))


#0 breaks this. needs to be 1 iter. This should constrain it to be in a common shape space with that mean. 
rnai.shape <- gpagen(gdf$landmarks, ProcD = T, max.iter = 1)

#the small one is the common mean. Is that ok? Double check with Ian about the scaling here. 
hist(rnai.shape$Csize)


#those are wings. 
plot(rnai.shape, links = wing15.links)

#This looks pretty normal to me. The top ones are all loco RNAi (which have a diffrent shape so?)
plotOutliers(rnai.shape$coords)

#Now I want to do some PCA stuff but also want to do this outside of geomorph because I hate base plotting. 

wings2 <- as.data.frame(two.d.array(rnai.shape$coords))
#the first line is the common mean that I don't actually care about 
wings2 <- wings2[2:nrow(wings2),]
wings2$CS <- rnai.shape$Csize[-1]
wings2$cross <- as.factor(clean.dat$cross)
#I want to set tripC as the default. 
wings2$cross <- relevel(wings2$cross, ref = "nb_TripC")
wings2$rep <- as.factor(clean.dat$rep)
wings2$sex <- as.factor(clean.dat$sex)

#numbers 
#min 12, max 20 (for vial within sex)
with(wings2, table(cross, rep, sex))

#some more general data crap 
wings2$logCS <- log2(wings2$CS)
wings2$logCS_c <- wings2$logCS - mean(wings2$logCS)

############PCA to look at stuff and make some nice plots#############

pca <- prcomp(wings2[,1:30])

#PC1 has a lot of variance. But I have seen this before (and it is 15 pt shape)
summary(pca)

wings.pca <- data.frame(wings2, pca$x)

#PC1 is really loco. Sex looks like PC2 
#png("../Figures/RNAiShape_PCA12.png", width = 1500, height = 1200, units = "px", res = 300)
pca.plot.12 <- ggplot(wings.pca, aes(x = PC1, y = PC2, col = cross, shape = sex)) + 
  geom_point(alpha = 0.3) + 
  theme_classic() +
  theme(text = element_text(size=12),
        axis.text.x= element_text(size=10), 
        axis.text.y= element_text(size=10), 
        legend.position = "none")

#dev.off()


pca.plot.34 <- ggplot(wings.pca, aes(x = PC3, y = PC4, col = cross, shape = sex)) + 
  geom_point(alpha = 0.3) + 
  theme_classic() +
  theme(text = element_text(size=12),
        axis.text.x= element_text(size=10), 
        axis.text.y= element_text(size=10))

  
png("../Figures/RNAiShape_PCA.png", width = 3000, height = 1200, units = "px", res = 300)
plot_grid(pca.plot.12, pca.plot.34, rel_widths = c(0.8, 1))
dev.off()

#just easier to see. 
ggplot(wings.pca, aes(x = PC1, y = PC2, col = sex)) + 
  geom_point(alpha = 0.3) + 
  theme_classic() 


#not much here. but PC1 and 2 have a lot of variance. 
ggplot(wings.pca, aes(x = PC3, y = PC4, col = cross)) + 
  geom_point(alpha = 0.3)

#PC3 actually really looks like sex. 
ggplot(wings.pca, aes(x = PC3, y = PC4, col = sex)) + 
  geom_point(alpha = 0.3)

#Is PC1 just size? Little bit there but its not explaining the loco and others differences. 

#png("../Figures/RNAiShape_CScor.png", width = 1500, height = 1200, units = "px", res = 300)
CSpcaPlpt <- ggplot(wings.pca, aes(x = CS, y = PC1, col = cross, shape = sex)) + 
  geom_point(alpha = 0.3) + 
  theme_classic() +
  theme(text = element_text(size=12),
        axis.text.x= element_text(size=10), 
        axis.text.y= element_text(size=10))

#dev.off()

#adding the two plots together because they sort of tell a story. 

# png("../Figures/RNAiShape_pcaAndSize.png", width = 3000, height = 1200, units = "px", res = 300)
# plot_grid(pca.plot, CSpcaPlpt, labels = c("A", "B"), rel_widths = c(0.9, 1))
# 
# dev.off()


#Back to geomorph for some models. 
cord <- as.matrix(wings2[,1:30])
shape <- arrayspecs(cord, 15, 2)

wings <- geomorph.data.frame(shape = shape,
                             CS = wings2$CS, 
                             cross = wings2$cross,
                             rep = wings2$rep,
                             sex = wings2$sex, 
                             logCS = wings2$logCS,
                             logCS_c = wings2$logCS_c)

#replicate vial effect is the final term  
shape.mod <-  procD.lm(shape ~ 1 + logCS_c*sex*cross + cross:sex:rep , 
                  SS.type = "II", data = wings)

summary(shape.mod)

anova(shape.mod)

#I don't actually know how to write this better but I am sure there is a way to just drop the term I care about. 
shape.null <-procD.lm(shape ~ 1 + logCS_c + sex + logCS_c:cross + sex:cross + cross:sex:logCS_c + cross:sex:rep , data = wings)

#defult is PD for this function.   
PD <- pairwise(shape.mod, 
               fit.null = shape.null, 
               groups = wings$cross)

#Did I do these distances right? This seems to make sence. 
summary.pairwise(PD)
#interesting contrasts:
#nb_TripC:nb_btn   0.014312379 0.01470722  1.5322852  0.063
#nb_TripC:nb_ef6a  0.022516797 0.01994418  2.1759296  0.011
#nb_TripC:nb_loco  0.039512836 0.03930319  1.6930571  0.047
#nb_TripC:nb_wake  0.012129174 0.01502198  0.6435723  0.252

#this one has sc in the background. 
#nb_loco:w_loco    0.059467559 0.04837985  3.7707674  0.001
#nb_TripC:w_loco   0.036771224 0.02683058  3.7954023  0.001

#this is a homebrew with a totally unknown background. 
#nb_wge:w_wge      0.021554636 0.01753765  3.0894209  0.002
#nb_TripC:nb_wge   0.012390613 0.01449245  0.9802418  0.160

#controls:
#w_loco:w_wge      0.013561398 0.01446907  1.2644659  0.098
#nb_TripC:w_loco   0.036771224 0.02683058  3.7954023  0.001
#nb_TripC:w_wge    0.028882254 0.02112015  3.3409167  0.002

#all of these are super high and I don't really care about this. 
summary.pairwise(PD, test.type = "VC")

##############are the differences just all allometry and boring?########

allo.mod <- procD.lm(shape ~ 1 + cross*logCS_c*sex + cross:sex:rep, 
                     SS.type = "II", data = wings)

allo.mod.null <- procD.lm(shape ~ 1 + cross*sex + logCS_c*sex + logCS_c*cross + cross:sex:rep, 
                          SS.type = "II", data = wings)


allo.pair <- pairwise(allo.mod, 
                      fit.null = allo.mod.null, 
                      groups = wings$cross)

#I mean a couple of these are slightly diffrent. Colouring these would be helpful but fuck base R. Stats are below and make me care less about this.
plotAllometry(allo.mod, wings$CS, wings$cross)

#everything has a really high cor. This means that allometry is not just driving this (as I showed with the PCA above)
summary.pairwise(allo.pair, test.type = "VC")

#############Now I want to get the effect vectors for these###########

#Because I am going to compare to males and females separate, I want to also fit the model separate

split.wings <- coords.subset(wings$shape, wings$sex)

names(split.wings)

#Now I have made this, I am going to have to find a way to make this model. I think the easiest way to do this is to make this a 

male.wings <- geomorph.data.frame(shape = split.wings[["m"]],
                             cross = wings2[wings2$sex == "m",]$cross,
                             rep = wings2[wings2$sex == "m",]$rep,
                             logCS_c = wings2[wings2$sex == "m",]$logCS_c)

male.mod <- procD.lm(shape ~ 1 + logCS_c*cross + cross:rep, 
                     data = male.wings)

summary(male.mod)
anova(male.mod)

#the model fits. But not really what I want? 
male.mod[["coefficients"]][, 1:4]

male.mod.null <- procD.lm(shape ~ 1 + logCS_c + logCS_c:cross + cross:rep, 
                     data = male.wings)

male.pair <-  pairwise(fit = male.mod, fit.null = male.mod.null,
                       groups = male.wings$cross)


#How to get the vector estimates 
#This also does change the answers to make more sense if anything. 
summary(male.pair, show.vectors = TRUE)

#grabbing the estimates for coralations. 
crap.m <- data.frame(male.pair[["LS.means"]][[1]])

tripC.m <- crap.m[1,]
btn.m <- crap.m[2,]
ef6a.m <- crap.m[3,]
loco.m <- crap.m[4,]
takl2.m <- crap.m[5,]
wake.m <- crap.m[6,]
wge.m <- crap.m[7,]
w_loco.m <- crap.m[8,]
w_wge.m <- crap.m[9,]

loco_corrected.m <- loco.m - w_loco.m
wge_corrected.m <- wge.m - w_wge.m

btn_corrected.m <- btn.m - tripC.m
ef6a_corrected.m <- ef6a.m - tripC.m
takl2_corrected.m <- takl2.m - tripC.m
wake_corrected.m <- wake.m - tripC.m

#quickly, looking at the magnitude. Take away here is that loco is bigger than the rest but I knew that. 
sqrt(sum((btn.m - tripC.m)^2)) #0.01262751
sqrt(sum((ef6a.m - tripC.m)^2)) #0.01579819
sqrt(sum((loco.m - tripC.m)^2)) #0.04694306
sqrt(sum((wake.m - tripC.m)^2)) #0.008896102
sqrt(sum((takl2.m- tripC.m)^2)) #0.01756147
sqrt(sum((wge.m - tripC.m)^2)) #0.009724824

sqrt(sum(wge_corrected.m^2)) #0.02017236
sqrt(sum(loco_corrected.m^2)) #0.06168519

#and the females. 
female.wings <- geomorph.data.frame(shape = split.wings[["f"]],
                                  cross = wings2[wings2$sex == "f",]$cross,
                                  rep = wings2[wings2$sex == "f",]$rep,
                                  logCS_c = wings2[wings2$sex == "f",]$logCS_c)

female.mod <- procD.lm(shape ~ 1 + logCS_c*cross + cross:rep, 
                     data = female.wings)

summary(female.mod)
anova(female.mod)

#the model fits. But not really what I want? 
female.mod[["coefficients"]][, 1:4]

female.mod.null <- procD.lm(shape ~ 1 + logCS_c + logCS_c:cross + cross:rep, 
                          data = female.wings)

female.pair <-  pairwise(fit = female.mod, fit.null = female.mod.null,
                       groups = female.wings$cross)


#How to get the vector estimates 
#Why do the females always look so shitty. 
summary(female.pair, show.vectors = TRUE)

#grabbing the estimates for coralations. 
crap.f <- data.frame(female.pair[["LS.means"]][[1]])

tripC.f <- crap.f[1,]
btn.f <- crap.f[2,]
ef6a.f <- crap.f[3,]
loco.f <- crap.f[4,]
takl2.f <- crap.f[5,]
wake.f <- crap.f[6,]
wge.f <- crap.f[7,]
w_loco.f <- crap.f[8,]
w_wge.f <- crap.f[9,]

loco_corrected.f <- loco.f - w_loco.f
wge_corrected.f <- wge.f - w_wge.f

btn_corrected.f <- btn.f - tripC.f
ef6a_corrected.f <- ef6a.f - tripC.f
takl2_corrected.f <- takl2.f - tripC.f
wake_corrected.f <- wake.f - tripC.f

#quickly, looking at the magnitude. same as before. 
sqrt(sum(btn_corrected.f^2)) #0.01348776
sqrt(sum(ef6a_corrected.f^2)) #0.01348776
sqrt(sum(loco_corrected.f^2)) #0.06383939
sqrt(sum(wake_corrected.f^2)) #0.01036505
sqrt(sum(takl2_corrected.f^2)) #0.02118132
sqrt(sum(wge_corrected.f^2)) #0.02120699


#Now I want to compare these vectors with the parental africans. 
kp.raw <- read.delim("../Data/KP_purelines_15lmk.txt", header = F)

head(kp.raw)

names(kp.raw) <- c("ID", "x1", "y1", "x2", "y2", "x3", "y3","x4", "y4", "x5", "y5","x6", "y6", "x7", "y7","x8", "y8", "x9", "y9", "x10", "y10","x11", "y11","x12", "y12","x13", "y13", "x14", "y14", "x15", "y15")


clean.kp <- (kp.raw %>%
               separate(ID, into = c("perp", "line", "slide", "sex", "fly")))

#superimposing. 
cord <- as.matrix(clean.kp[,6:35])
shape <- arrayspecs(cord, 15, 2)

shape2 <- abind(common.mean, shape, along=3)

gdf.pure <- geomorph.data.frame(landmarks = shape2,
                           line = c(NA, clean.kp$line),
                           sex = c(NA, clean.kp$sex))

#0 breaks this. needs to be 1 iter. This should constrain it to be in a common shape space with that mean. 
pure.shape <- gpagen(gdf.pure$landmarks, ProcD = T, max.iter = 1)

#same as before. 
hist(pure.shape$Csize)

#those are wings still. I have checked these data so many times. 
plot(pure.shape, links = wing15.links)

#taking it all out to make my vecs. 
pure.wings <- as.data.frame(two.d.array(pure.shape$coords))
#the first line is the common mean that I don't actually care about 
pure.wings <- pure.wings[2:nrow(pure.wings),]
pure.wings$CS <- pure.shape$Csize[-1]
pure.wings$line <- as.factor(clean.kp$line)
pure.wings$population <-as.factor(ifelse(grepl("Z", pure.wings$line), "LA", "HA"))
pure.wings$sex <- as.factor(clean.kp$sex)
# pure.wings$fly <- clean.kp$fly
# pure.wings$slide <- clean.kp$slide


#some more general data crap 
pure.wings$logCS <- log2(pure.wings$CS)
pure.wings$logCS_c <- pure.wings$logCS - mean(pure.wings$logCS)

#writing this out now so I don't have to keep doing this. 
#write.csv(pure.wings, "../Data/pureWings_15lmk_superimposed.csv", row.names = FALSE)

# write.csv(pure.wings, "../Data/pureWings_15lmk_superimposed_withFly.csv", row.names = FALSE)

#males first. 
pure.males <- filter(pure.wings, sex == "M")

pure.males.fit <- data.frame(pure.males, lm(as.matrix(pure.males[,1:30])~ logCS_c, data = pure.males)$resid)

males.ha <- colMeans(pure.males.fit[pure.males.fit$population == "HA", 37:66])
males.la <- colMeans(pure.males.fit[pure.males.fit$population == "LA", 37:66])

males.diff <- males.ha - males.la

#That is unexpected. 
sqrt(sum(males.diff^2)) #0.005161079

#these are something. Some directions shared. not all. Much lower than I got with the other stuff (but I also estimated that diff)
cor(males.diff, t(btn_corrected.m)) #-0.008327475
cor(males.diff, t(ef6a_corrected.m)) #-0.4009203
cor(males.diff, t(loco_corrected.m)) #0.4735495
cor(males.diff, t(takl2_corrected.m)) #-0.1542452
cor(males.diff, t(wge_corrected.m)) #0.4708442
cor(males.diff, t(wake_corrected.m)) #-0.2079081

#Now with females. 

pure.females <- filter(pure.wings, sex == "F")

pure.females.fit <- data.frame(pure.females, lm(as.matrix(pure.females[,1:30])~ logCS_c, data = pure.females)$resid)

females.ha <- colMeans(pure.females.fit[pure.females.fit$population == "HA", 37:66])
females.la <- colMeans(pure.females.fit[pure.females.fit$population == "LA", 37:66])

females.diff <- females.ha - females.la

#That is unexpected. 
sqrt(sum(females.diff^2)) #0.005106257

cor(females.diff, t(btn_corrected.f)) #-0.5718972
cor(females.diff, t(ef6a_corrected.f)) #-0.8263311
cor(females.diff, t(loco_corrected.f)) #0.7315722
cor(females.diff, t(takl2_corrected.f)) #-0.6425486
cor(females.diff, t(wge_corrected.f)) #0.5654785


#Now to compare males and female effects because those seem to be diff. 
cor(t(btn_corrected.m), t(btn_corrected.f)) #0.5908595
cor(t(ef6a_corrected.m), t(ef6a_corrected.f)) #0.6519004
cor(t(loco_corrected.m), t(loco_corrected.f)) #0.9861952
cor(t(takl2_corrected.m), t(takl2_corrected.f)) #0.8393971
cor(t(wge_corrected.m), t(wge_corrected.f)) #0.8965091

cor(females.diff, males.diff) #0.9203572

#so all of these are generally in the same direction ish. but maybe not the particularly interesting directions?
#where to go from here?


#But I do want to make a figure with all of these effects now!!! 
#First, I need to take all the effect vecs and turn them into 2D matricies. I think that I can make them a list and it will just work better that way. 

#I need a vector mean shape for this. I also think it needs to be a geomorph object. 
m <- mshape(pure.shape$coords)

#function to turn vector into 2d format 
vec2matrix <- function(x) {
  y <- unlist(x)
  m - matrix(y, nrow = 15, ncol = 2, byrow = TRUE)
}


eff.mat <- matrix(unlist(loco_corrected.m), nrow = 15, ncol = 2, byrow = TRUE)

crap <- m - eff.mat
poop <- m - matrix(unlist(loco_corrected.m), nrow = 15, ncol = 2, byrow = TRUE)


m.eff.vecs <- list(males.diff, loco_corrected.m, wge_corrected.m, btn_corrected.m, ef6a_corrected.m, takl2_corrected.m, wake_corrected.m)

#names. this was for troubleshooting but now I think this is actually worth keeping because its the order of effect mats in the array. 
m.effs <- c("loco_corrected.m", "wge_corrected.m", "btn_corrected.m", "ef6a_corrected.m", "takl2_corrected.m", "wake_corrected.m")

m.eff.mats <- sapply(m.eff.vecs, vec2matrix, simplify = "array")

#for me for nice plotting 
source("~/Dropbox/KatiePelletier/KP_geomorphwingfunctions.R")

#Why is this so small? None of these capture the posterior cross vein shift. 

png("../Figures/highLowMaleEff_10x.png")
plotRefToTarget(m, m.eff.mats[,,1], 
                method = "points",
                links = wing15.links, mag = 10, 
                gridPars=wing.spec)
dev.off()


#That is a much more narrow wing
png("../Figures/locoMaleEff_1.5x.png")
plotRefToTarget(m, m.eff.mats[,,2], 
                method = "points",
                links = wing15.links, mag =1.5, 
                gridPars=wing.spec)
dev.off()

png("../Figures/wgeMaleEff_3x.png")
plotRefToTarget(m, m.eff.mats[,,3], 
                method = "points",
                links = wing15.links, mag = 3, 
                gridPars=wing.spec)
dev.off()

png("../Figures/btnMaleEff_5x.png")
plotRefToTarget(m, m.eff.mats[,,4], 
                method = "points",
                links = wing15.links, mag = 5, 
                gridPars=wing.spec)
dev.off()

#small p-vein shift. maybe
png("../Figures/ef6aMaleEff_3x.png")
plotRefToTarget(m, m.eff.mats[,,5], 
                method = "points",
                links = wing15.links, mag = 3, 
                gridPars=wing.spec)
dev.off()

png("../Figures/takl2MaleEff_3x.png")
plotRefToTarget(m, m.eff.mats[,,6], 
                method = "points",
                links = wing15.links, mag = 3, 
                gridPars=wing.spec)
dev.off()

#truly no change here
png("../Figures/wakeMaleEff_5x.png")
plotRefToTarget(m, m.eff.mats[,,7], 
                method = "points",
                links = wing15.links, mag = 5, 
                gridPars=wing.spec)
dev.off()

#females. 
eff.mat.f <- matrix(unlist(loco_corrected.f), nrow = 15, ncol = 2, byrow = TRUE)

f.eff.vecs <- list(females.diff, loco_corrected.f, wge_corrected.f, btn_corrected.f, ef6a_corrected.f, takl2_corrected.f, wake_corrected.f)

#names. this was for troubleshooting but now I think this is actually worth keeping because its the order of effect mats in the array. 
f.effs <- c("loco_corrected.f", "wge_corrected.f", "btn_corrected.f", "ef6a_corrected.f", "takl2_corrected.f", "wake_corrected.f")

f.eff.mats <- sapply(f.eff.vecs, vec2matrix, simplify = "array")

#Why is this so small? None of these capture the posterior cross vein shift. 
plotRefToTarget(m, f.eff.mats[,,1], 
                method = "points",
                links = wing15.links, mag = 10, 
                gridPars=wing.spec)


#That is a much more narrow wing
plotRefToTarget(m, f.eff.mats[,,2], 
                method = "points",
                links = wing15.links, mag =1.5, 
                gridPars=wing.spec)

plotRefToTarget(m, f.eff.mats[,,3], 
                method = "points",
                links = wing15.links, mag = 3, 
                gridPars=wing.spec)

plotRefToTarget(m, f.eff.mats[,,4], 
                method = "points",
                links = wing15.links, mag = 3, 
                gridPars=wing.spec)

#small p-vein shift. maybe 
plotRefToTarget(m, f.eff.mats[,,5], 
                method = "points",
                links = wing15.links, mag = , 
                gridPars=wing.spec)

plotRefToTarget(m, f.eff.mats[,,6], 
                method = "points",
                links = wing15.links, mag = 3, 
                gridPars=wing.spec)

#truly no change here
plotRefToTarget(m, f.eff.mats[,,7], 
                method = "points",
                links = wing15.links, mag = 5, 
                gridPars=wing.spec)


#The super small effect size makes me worried. 
#I wonder if this is because of the model. going to look at everything with CS included and then I am also going to make the vector like I actually do for selection (probably a better comparison)

males.ha.lm <- colMeans(pure.males.fit[pure.males.fit$population == "HA", 1:30])
males.la.lm <- colMeans(pure.males.fit[pure.males.fit$population == "LA", 1:30])

males.diff.lm <- males.ha.lm - males.la.lm

#this is much better 
sqrt(sum(males.diff.lm^2)) #0.01929885

#these are something. Some directions shared. not all. Much lower than I got with the other stuff (but I also estimated that diff)
cor(males.diff.lm, t(btn_corrected.m)) #-0.09598382
cor(males.diff.lm, t(ef6a_corrected.m)) #-0.4555911
cor(males.diff.lm, t(loco_corrected.m)) #-0.3234394
cor(males.diff.lm, t(takl2_corrected.m)) #0.2897442
cor(males.diff.lm, t(wge_corrected.m)) #-0.4261016

#and with the PCA method that I use for my other vector. The problem here is that because I drop PC1, I am not TOTALLY sure how to backtransform? maybe take the mean? I need to ask Ian this. 

names(pure.males)

pure.pca.m <- prcomp(pure.males[c(1:30, 35)])

summary(pure.pca.m)

pure.wings.pca.m <- data.frame(pure.males, pure.pca.m$x)

# all good. 
ggplot(pure.wings.pca.m, aes(x = logCS, y = PC1, col = line)) + 
  geom_point(alpha = 0.5)

#This should be moved somewhere else
#But this looks the same as other ones I've looked at. I should really do some more close looks at the high low diffrences. I actually have that from the geomorph class and should just pull that back up and compare distances in the full spline and 15pt shapes. 
ggplot(pure.wings.pca.m, aes(x = PC2, y = PC3, col = line)) + 
  geom_point(alpha = 0.5)


#Now I want to get the diff vector, dropping PC1.Taking the rest avalible dimentions (30 + 1 (size) - 4 (superimposition))
ha.m.pcs <- colMeans(pure.wings.pca.m[pure.wings.pca.m$population == "HA",38:64])
la.m.pcs <- colMeans(pure.wings.pca.m[pure.wings.pca.m$population == "LA",38:64])

male.diff.pcs <- ha.m.pcs - la.m.pcs

#about the same. I should look back at older code for this. 
#not actually worth doing the rest until I look at the full spline vs this 
sqrt(sum(male.diff.pcs^2)) #0.005152562

#making a nice RNAi figure 
library(cowplot)
library(magick)

highlow.raw <- image_read("../Figures/highLowMaleEff_10x.png")
highlow.crop <- image_trim(highlow.raw)
highlow <- ggdraw() + draw_image(highlow.crop)

loco.raw <- image_read("../Figures/locoMaleEff_1.5x.png")
loco.crop <- image_trim(loco.raw)
loco <- ggdraw() + draw_image(loco.crop)

wge.raw <- image_read("../Figures/wgeMaleEff_3x.png")
wge.crop <- image_trim(wge.raw)
wge <- ggdraw() + draw_image(wge.crop)

btn.raw <- image_read("../Figures/btnMaleEff_5x.png")
btn.crop <- image_trim(btn.raw)
btn <- ggdraw() + draw_image(btn.crop)

ef6a.raw <- image_read("../Figures/ef6aMaleEff_3x.png")
ef6a.crop <- image_trim(ef6a.raw)
ef6a <- ggdraw() + draw_image(ef6a.crop)

takl2.raw <- image_read("../Figures/takl2MaleEff_3x.png")
takl2.crop <- image_trim(takl2.raw)
takl2 <- ggdraw() + draw_image(takl2.crop)

wake.raw <- image_read("../Figures/wakeMaleEff_5x.png")
wake.crop <- image_trim(wake.raw)
wake <- ggdraw() + draw_image(wake.crop)

wing.labs <- c("loco", "wge", "btn","ef6a",
               "takl2", "wake")


png("../Figures/RNAiEffectsAndCor.png", width =2000, height = 1000, units = "px",res = 300)
plot_grid(loco, wge, btn, ef6a, takl2, wake, 
          ncol = 3, nrow = 2, 
          labels = wing.labs, label_fontface = "italic")
dev.off()

####################################################
#bootstrap CI for cor
###################################################
#using an old tutorial ID wrote to work from here. 

#Question for ID: because I don't estimate the high low from the model, what does that mean?
#clean dat is the RNAi data. 
rnai.males <- filter(wings2, sex == "m")
#pure.males is all the males 

#So I want to estimate the effect of RNAi knockdown from the model and also a second model that gets the sex effect. 

#the high-low vec, corrected for allometry would be just the population coef here?
#ignoring line level effects for now. Will meet with ID to talk more.
#lm should give same answers

levels(pure.males$population)
pure.males$population <- relevel(pure.males$population, ref = "LA")

#centre logCS 
mlm.mod1 <- manova(as.matrix(pure.males[,1:30]) ~ logCS_c*population, 
                   data = pure.males)

mlm.mod1_alt <- lm(as.matrix(pure.males[,1:30]) ~ logCS_c*population, 
                   data = pure.males)

coef1 <- mlm.mod1$coef[3,]

#also ignoring replicate vial for now to make this easier.

#I want to do this with just one at a time. 
rnai.males.btn <- droplevels(rnai.males[rnai.males$cross == "nb_TripC" | rnai.males$cross == "nb_btn", ])
levels(rnai.males.btn$cross)


mlm.mod2 <- manova(as.matrix(rnai.males.btn[,1:30]) ~ logCS_c*cross + cross:rep, 
                   data = rnai.males.btn)

coef(mlm.mod2)

#So for some of these, it is easy because the base level is the trip control but for some I have to find the difference. Going to start with btn because its the easiest choice as its a vector I can just grab. 
coef2 <- mlm.mod2$coef[3,]

btn.coef <- mlm.mod2$coef[3,]

#this is in the lab dropbox functions. 
#angs <- ang.vec.abs(coef1, coef2)
obs.cor.btn <- cor(coef1, coef2)
#MP code also gets the mag of the vecs, but I only really care about the direction because of course the mag is diffrent and not in an interesting way. 

#Now I want to make this a function. 
#I need to modify this to take arguments for the line to pull out of the second model. 
#this works if I have about equal numbers between the control and the test 
boot_trial_rnaiEff <- function(data1 = data.subset1, data2 = data.subset2, vars = 1:30) {
  data.sample1 <- data1[sample(x = nrow(data1), size = nrow(data1), replace = T), ]
  
  mlm.mod1 <- manova(as.matrix(data.sample1[,vars]) ~ logCS_c*population, 
                     data = data.sample1)
  coef1 <-    mlm.mod1$coef[3,]  
  
  data.sample2 <- data2[sample(x = nrow(data2), size = nrow(data2), replace = T), ]
  
  mlm.mod2 <- manova(as.matrix(data.sample2[,vars]) ~ 1 + logCS_c*cross + cross:rep, 
                     data = data.sample2)
  coef2 <-    mlm.mod2$coef[3,]
  angs <- cor(coef1, coef2)
  
  return(angs)
  #cor(coef1, coef2)
}

#testing. Yeah. that seems to work. 
boot_trial_rnaiEff(pure.males, rnai.males.btn)
boot_trial_rnaiEff(pure.males, rnai.males.btn)
boot_trial_rnaiEff(pure.males, rnai.males.btn)

#now to replicate.
boot_samples_btn_males <- replicate(1000, boot_trial_rnaiEff(pure.males, rnai.males.btn))

#ask Ian what margin argument means. 
#these are giant. 
apply(boot_samples_btn_males, MARGIN = 1, quantile, probs = c(0.025, 0.5, 0.975))

obs.cor.btn
quantile(boot_samples_btn_males, probs = c(0.025, 0.5, 0.975))

apply(boot_samples_btn_males, MARGIN = 1, sd)

sd(boot_samples_btn_males)

apply(boot_samples_btn_males, MARGIN = 1, mean)

#should be similar or else we need to do bias corrections. 
obs.cor.btn
mean(boot_samples_btn_males)

plot(density(boot_samples_btn_males),
     xlab = "Vector Correlation", main = "")
abline(v = obs.cor.btn, col = "red", lwd = 2)


#######
#I want to do this with just one at a time. 
rnai.males.ef6a <- droplevels(rnai.males[rnai.males$cross == "nb_TripC" | rnai.males$cross == "nb_ef6a", ])
levels(rnai.males.ef6a$cross)

boot_samples_ef6a_males <- replicate(1000, boot_trial_rnaiEff(pure.males, rnai.males.ef6a))

ef6a.coef <- manova(as.matrix(rnai.males.ef6a[,1:30]) ~ logCS_c*cross + cross:rep, 
       data = rnai.males.ef6a)$coef[3,]

obs.cor.ef6a <- cor(ef6a.coef, coef1)

quantile(boot_samples_ef6a_males, probs = c(0.025, 0.5, 0.975))
sd(boot_samples_ef6a_males)
mean(boot_samples_ef6a_males)
obs.cor.ef6a


plot(density(boot_samples_ef6a_males),
     xlab = "Vector Correlation", main = "")
abline(v = obs.cor.ef6a, col = "red", lwd = 2)


#######
#I want to do this with just one at a time. 
rnai.males.takl2 <- droplevels(rnai.males[rnai.males$cross == "nb_TripC" | rnai.males$cross == "nb_takl2", ])
levels(rnai.males.takl2$cross)

boot_samples_takl2_males <- replicate(1000, boot_trial_rnaiEff(pure.males, rnai.males.takl2))

takl2.coef <- manova(as.matrix(rnai.males.takl2[,1:30]) ~ logCS_c*cross + cross:rep, 
                    data = rnai.males.takl2)$coef[3,]

obs.cor.takl2 <- cor(takl2.coef, coef1)

quantile(boot_samples_takl2_males, probs = c(0.05, 0.5, 0.95))
sd(boot_samples_takl2_males)
mean(boot_samples_takl2_males)
obs.cor.takl2


plot(density(boot_samples_takl2_males),
     xlab = "Vector Correlation", main = "")
abline(v = obs.cor.takl2, col = "red", lwd = 2)


#######
#I want to do this with just one at a time. 
rnai.males.wake <- droplevels(rnai.males[rnai.males$cross == "nb_TripC" | rnai.males$cross == "nb_wake", ])
levels(rnai.males.wake$cross)

boot_samples_wake_males <- replicate(1000, boot_trial_rnaiEff(pure.males, rnai.males.wake))

wake.coef <- manova(as.matrix(rnai.males.wake[,1:30]) ~ logCS_c*cross + cross:rep, 
                     data = rnai.males.wake)$coef[3,]

obs.cor.wake <- cor(wake.coef, coef1)

quantile(boot_samples_wake_males, probs = c(0.025, 0.5, 0.975))
sd(boot_samples_wake_males)
mean(boot_samples_wake_males)
obs.cor.wake


plot(density(boot_samples_wake_males),
     xlab = "Vector Correlation", main = "")
abline(v = obs.cor.wake, col = "red", lwd = 2)

#######
#I want to do this with just one at a time. 
rnai.males.wge <- droplevels(rnai.males[rnai.males$cross == "nb_wge" | rnai.males$cross == "w_wge", ])
levels(rnai.males.wge$cross)
rnai.males.wge$cross <- relevel(rnai.males.wge$cross, ref = "w_wge")
levels(rnai.males.wge$cross)

boot_samples_wge_males <- replicate(1000, boot_trial_rnaiEff(pure.males, rnai.males.wge))

wge.coef <- manova(as.matrix(rnai.males.wge[,1:30]) ~ logCS_c*cross + cross:rep, 
                    data = rnai.males.wge)$coef[3,]

obs.cor.wge <- cor(wge.coef, coef1)

quantile(boot_samples_wge_males, probs = c(0.05, 0.5, 0.95))
sd(boot_samples_wge_males)
mean(boot_samples_wge_males)
obs.cor.wge


plot(density(boot_samples_wge_males),
     xlab = "Vector Correlation", main = "")
abline(v = obs.cor.wge, col = "red", lwd = 2)

#######
#I want to do this with just one at a time. 
rnai.males.loco <- droplevels(rnai.males[rnai.males$cross == "nb_loco" | rnai.males$cross == "w_loco", ])
levels(rnai.males.loco$cross)
rnai.males.loco$cross <- relevel(rnai.males.loco$cross, ref = "w_loco")
levels(rnai.males.loco$cross)

boot_samples_loco_males <- replicate(1000, boot_trial_rnaiEff(pure.males, rnai.males.loco))

loco.coef <- manova(as.matrix(rnai.males.loco[,1:30]) ~ logCS_c*cross + cross:rep, 
                    data = rnai.males.loco)$coef[3,]

obs.cor.loco <- cor(loco.coef, coef1)

quantile(boot_samples_loco_males, probs = c(0.05, 0.5, 0.95))
sd(boot_samples_loco_males)
mean(boot_samples_loco_males)
obs.cor.loco


#this makes sense becasuse the loco wings are fatter. 
plot(density(boot_samples_loco_males),
     xlab = "Vector Correlation", main = "")
abline(v = obs.cor.loco, col = "red", lwd = 2)


#Now I want to plot these density plots with the most interesting shape changes. 
#Making a df to plot this all. 

gene.vec <- c("btn", "ef6a", "loco", "takl2", "wake", "wge")

obs.cor.all <- c(obs.cor.btn, obs.cor.ef6a, obs.cor.loco, obs.cor.takl2, obs.cor.wake, obs.cor.wge)

obs.cor.all2 <- data.frame(gene.vec, obs.cor.all)
names(obs.cor.all2) <- c("Gene", "observed")


boot.df <- data.frame(boot_samples_btn_males, boot_samples_ef6a_males,
                      boot_samples_loco_males, boot_samples_takl2_males, 
                      boot_samples_wake_males, boot_samples_wge_males)

names(boot.df) <- gene.vec

boot.df.long <- pivot_longer(boot.df, everything(), names_to = "Gene", values_to = "cor")

corr.dist.fig <- ggplot(boot.df.long, aes(x = Gene, y = cor)) + 
  stat_eye(alpha = 0.6, 
           p_limits = c(0.05, 0.95), 
           show_interval = FALSE,show_point = FALSE) + 
  geom_point(data = obs.cor.all2, aes(x = Gene, y = observed), col = "red") +
  geom_hline(yintercept = 0, lwd = 1, lty = 2, alpha = 0.3) + 
  theme_classic() + 
  xlab("") + 
  ylab("corralation with altitude effect vector") +  
  theme(axis.text.x = element_text(face = "italic"), text = element_text(size = 20)) 

alt.pan <- plot_grid(highlow, labels = "Altitude Effect", label_fontface = "italic")

wings.pan <- plot_grid(alt.pan, loco, takl2, wge, ncol = 1, labels = c("", "loco", "takl2", "wge"), label_fontface = "bold.italic", hjust = -2)
  

png("../Figures/RNAi_effectFigWithdist.png", height = 2000, width = 4000, units = "px", res = 300)
plot_grid(corr.dist.fig, wings.pan, nrow = 1)
dev.off()


#I also want to compare the RNAi eff vectors to eachother for the supplement (and to look for shared directions of effects to support the more than one gene argument.)

#need a new boot function 
boot_trial_rnaiTOeachother <- function(data1 = data.subset1, data2 = data.subset2, vars = 1:30) {
  data.sample1 <- data1[sample(x = nrow(data1), size = nrow(data1), replace = T), ]
  
  mlm.mod1 <- manova(as.matrix(data.sample1[,vars]) ~ 1 + logCS_c*cross + cross:rep, 
                     data = data.sample1)
  coef1 <-    mlm.mod1$coef[3,]  
  
  data.sample2 <- data2[sample(x = nrow(data2), size = nrow(data2), replace = T), ]
  
  mlm.mod2 <- manova(as.matrix(data.sample2[,vars]) ~ 1 + logCS_c*cross + cross:rep, 
                     data = data.sample2)
  coef2 <-    mlm.mod2$coef[3,]
  angs <- cor(coef1, coef2)
  
  return(angs)
  #cor(coef1, coef2)
}

#thankfully I have the DF already made to do this all with and can just recalculate with what I need 

boot_samples_locoWge <- replicate(1000, boot_trial_rnaiTOeachother(rnai.males.loco, rnai.males.wge))

obs.cor.locoWge <- cor(loco.coef, wge.coef)

quantile(boot_samples_locoWge, probs = c(0.05, 0.5, 0.95))
sd(boot_samples_locoWge)
mean(boot_samples_locoWge)
obs.cor.locoWge

#################################
boot_samples_locoBtn <- replicate(1000, boot_trial_rnaiTOeachother(rnai.males.loco, rnai.males.btn))

obs.cor.locoBtn <- cor(loco.coef, btn.coef)

quantile(boot_samples_locoBtn, probs = c(0.05, 0.5, 0.95))
sd(boot_samples_locoBtn)
mean(boot_samples_locoBtn)
obs.cor.locoBtn


#################################
boot_samples_wgeBtn <- replicate(1000, boot_trial_rnaiTOeachother(rnai.males.wge, rnai.males.btn))

obs.cor.wgeBtn <- cor(wge.coef, btn.coef)

quantile(boot_samples_wgeBtn, probs = c(0.05, 0.5, 0.95))
sd(boot_samples_wgeBtn)
mean(boot_samples_wgeBtn)
obs.cor.wgeBtn

#################################
boot_samples_ef6aLoco <- replicate(1000, boot_trial_rnaiTOeachother(rnai.males.ef6a, rnai.males.loco))

obs.cor.ef6aLoco <- cor(ef6a.coef, loco.coef)

quantile(boot_samples_ef6aLoco, probs = c(0.05, 0.5, 0.95))
sd(boot_samples_ef6aLoco)
mean(boot_samples_ef6aLoco)
obs.cor.ef6aLoco

#################################
boot_samples_ef6aWge <- replicate(1000, boot_trial_rnaiTOeachother(rnai.males.ef6a, rnai.males.wge))

obs.cor.ef6aWge <- cor(ef6a.coef, wge.coef)

quantile(boot_samples_ef6aWge, probs = c(0.05, 0.5, 0.95))
sd(boot_samples_ef6aWge)
mean(boot_samples_ef6aWge)
obs.cor.ef6aWge

#################################
boot_samples_ef6abtn <- replicate(1000, boot_trial_rnaiTOeachother(rnai.males.ef6a, rnai.males.btn))

obs.cor.ef6abtn  <- cor(ef6a.coef, btn.coef)

quantile(boot_samples_ef6abtn, probs = c(0.05, 0.5, 0.95))
sd(boot_samples_ef6abtn)
mean(boot_samples_ef6abtn)
obs.cor.ef6abtn

#################################
boot_samples_takl2Loco <- replicate(1000, boot_trial_rnaiTOeachother(rnai.males.takl2, rnai.males.loco))

obs.cor.takl2Loco  <- cor(takl2.coef, loco.coef)

quantile(boot_samples_takl2Loco, probs = c(0.05, 0.5, 0.95))
sd(boot_samples_takl2Loco)
mean(boot_samples_takl2Loco)
obs.cor.takl2Loco

#################################
boot_samples_takl2wge <- replicate(1000, boot_trial_rnaiTOeachother(rnai.males.takl2, rnai.males.wge))

obs.cor.takl2wge  <- cor(takl2.coef, wge.coef)

quantile(boot_samples_takl2wge, probs = c(0.05, 0.5, 0.95))
sd(boot_samples_takl2wge)
mean(boot_samples_takl2wge)
obs.cor.takl2wge

#################################
boot_samples_takl2btn <- replicate(1000, boot_trial_rnaiTOeachother(rnai.males.takl2, rnai.males.btn))

obs.cor.takl2btn <- cor(takl2.coef, btn.coef)

quantile(boot_samples_takl2btn, probs = c(0.05, 0.5, 0.95))
sd(boot_samples_takl2btn)
mean(boot_samples_takl2btn)
obs.cor.takl2btn

#################################
boot_samples_takl2ef6a <- replicate(1000, boot_trial_rnaiTOeachother(rnai.males.takl2, rnai.males.ef6a))

obs.cor.takl2ef6a <- cor(takl2.coef, ef6a.coef)

quantile(boot_samples_takl2ef6a, probs = c(0.05, 0.5, 0.95))
sd(boot_samples_takl2ef6a)
mean(boot_samples_takl2ef6a)
obs.cor.takl2ef6a

#################################
boot_samples_wakeLoco <- replicate(1000, boot_trial_rnaiTOeachother(rnai.males.wake, rnai.males.loco))

obs.cor.wakeLoco <- cor(wake.coef, loco.coef)

quantile(boot_samples_wakeLoco, probs = c(0.05, 0.5, 0.95))
sd(boot_samples_wakeLoco)
mean(boot_samples_wakeLoco)
obs.cor.wakeLoco

#################################
boot_samples_wakeWge <- replicate(1000, boot_trial_rnaiTOeachother(rnai.males.wake, rnai.males.wge))

obs.cor.wakeWge <- cor(wake.coef, wge.coef)

quantile(boot_samples_wakeWge, probs = c(0.05, 0.5, 0.95))
sd(boot_samples_wakeWge)
mean(boot_samples_wakeWge)
obs.cor.wakeWge

#################################
boot_samples_wakeBtn <- replicate(1000, boot_trial_rnaiTOeachother(rnai.males.wake, rnai.males.btn))

obs.cor.wakeBtn <- cor(wake.coef, btn.coef)

quantile(boot_samples_wakeBtn, probs = c(0.05, 0.5, 0.95))
sd(boot_samples_wakeBtn)
mean(boot_samples_wakeBtn)
obs.cor.wakeBtn

#################################
boot_samples_wakeEf6a <- replicate(1000, boot_trial_rnaiTOeachother(rnai.males.wake, rnai.males.ef6a))

obs.cor.wakeEf6a <- cor(wake.coef, ef6a.coef)

quantile(boot_samples_wakeEf6a, probs = c(0.05, 0.5, 0.95))
sd(boot_samples_wakeEf6a)
mean(boot_samples_wakeEf6a)
obs.cor.wakeEf6a

#################################
boot_samples_wakeLoco <- replicate(1000, boot_trial_rnaiTOeachother(rnai.males.wake, rnai.males.loco))

obs.cor.wakeLoco <- cor(wake.coef, loco.coef)

quantile(boot_samples_wakeLoco, probs = c(0.05, 0.5, 0.95))
sd(boot_samples_wakeLoco)
mean(boot_samples_wakeLoco)
obs.cor.wakeLoco


