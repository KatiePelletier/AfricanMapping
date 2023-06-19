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


panel.cor <- function(x, y) {
  r <- cor(x, y)
  par( usr =c(0, 1, 0, 1))
  Cor <- formatC(c(r, 0.123456789), format = "f", digits=2)[1]
  text(0.5, 0.5, paste(Cor), cex=1.5)
}
#I have these superimposed in other scripts so I'm just reading them all in here

zi192ef81 <- read.csv("../Data/zi192ef81_15pt_superimposed.csv")
zi192ef81$line <- "Zi192Ef81"

zi192ef43 <- read.csv("../Data/zi192ef43_15pt_superimposed.csv")
#there is a problem with the labels here. 
zi192ef43$line <- "Zi192Ef43"

zi192ef96 <- read.csv("../Data/temp_zi192ef96_15pt_superimposed.csv")
zi192ef96$line <- "Zi192Ef96"

parents <- read.csv("../Data/pureWings_15lmk_superimposed.csv")

#Main questions here are: 
# 1) Do these cross populations have similar matrix shapes or are we getting diffrent shapes with each cross (because unique genetics). This can be used for inference about the fact that this is polygenic (and epistasis). I can also look at matrix eccectricinty here because large alleles or many corralated alleles should reduce this measure while many alleles in many directions would increase this measure. This can't all be genetic either because I think there might be a lot of E in these (hence the small G mag)
# 2) Does every cross have the same amount of variation? (size of the matrix). How does this relate to similarity for the parents?
# 3) how does PC1 compare between crosses? Does this always allign with the change vector? This is also linked with question 1. 

#I think I want to use only one sex for right now. Males first.
zi192ef96 <- dplyr::select(zi192ef96, -c("perp", "col", "plate", "fly"))

cross.shape <- rbind(zi192ef43, zi192ef81, zi192ef96)
str(cross.shape)

cross.male <- filter(cross.shape, sex == "m")
pure.male <- filter(parents, sex == "M")

all.male <- rbind(cross.male, pure.male)

save(all.male, file = "../Data/allMalesForThesisAnalysis.rda")

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
#I want to do this with the parents included so this is easier to understand
#first, I want to remove the other zambian lines
all.male.zi192 <- (all.male %>% filter(line != "ZI251") %>%
                     filter(line != "ZI418"))



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
cross.males.alloNull <- procD.lm(shape ~ logCS + line, 
                                   data = gdf.cross.males)



#png("../Figures/CrossAlometryPlot.png")
plotAllometry(cross.males.fullMod,
              size = all.male$logCS,
              logsz = FALSE,
              col = c("grey", "grey", "grey", "black", "black","black","red", "green", "blue")[as.factor(all.male$line)])
#dev.off()

crap_ian <- plotAllometry(cross.males.fullMod,
              size = all.male$logCS,
              logsz = FALSE,
              method = "RegScore",
              col = c("grey", "grey", "grey", "black", "black","black","red", "green", "blue")[as.factor(all.male$line)])

#I think that this is the line I want in my df?
allo.pred.line <- crap_ian$PredLine

allo.pred.dat <- data.frame(allo.pred.line, all.male$logCS, all.male$line)
names(allo.pred.dat) <- c("PredictedShape", "logCS", "Genotype")

allo.score.dat <- data.frame(crap_ian$RegScore, all.male$logCS, all.male$line)
names(allo.score.dat) <- c("ShapeScore", "logCS", "Genotype")


reg.allo.slopes <- ggplot(allo.pred.dat, aes(x = logCS, y = PredictedShape, col = Genotype)) +
  geom_point(alpha = 0.3) +
  #scale_colour_manual(values = c("black", "grey45", "grey75", "orange", "orange","orange","blue4", "blue3", "blue")) +
  theme_classic() +
  theme(legend.position = "none") +
  ylab("Predicted shape from shape ~ size*line regression") +
  xlab("Size")



reg.score.plot <- ggplot(allo.score.dat, aes(x = logCS, y = ShapeScore, col = Genotype)) +
  geom_point(alpha = 0.2) +
  #scale_colour_manual(values = c("black", "grey45", "grey75",  "orange", "orange","orange", "blue4", "blue3", "blue")) +
  theme_classic() +
  ylab("ShapeScore from shape ~ size*line regression") +
  xlab("Size")

png("../Figures/AllometricRelationships.png", width =2500, height = 1000, units = "px",res = 300)
plot_grid(reg.allo.slopes, reg.score.plot, labels = c("A", "B"), rel_widths = c(0.9, 1))
dev.off()

#these are the actual vector cors
allo.pair <- pairwise(cross.males.fullMod, cross.males.alloNull,
                      covariate = all.male.zi192$logCS,
                      groups = gdf.cross.males$line)

#This makes sense. Want to get the upper and lower bounds of these CIs
summary(allo.pair, test.type = "VC", angle.type = "deg")
#The pure lines have something going on as I have suspected.
summary(allo.pair, test.type = "DL") # magnitude of allometry vectors

#I want to get the resampling of these to get CI.
names(allo.pair)

crap_crap <- allo.pair$slopes.vec.cor[[1]] # I think observed
crap_crap_crap <- allo.pair$slopes.vec.cor[[2]]  # first permutation

#Now I need to write a function to extract all the permutations that I am interested from this.

#not particularly helpful. Not going to lie. Maybe I can just get the quantiles out of here and skip all the other steps
#this makes a data frame where the rows are observations and the cols are each cor vec.
ef43cors <- data.frame(t(sapply(allo.pair$slopes.vec.cor, function(x) x[1,])))
quantile(ef43cors$EF81[2:1000], probs = c(0.05, 0.95))
ef43cors$EF81[1] #observed
quantile(ef43cors$EF96[2:1000], probs = c(0.05, 0.95))
ef43cors$EF96[1]
quantile(ef43cors$ZI192[2:1000], probs = c(0.05, 0.95))
ef43cors$ZI192[1]
quantile(ef43cors$Zi192Ef43[2:1000], probs = c(0.05, 0.95))
ef43cors$Zi192Ef43[1]
quantile(ef43cors$Zi192Ef81[2:1000], probs = c(0.05, 0.95))
ef43cors$Zi192Ef81[1]
quantile(ef43cors$Zi192Ef96[2:1000], probs = c(0.05, 0.95))
ef43cors$Zi192Ef96[1]

ef81cors <- data.frame(t(sapply(allo.pair$slopes.vec.cor, function(x) x[2,])))
# quantile(ef81cors$EF43[2:1000], probs = c(0.05, 0.95))
# ef81cors$EF43[1] #observed
quantile(ef81cors$EF96[2:1000], probs = c(0.05, 0.95))
ef81cors$EF96[1]
quantile(ef81cors$ZI192[2:1000], probs = c(0.05, 0.95))
ef81cors$ZI192[1]
quantile(ef81cors$Zi192Ef43[2:1000], probs = c(0.05, 0.95))
ef81cors$Zi192Ef43[1]
quantile(ef81cors$Zi192Ef81[2:1000], probs = c(0.05, 0.95))
ef81cors$Zi192Ef81[1]
quantile(ef81cors$Zi192Ef96[2:1000], probs = c(0.05, 0.95))
ef81cors$Zi192Ef96[1]

ef96cors <- data.frame(t(sapply(allo.pair$slopes.vec.cor, function(x) x[3,])))
# quantile(ef96cors$EF43[2:1000])
# ef96cors$EF43[1] #observed
# quantile(ef96cors$EF81[2:1000])
# ef96cors$EF81[1]
quantile(ef96cors$ZI192[2:1000], probs = c(0.05, 0.95))
ef96cors$ZI192[1]
quantile(ef96cors$Zi192Ef43[2:1000], probs = c(0.05, 0.95))
ef96cors$Zi192Ef43[1]
quantile(ef96cors$Zi192Ef81[2:1000], probs = c(0.05, 0.95))
ef96cors$Zi192Ef81[1]
quantile(ef96cors$Zi192Ef96[2:1000], probs = c(0.05, 0.95))
ef96cors$Zi192Ef96[1]

zi192cors <- data.frame(t(sapply(allo.pair$slopes.vec.cor, function(x) x[4,])))
# quantile(zi192cors$EF43[2:1000])
# zi192cors$EF43[1] #observed
# quantile(zi192cors$EF96[2:1000])
# zi192cors$EF96[1]
# quantile(zi192cors$EF81[2:1000])
# zi192cors$EF81[1]
quantile(zi192cors$Zi192Ef43[2:1000], probs = c(0.05, 0.95))
zi192cors$Zi192Ef43[1]
quantile(zi192cors$Zi192Ef81[2:1000], probs = c(0.05, 0.95))
zi192cors$Zi192Ef81[1]
quantile(zi192cors$Zi192Ef96[2:1000], probs = c(0.05, 0.95))
zi192cors$Zi192Ef96[1]


zi192ef43cors <- data.frame(t(sapply(allo.pair$slopes.vec.cor, function(x) x[5,])))
# quantile(zi192ef43cors$EF43[2:1000])
# zi192ef43cors$EF43[1] #observed
# quantile(zi192ef43cors$EF96[2:1000])
# zi192ef43cors$EF96[1]
# quantile(zi192ef43cors$EF81[2:1000])
# zi192ef43cors$EF81[1]
# quantile(zi192ef43cors$ZI192[2:1000])
# zi192ef43cors$ZI192[1]
quantile(zi192ef43cors$Zi192Ef81[2:1000], probs = c(0.05, 0.95))
zi192ef43cors$Zi192Ef81[1]
quantile(zi192ef43cors$Zi192Ef96[2:1000], probs = c(0.05, 0.95))
zi192ef43cors$Zi192Ef96[1]

zi192ef81cors <- data.frame(t(sapply(allo.pair$slopes.vec.cor, function(x) x[6,])))
# quantile(zi192ef81cors$EF43[2:1000])
# zi192ef81cors$EF43[1] #observed
# quantile(zi192ef81cors$EF96[2:1000])
# zi192ef81cors$EF96[1]
# quantile(zi192ef81cors$EF81[2:1000])
# zi192ef81cors$EF81[1]
# quantile(zi192ef81cors$ZI192[2:1000])
# zi192ef81cors$ZI192[1]
# quantile(zi192ef81cors$Zi192Ef43[2:1000])
# zi192ef81cors$Zi192Ef43[1]
quantile(zi192ef81cors$Zi192Ef96[2:1000], probs = c(0.05, 0.95))
zi192ef81cors$Zi192Ef96[1]


zi192ef96cors <- data.frame(t(sapply(allo.pair$slopes.vec.cor, function(x) x[7,])))
# quantile(zi192ef96cors$EF43[2:1000])
# zi192ef96cors$EF43[1] #observed
# quantile(zi192ef96cors$EF96[2:1000])
# zi192ef96cors$EF96[1]
# quantile(zi192ef96cors$EF81[2:1000])
# zi192ef96cors$EF81[1]
# quantile(zi192ef96cors$ZI192[2:1000])
# zi192ef96cors$ZI192[1]
# quantile(zi192ef96cors$Zi192Ef81[2:1000])
# zi192ef96cors$Zi192Ef81[1]
# quantile(zi192ef96cors$Zi192Ef43[2:1000])
# zi192ef96cors$Zi192Ef43[1]

#Because the geomorph null is cor == 0, which is not very helpful for allometry vectors. I want get bootstrapped estimates for the effects.
#Also going to get magnitude while I'm at it because I will have to talk about why its so much bigger in the pure lines.
#for this, I want to use the code I adapted from MP and MP's code to do this

#setting up the model first.
all.male$population <- relevel(as.factor(all.male$population), ref  = "LA")

mlm.mod_test <- lm(as.matrix(all.male[,1:30]) ~ logCS_c*population + population:line,
                   data = all.male)

coef(mlm.mod_test)


#size term
obs.la.allo <- coef(mlm.mod_test)[2,]
#size term plus interaction term
obs.cross.allo <- coef(mlm.mod_test)[2,] + coef(mlm.mod_test)[5,]
#size term plus interaction term
obs.ha.allo <- coef(mlm.mod_test)[2,] + coef(mlm.mod_test)[6,]

#the observed magnitude of allometry is just the length of these vectors
sqrt(sum(obs.la.allo^2)) #0.1123534
sqrt(sum(obs.cross.allo ^2)) #0.06839164
sqrt(sum(obs.ha.allo^2)) # 0.1502244


VecCorBoot_allo <- function(data, indices, vars = 1:30) {
  d <- data[indices,] # allows boot to select sample
  mlm.mod1 <- lm(as.matrix(d[,vars]) ~ logCS_c*population + population:line,
                 data = d)
  coef1 <- mlm.mod1$coef[2,]
  coef2 <- mlm.mod1$coef[2,] + mlm.mod1$coef[5,]
  coef3 <- mlm.mod1$coef[2,] + mlm.mod1$coef[6,]
  cor1.2 <- cor(coef1, coef2)
  cor1.3 <- cor(coef1, coef3)
  cor2.3 <- cor(coef2, coef3)

  mag1 <- sqrt(sum(coef1 ^2))
  mag2 <- sqrt(sum(coef2 ^2))
  mag3 <- sqrt(sum(coef3 ^2))

  return(c(cor1.2, cor1.3, cor2.3, mag1, mag2, mag3))
}

#checking to make sure it runs and matches what I think it should.
VecCorBoot_allo(all.male.zi192)


#I think I want to sample within populatuion. about equal within population.
VecBoot_results <- boot(data = all.male.zi192,
                        statistic = VecCorBoot_allo, R = 2000,
                        strata = as.factor(all.male.zi192$population))


#this is a matrix of the results where here line is an observaton and the cols are the 6 statistics. in the order (cor1.2, cor1.3, cor2.3, mag1, mag2, mag3)
vec.boot.mat <- data.frame(VecBoot_results[["t"]])

names(vec.boot.mat) <- c("LA.cross", "LA.HA", "HA.cross", "LA", "cross", "HA")

cor.boot.long <- (vec.boot.mat %>%
                    dplyr::select(LA.cross, LA.HA, HA.cross) %>%
                    pivot_longer(everything(), names_to = "comparison",
                                 values_to = "stat"))

cor.obs <- t(vec.boot.mat[1,1:3])
cor.obs  <- data.frame(cor.obs, rownames(cor.obs))
names(cor.obs) <- c("stat", "comparison")


#the first line is the observation and then the rest are permutations.
#1 = LA, 2 = cross, 3 = HA

#PROBLEMATIC
cor.plot <- ggplot(cor.boot.long , aes(x = comparison, y = stat)) +
  geom_violin(fill = "grey", alpha = 0.6) +
  geom_point(data = cor.obs, col = "red") +
  xlab("") +
  ylab("Correlation") +
  theme_classic() +
  theme(text = element_text(size = 20))

vec.boot.mat[1,]


mag.boot.long <- (vec.boot.mat %>%
                    dplyr::select(LA, cross, HA) %>%
                    pivot_longer(everything(), names_to = "group",
                                 values_to = "stat"))
mag.obs <- t(vec.boot.mat[1,4:6])
mag.obs  <- data.frame(mag.obs, rownames(mag.obs))
names(mag.obs) <- c("stat", "group")

#the first line is the observation and then the rest are permutations.
#1 = LA, 2 = cross, 3 = HA

quantile(vec.boot.mat[,4], probs= c(0.05, 0.95))
quantile(vec.boot.mat[,5], probs= c(0.05, 0.95))
quantile(vec.boot.mat[,6], probs= c(0.05, 0.95))

mag.plot <- ggplot(mag.boot.long , aes(x = group, y = stat)) +
  geom_violin(fill = "grey", alpha = 0.6) +
  geom_point(data = mag.obs, col = "red") +
  xlab("") +
  ylab("Vector Magnitude") +
  theme_classic() +
  theme(text = element_text(size = 20))


vec.boot.mat[1,]

png("../Figures/allometryVecPlots.png")
plot_grid(cor.plot, mag.plot, labels = c("A", "B"))
dev.off()



# check the sd for size for each group

# From a previous run without the parents. 
# r      angle  UCL (95%)          Z Pr > angle
# Zi192Ef43:Zi192Ef81 0.9999422 0.01075262 0.01140389 -0.2833331      0.614
# Zi192Ef43:Zi192Ef96 0.9998948 0.01450527 0.01532752 -0.6839989      0.762
# Zi192Ef81:Zi192Ef96 0.9999159 0.01296874 0.01326854  0.5663238      0.294

#Or I can do this outside of geomorph. Following an old tutorial ID wrote for MP and trying to. It looks to me like Maria used lm and not PCs for this. 

###############################################################################
#plotting allometry with ID functions. 

#for rounding on the cor

#for getting lines 
wing_links2 <- c(1, 7, 
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

wing_links_mat <- matrix(wing_links2, 
                         byrow = TRUE,
                         nrow = (length(wing_links2)/2),
                         ncol = 2)

wing_links_mat <- cbind(wing_links_mat, NA)

wing_links_flat <- as.vector(t(wing_links_mat), 
                             mode = "integer")

index_x <- seq.int(from = 1, to = 29, by = 2)
index_y <- seq.int(from = 2, to = 30, by = 2)

line_index_x <- (wing_links_flat*2 - 1) # landmark number, indexing column number in matrix for x coordinate

line_index_y <- wing_links_flat*2

#plotting within gentotype for the largest and smallest.
#Setting the plot colours. Right now, red blue and gray. 
point_colours_anisotropy <- c(rgb(0,0,1, 0.2),
                              rgb(0.6,0.6,0.6, 0.05),
                              rgb(1,0,0, 0.2))

#taking the data and making the smallest, largest and middlest groups. 
zi192ef81.males <- filter(cross.male, line == "Zi192Ef81")
zi192ef81.males.lmk <- as.matrix(zi192ef81.males[,1:30])
zi192ef81.males_logCS <- zi192ef81.males$logCS

zi192ef81.males_logCS_breaks <- cut(zi192ef81.males_logCS, include.lowest = T,
                       breaks = quantile(zi192ef81.males_logCS, 
                                         probs = c(0, 0.1,0.9, 1)),
                       labels = c("small","medium","high"))

png("../Figures/zi192ef81_anisotropyPlot.png")
matplot(x = t(zi192ef81.males.lmk[,line_index_x]),
        y = t(zi192ef81.males.lmk[,line_index_y]),
        lty = 1, type = "l", lwd = 0.3,
        xlab = "", ylab = "", asp = 1,
        frame.plot = F, ann = F, axes = F,
        col = point_colours_anisotropy[zi192ef81.males_logCS_breaks])

matpoints(x = t(zi192ef81.males.lmk[,index_x]),
          y = t(zi192ef81.males.lmk[,index_y]),
          pch = c(20, 20, 20)[zi192ef81.males_logCS_breaks], 
          cex = c(0.5, 0.5, 0.5)[zi192ef81.males_logCS_breaks],
          col = point_colours_anisotropy[zi192ef81.males_logCS_breaks],
          xlab = "", ylab = "", asp = 1,
          frame.plot = F, ann = F, axes = F)
dev.off()

#this is a smaller number of wings. 
zi192ef96.males <- filter(cross.male, line == "Zi192Ef96")
zi192ef96.males.lmk <- as.matrix(zi192ef96.males[,1:30])
zi192ef96.males_logCS <- zi192ef96.males$logCS

png("../Figures/zi192ef96_anisotropyPlot.png")
zi192ef96.males_logCS_breaks <- cut(zi192ef96.males_logCS, include.lowest = T,
                                    breaks = quantile(zi192ef96.males_logCS, 
                                                      probs = c(0, 0.1,0.9, 1)),
                                    labels = c("small","medium","high"))

matplot(x = t(zi192ef96.males.lmk[,line_index_x]),
        y = t(zi192ef96.males.lmk[,line_index_y]),
        lty = 1, type = "l", lwd = 0.3,
        xlab = "", ylab = "", asp = 1,
        frame.plot = F, ann = F, axes = F,
        col = point_colours_anisotropy[zi192ef96.males_logCS_breaks])

matpoints(x = t(zi192ef96.males.lmk[,index_x]),
          y = t(zi192ef96.males.lmk[,index_y]),
          pch = c(20, 20, 20)[zi192ef96.males_logCS_breaks], 
          cex = c(0.5, 0.5, 0.5)[zi192ef96.males_logCS_breaks],
          col = point_colours_anisotropy[zi192ef96.males_logCS_breaks],
          xlab = "", ylab = "", asp = 1,
          frame.plot = F, ann = F, axes = F)

dev.off()

zi192ef43.males <- filter(cross.male, line == "Zi192Ef43")
zi192ef43.males.lmk <- as.matrix(zi192ef43.males[,1:30])
zi192ef43.males_logCS <- zi192ef43.males$logCS

zi192ef43.males_logCS_breaks <- cut(zi192ef43.males_logCS, include.lowest = T,
                                    breaks = quantile(zi192ef43.males_logCS, 
                                                      probs = c(0, 0.1,0.9, 1)),
                                    labels = c("small","medium","high"))

png("../Figures/zi192ef43_anisotropyPlot.png")

matplot(x = t(zi192ef43.males.lmk[,line_index_x]),
        y = t(zi192ef43.males.lmk[,line_index_y]),
        lty = 1, type = "l", lwd = 0.3,
        xlab = "", ylab = "", asp = 1,
        frame.plot = F, ann = F, axes = F,
        col = point_colours_anisotropy[zi192ef43.males_logCS_breaks])

matpoints(x = t(zi192ef43.males.lmk[,index_x]),
          y = t(zi192ef43.males.lmk[,index_y]),
          pch = c(20, 20, 20)[zi192ef43.males_logCS_breaks], 
          cex = c(0.5, 0.5, 0.5)[zi192ef43.males_logCS_breaks],
          col = point_colours_anisotropy[zi192ef43.males_logCS_breaks],
          xlab = "", ylab = "", asp = 1,
          frame.plot = F, ann = F, axes = F)
dev.off()

#going to make a figure that includes this and the allometry vectors for the paper. 

zi192ef96.raw <- image_read("../Figures/zi192ef96_anisotropyPlot.png")
zi192ef96.crop <- image_trim(zi192ef96.raw)
zi192ef96.aniso <- ggdraw() + draw_image(zi192ef96.crop)

zi192ef81.raw <- image_read("../Figures/zi192ef81_anisotropyPlot.png")
zi192ef81.crop <- image_trim(zi192ef81.raw)
zi192ef81.aniso <- ggdraw() + draw_image(zi192ef81.crop)

zi192ef43.raw <- image_read("../Figures/zi192ef43_anisotropyPlot.png")
zi192ef43.crop <- image_trim(zi192ef43.raw)
zi192ef43.aniso <- ggdraw() + draw_image(zi192ef43.crop)

wing.blurs <- plot_grid(zi192ef43.aniso, zi192ef81.aniso, zi192ef96.aniso, 
                        nrow = 1, 
                        labels = c("Zi192 x Ef43", "Zi192 x Ef81", "Zi192 x Ef96"))

allo.plots <- plot_grid(reg.allo.slopes, reg.score.plot, labels = c("A", "B"), rel_widths = c(0.85, 1))

png("../Figures/allometryCrosses.png", width = 2500, height = 1700, units = "px", res= 300)
plot_grid(allo.plots, wing.blurs, nrow =2, rel_heights = c(1, 0.5))
dev.off()


###############################################################################
#Is the size actually significantly different?

#first, I want to remove the other zambian lines
all.male.zi192 <- (all.male %>% filter(line != "ZI251") %>%
                     filter(line != "ZI418"))

#these numbers make sense. But the model fits are wild, probably because I have a shitty plan there. 
(all.male.zi192 %>% group_by(line) %>% summarise(mean_size = mean(CS)))

#I think I want to add estimates with CI on top of this?
ggplot(all.male.zi192, aes(x = line, y = logCS, col = population)) + 
  geom_violin()

with(all.male.zi192, table(line))


#this treats everything as an independent genotype, which they are not, because there is the added 
#Is this fine when there is so much more data for mine?
size.cross.mod <- lm(CS ~ line, data = all.male.zi192)

summary(size.cross.mod)
anova(size.cross.mod)

#the answer is yes
plot(emmeans(size.cross.mod, ~ line))
#because I don't want to scroll up. I know this is also in the summary. This is more here as a place holder for when I have more crosses and don't feel like doing this. 
pairs(emmeans(size.cross.mod, ~line))


size.est <- data.frame(emmeans(size.cross.mod, ~ line))

plot.order <- c("ZI192", "Zi192Ef43", "EF43", "Zi192Ef81", "EF81", "Zi192Ef96", "EF96")

all.male.zi192$line <- factor(all.male.zi192$line, levels = plot.order)
size.est$line <- factor(size.est$line, levels = plot.order)

cross.size.plot <- ggplot(all.male.zi192, aes(x = line, y = CS)) + 
  geom_violin(aes(fill = population)) + 
  scale_fill_manual(values=c('blue', 'grey66', 'grey32')) +
  geom_errorbar(data = size.est, 
                aes(y = emmean, ymin = lower.CL, ymax = upper.CL, width = 0.25, 
                    alpha = 0.3)) + 
  geom_point(data = size.est, aes(y = emmean, 
             alpha = 0.3)) + 
  theme_classic() + 
  theme(legend.position = "none") + 
  ylab("Size") + 
  xlab("")

# png("../Figures/size_est_parentsNext.png")
# cross.size.plot
# dev.off()

###################################################

#back to what I was doing before this CHAOS CODE break. 
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
  geom_quasirandom() +
  labs(x = "", y = "Phenotypic Integration") 

ggplot(mat_stats, aes(y =  PC1, x = population )) +
  geom_quasirandom() +
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
  geom_violin(alpha = 0.5) +  
  geom_point(data = mat_stats, aes(x = line, y = PC1)) + 
  theme_classic() + 
  ylab("Propotrion Variance PC1") + 
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
  geom_violin(alpha = 0.5) +  
  geom_point(data = mat_stats, aes(x = line, y = ICV)) + 
  theme_classic() + 
  ylab("Eccentricity") + 
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
  geom_violin(alpha = 0.5) + 
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
  geom_violin(alpha = 0.5) +  
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






