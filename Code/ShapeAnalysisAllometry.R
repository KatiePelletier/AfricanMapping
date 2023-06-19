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

#save(all.male, file = "../Data/allMalesForThesisAnalysis.rda")


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

# cross.males.morphoNull <- procD.lm(shape ~ logCS + line:logCS, 
#                                    data = gdf.cross.males)

#Cs:line is actually the allo null, as pointed out by ID
cross.males.alloNull <- procD.lm(shape ~ logCS + line, 
                                 data = gdf.cross.males)


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

anova(cross.males.fullMod)

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

#This is the model I really want to fit but it wont let me. 
# lmer(as.matrix(all.male[,1:30]) ~ logCS_c*population + (1|population:line),
#      data = all.male)

VecCorBoot_allo <- function(data, indices, vars = 1:30) {
  d <- data[indices,] # allows boot to select sample
  mlm.mod1 <- lm(as.matrix(d[,vars]) ~ logCS_c*population ,
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
VecBoot_results <- boot(data = all.male,
                        statistic = VecCorBoot_allo, R = 2000,
                        strata = as.factor(all.male$population))


#this is a matrix of the results where here line is an observaton and the cols are the 6 statistics. in the order (cor1.2, cor1.3, cor2.3, mag1, mag2, mag3)
vec.boot.mat <- data.frame(VecBoot_results[["t"]])


#the first line is the observation and then the rest are permutations.
#1 = LA, 2 = HA, 3 = HA
#cors 
quantile(vec.boot.mat[,1], probs= c(0.05, 0.95))
quantile(vec.boot.mat[,2], probs= c(0.05, 0.95))
quantile(vec.boot.mat[,3], probs= c(0.05, 0.95))

#mags 
quantile(vec.boot.mat[,4], probs= c(0.05, 0.95))
quantile(vec.boot.mat[,5], probs= c(0.05, 0.95))
quantile(vec.boot.mat[,6], probs= c(0.05, 0.95))

#I want to make a 

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
  stat_eye(fill = "grey", alpha = 0.6, 
                    p_limits = c(0.05, 0.95), 
           show_interval = FALSE,show_point = FALSE) +
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
  #theme(legend.position = "none") +
  ylab("Predicted shape from shape ~ size*line regression") +
  xlab("Size")



reg.score.plot <- ggplot(allo.score.dat, aes(x = logCS, y = ShapeScore, col = Genotype)) +
  geom_point(alpha = 0.2) +
  #scale_colour_manual(values = c("black", "grey45", "grey75",  "orange", "orange","orange", "blue4", "blue3", "blue")) +
  theme_classic() +
  ylab("ShapeScore from shape ~ size*line regression") +
  xlab("Size")


mag.plot <- ggplot(mag.boot.long , aes(x = group, y = stat)) +
  stat_eye(fill = "grey", alpha = 0.6, 
           p_limits = c(0.05, 0.95), 
           show_interval = FALSE,show_point = FALSE) +
  geom_point(data = mag.obs, col = "red") +
  xlab("") +
  ylab("Vector Magnitude") +
  theme_classic() +
  theme(text = element_text(size = 20))


vec.boot.mat[1,]

png("../Figures/allometryVecPlots.png", width = 2000, height = 1000, units = "px", res = 300)
plot_grid(cor.plot, mag.plot, labels = c("A", "B"))
dev.off()


###############################################################################
#plotting allometry with ID functions. 

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


reg.allo.slopes2 <- reg.allo.slopes + ylab("PC1 of Fitted Values")
reg.score.plot2 <- reg.score.plot + theme(legend.position = "none") + ylab("Shape Score")

plot_grid(reg.allo.slopes2, reg.score.plot2, wing.blurs, nrow =3, rel_heights = c(1, 1, 0.5), labels = c("A", "B", "C"))



png("../Figures/allometryCrosses.png", width = 3000, height = 2500, units = "px", res= 300)
plot_grid(reg.allo.slopes2, reg.score.plot2, wing.blurs, nrow =3, rel_heights = c(1, 1, 0.5), labels = c("A", "B", "C"))
dev.off()

