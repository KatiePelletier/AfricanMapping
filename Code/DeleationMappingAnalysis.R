#Looking at the deleation mapping analysis for the peak on chromosome 3R
library(tidyverse)
library(abind)
library(geomorph)
library(broom)
library(ggdist)
library(distributional)
library(broom.mixed)
library(car)
library(cowplot)
library(emmeans)
library(lme4)

wing15.links <- c(1, 7, 
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


projFunction <- function(x, y) {
  scalarProj <- (x %*% y) / norm(y, type = "2")
  return(scalarProj)
}

###########################################
#cleaning 
###########################################
raw.dat <- read.delim("../Data/delMappling15lmkNEW.txt", header = F)

colnames(raw.dat) <- c("image", "x1", "y1", "x2", "y2", "x3", "y3", "x4", "y4","x5", "y5","x6", "y6","x7", "y7","x8", "y8","x9", "y9","x10", "y10","x11", "y11", "x12", "y12","x13", "y13", "x14", "y14","x15", "y15")


#cleaning up the image id
raw.dat <- separate(raw.dat, image, into = c("perp", "delLine", "background", "rep", "sex", "fly"))

#Now I need to superimpose everything.
#again using the same common mean I have been using all along.
common.mean <- read.csv("../Data/afMap_commonMean.csv")

cord <- as.matrix(raw.dat[,7:36])
shape <- arrayspecs(cord, 15, 2)

shape2 <- abind(common.mean, shape, along=3)

gdf <- geomorph.data.frame(landmarks = shape2)

#0 breaks this. needs to be 1 iter. This should constrain it to be in a common shape space with that mean.
del.shape <- gpagen(gdf$landmarks, ProcD = T, max.iter = 1)

#the small one is the common mean. Is that ok? Double check with Ian about the scaling here.
#looks normal enough.
hist(del.shape$Csize)

#that is a wing.
plot(del.shape, links = wing15.links)

#outliers.
#284 is bad. 556 is the one where I swapped the points and I was worried about
#the rest look fine 
plotOutliers(del.shape$coords)
#check wings
#plotOutliers(del.shape$coords, inspect.outliers = TRUE)

#plotOutliers(del.shape$coords[,,-c(284, 556)], inspect.outliers = TRUE)
plotOutliers(del.shape$coords[,,-c(284, 556)])

crap <- as.data.frame(two.d.array(del.shape$coords))
#the first line is the common mean that I don't actually care about
wings.del <- crap[2:nrow(crap),]
wings.del$CS <- del.shape$Csize[-1]
wings.del$delLine <- raw.dat$delLine
wings.del$background <- raw.dat$background
wings.del$sex <- raw.dat$sex
wings.del$rep <- raw.dat$rep

#now I need to drop out the bad wings. 
#subtract 1 from index because it is the common mean in that array 283, 555

wings.del <- wings.del[-c(283, 555), ]

#got to fix these names so they are all the same.
with(wings.del, table(delLine, background))

wings.del$background <- gsub("eF96", "ef96", wings.del$background)
wings.del$background <- gsub("Zi192", "zi192", wings.del$background)
wings.del$background <- gsub("ef46", "ef43", wings.del$background)
wings.del$background <- gsub("Zi360", "zi360", wings.del$background)



#Some of these will really be a problem.
#These are all males too.
#ef96 DDC and zi251 ExC control numbers are a little low. 
with(wings.del, table(delLine, background))
#Talk to ID. I might want to ignore the vial level effects if I have to for this
#alot of zi251 vials failed. Makes sense its a shitty line
with(wings.del, table(delLine, background, rep))

#need a dictionary of which dels are in which population
del.dict <- read.delim("../Data/delLineDict.txt")

wings.del <- left_join(wings.del, del.dict, by = "delLine")

# write.csv(wings.del, "../Data/delMapping_15lmk_superimposed.csv", row.names = FALSE)

wings.del$logCS <- log2(wings.del$CS)
wings.del$logCS_c <- wings.del$logCS - mean(wings.del$logCS)


#grabbing this so that I can get the shape score. 
parents.wings <- read.csv("../Data/pureWings_15lmk_superimposed.csv")

#Joint pca of both datasets. 
#need to fix these col names
names(parents.wings) <- c(names(wings.del[,1:30]), names(parents.wings[,31:36]))

total.lm <- rbind(parents.wings[,c(1:30, 35)], wings.del[,c(1:30, 37)])

totalwing.pca <- prcomp(total.lm)

parents.wings.pcs <- data.frame(parents.wings, totalwing.pca$x[1:633,1:27])
del.wings.pcs <- data.frame(wings.del, totalwing.pca$x[634:2621,1:27])

#ok now that looks good. something was weird with the extraction before. 
plot(totalwing.pca$x[,1], c(parents.wings$logCS, wings.del$logCS))

#ok. This is doing what I think it is. 
ggplot(parents.wings.pcs, aes(x = logCS, y = PC1)) + geom_point()
ggplot(del.wings.pcs, aes(x = logCS, y = PC1)) + geom_point()

#getting the high-low vector using PC 2-27 
 
highmale <- colMeans((parents.wings.pcs %>% filter(sex == "M", population == "HA"))[,38:63])

lowmale <- colMeans((parents.wings.pcs %>% filter(sex == "M", population == "LA"))[,38:63])

altVec <- highmale - lowmale

del.wings.pcs$shapescore <- projFunction(as.matrix(del.wings.pcs[,40:65]), altVec)


wings.del$shapescore <- del.wings.pcs$shapescore
###########################################
#PCA Analysis of just the deletion lines for a quick look 
###########################################

del.wings.onlyPC <- data.frame(wings.del, prcomp(wings.del[,1:30])$x)

#there is some allometry in here.
#much more sep on background than delLine
#Overall, this is a lot less variation than I would expect. Maybe these all have no wing shape genes? FUCK. 
#There is some background level effects here. 
ggplot(del.wings.onlyPC, aes(x = CS, y = PC1, col = background)) + 
  geom_point(alpha = 0.5)

ggplot(del.wings.onlyPC, aes(x = PC1, y = PC2, 
                             shape = background, 
                             col = pannel)) + 
  geom_point(alpha = 0.5)

#There are things here but not the most interesting. I could always do this with an allometry correction. but I honestly don't care that much. 

###########################################
#DrosDel model
###########################################
#going to consider the two panels separate because there are different controls. 

drosDel.wings <- filter(wings.del, pannel == "DrosDel")
#releving to make the control the ref level. 

drosDel.wings$delLine <- relevel(as.factor(drosDel.wings$delLine), ref = "DDC")
levels(drosDel.wings$delLine)
drosDel.wings$background <- as.factor(drosDel.wings$background)
levels(drosDel.wings$background) #"ef81"  "ef96"  "zi192"

#random effect is at the level of vial replicate. 
#get rank deficient warning. 
drosDel.mod <- lmer(shapescore ~ (CS+delLine+background)^2 + 
                      (1|delLine:background:rep), data = drosDel.wings)

#I actually think this is part of the problem. Talk to ID?
#Is this lental?
with(drosDel.wings, table(delLine, background))

#these look super little to me?
summary(drosDel.mod)
Anova(drosDel.mod)

#that one missing zi192 one is super annoying. 
#contrasts below are way easier to understand. 
# plot(emmeans(drosDel.mod, ~ background | delLine))
# 
# estimates <- tidy(emmeans(drosDel.mod, ~  delLine | background), effects="fixed", conf.int = T)
# 
# emmip(drosDel.mod, ~  background | delLine, CI = TRUE)

contrasts <- tidy(contrast(emmeans(drosDel.mod, ~  delLine | background), "trt.vs.ctrl", conf.int = TRUE, conf.level = 0.95)) # contrasts between Pops for each treatment
plot(contrast(emmeans(drosDel.mod, ~  delLine | background), "trt.vs.ctrl")) +
  geom_vline(xintercept = 0, alpha = 0.25, size = 2, color = "red")

#Not a whole lot to see here?
contrasts$dummy <- with(contrasts, interaction(background, contrast, drop = T, sep = ":"))

#NA will fuck this 
contrasts <- contrasts[complete.cases(contrasts),]

crap <- (contrasts %>% 
           separate(dummy, into = c("dummy2", "control"), sep = " - ") %>% separate(dummy2, into = c("pop", "del"), sep = ":"))


#whomp whomp. These are giant tho so I guess it is what it is.
#This is a beautiful rainbow. But super hard to look at. 
ggplot(crap , aes(x = estimate, y = del, color = background)) +
  stat_halfeye(aes(xdist = dist_student_t(df = df,mu = estimate,
                                          sigma = std.error), color =background,
                   fill = background), slab_alpha = 0.25, interval_alpha = 0.55) +
  geom_vline(xintercept = 0, lwd = 1, lty = 2, alpha = 0.3) 

#means estimates 
means.drosDel <- tidy(emmeans(drosDel.mod, ~  delLine | background), conf.int = TRUE, conf.level = 0.95)

means.drosDel.control <- filter(means.drosDel, delLine == "DDC")
means.drosDel.dels <- filter(means.drosDel, delLine != "DDC")

ggplot(means.drosDel.dels, aes(x = estimate,  y= delLine)) + 
  geom_point() + 
  geom_linerange(aes(xmin = conf.low, xmax = conf.high)) + 
  geom_vline(data = means.drosDel.control, aes(xintercept = estimate), col = "red") + 
  geom_rect(data = means.drosDel.control, aes(xmin = conf.low, xmax = conf.high, ymin = 0, ymax = Inf), fill = "red", alpha = 0.25) + 
  facet_wrap(~background) + 
  xlab("Estimated Shape Score") + 
  ylab("Deleation Line") 

###########################################
#Exel model
###########################################

exel.wings <- filter(wings.del, pannel == "Exelixis")
#releving to make the control the ref level. 

exel.wings$delLine <- relevel(as.factor(exel.wings$delLine), ref = "ExC")
levels(exel.wings$delLine)
exel.wings$background <- as.factor(exel.wings$background)
levels(exel.wings$background) #"ef81"  "ef96"  "zi192"

#random effect is at the level of vial replicate. 
#get rank deficient warning. 
exel.mod <- lmer(shapescore ~ (CS+delLine+background)^2 + 
                   (1|delLine:background:rep), data = exel.wings)

#I actually think this is part of the problem. Talk to ID?
#Is this lental? what happened?
with(exel.wings, table(delLine, background))

#these look super little to me?
summary(exel.mod)
Anova(exel.mod)

#that one missing zi192 one is super annoying. 
#contrasts below are way easier to understand. 
# plot(emmeans(exel.mod, ~ background | delLine))
# 
# estimates <- tidy(emmeans(exel.mod, ~  delLine | background), effects="fixed", conf.int = T)
# 
# emmip(exel.mod, ~  background | delLine, CI = TRUE)

contrasts.exel <- tidy(contrast(emmeans(exel.mod, ~  delLine | background), "trt.vs.ctrl", conf.int = TRUE, conf.level = 0.95)) # contrasts between Pops for each treatment
plot(contrast(emmeans(exel.mod, ~  delLine | background), "trt.vs.ctrl")) +
  geom_vline(xintercept = 0, alpha = 0.25, size = 2, color = "red")

#Not a whole lot to see here?
contrasts.exel$dummy <- with(contrasts.exel, interaction(background, contrast, drop = T, sep = ":"))

#NA will fuck this 
contrasts.exel <- contrasts.exel[complete.cases(contrasts.exel),]

crap.exel <- (contrasts.exel %>% 
           separate(dummy, into = c("dummy2", "control"), sep = " - ") %>% separate(dummy2, into = c("pop", "del"), sep = ":"))


#whomp whomp. Literally nothing to see here. No changes to shapescore. 
#What if nothing matters at all? 
#7670 looks kind of cool. 

#change plot to look more like SM paer. 
ggplot(crap.exel , aes(x = estimate, y = del, color = background)) +
  stat_halfeye(aes(xdist = dist_student_t(df = df,mu = estimate,
                                          sigma = std.error), color =background,
                   fill = background), slab_alpha = 0.25, interval_alpha = 0.55) +
  geom_vline(xintercept = 0, lwd = 1, lty = 2, alpha = 0.3) 



# ggplot(crap.exel , aes(x = estimate, y = del, color = background)) +
#   stat_halfeye(aes(xdist = dist_student_t(df = df,mu = estimate,
#                                           sigma = std.error), color =background,
#                    fill = background), slab_alpha = 0.25, interval_alpha = 0.55) +
#   geom_vline(xintercept = 0, lwd = 1, lty = 2, alpha = 0.3) 
# tidy(exel.mod, effects = "fixed", conf.int=TRUE, conf.method="profile")


#means estimates 
means.exel <- tidy(emmeans(exel.mod, ~  delLine | background), conf.int = TRUE, conf.level = 0.95)

means.exel.control <- filter(means.exel, delLine == "ExC")
means.exel.dels <- filter(means.exel, delLine != "ExC")

ggplot(means.exel.dels, aes(x = estimate,  y= delLine)) + 
  geom_point() + 
  geom_linerange(aes(xmin = conf.low, xmax = conf.high)) + 
  geom_vline(data = means.exel.control, aes(xintercept = estimate), col = "red") + 
  geom_rect(data = means.exel.control, aes(xmin = conf.low, xmax = conf.high, ymin = 0, ymax = Inf), fill = "red", alpha = 0.25) + 
  facet_wrap(~background) + 
  xlab("Estimated Shape Score") + 
  ylab("Deleation Line") 



#After meeting with ID: make reaction norm plots of the iteraction between pop:del.
drosDel.wings$population <- ifelse(grepl("^z", drosDel.wings$background), "LA", "HA")

drosDel.wings$population <- relevel(factor(drosDel.wings$population), ref = "LA")
levels(drosDel.wings$population) #"LA" "HA"


#model fit
drosDel.wings$CS.c <- drosDel.wings$CS - mean(drosDel.wings$CS)

drosDel.mod2 <- lmer(shapescore ~ (CS.c+delLine+population)^2 + (1|background:population) + 
                   (1|delLine:background:rep), data = drosDel.wings)

Anova(drosDel.mod2)
Anova(exel.mod2)

# plot(allEffects(drosDel.mod2))
# 
# plot(predictorEffect("population", drosDel.mod2))
# 
# #popBydelEffs <- 


plot(emmeans(drosDel.mod2, ~ population | delLine))

pairs(emmeans(drosDel.mod2, ~ delLine | population))

interaction.est <- tidy(emmeans(drosDel.mod2, ~ population | delLine), conf.int = TRUE, conf.level = 0.95)

#now i think I need to plot these all out individually because I can't think of another way to do it. 

#first relevel to control and change DDC to control 
interaction.est$delLine <- gsub("DDC", "Control", interaction.est$delLine )
interaction.est$delLine <- relevel(as.factor(interaction.est$delLine), ref = "Control")

interaction.est$population <- relevel(as.factor(interaction.est$population), ref = "LA")

dd25015 <- interaction.est[interaction.est$delLine == "Control" |
                             interaction.est$delLine == "25015",]


dd25015plot_legend <- ggplot(dd25015, aes(x = delLine, y = estimate, shape = population)) + 
  geom_line(aes(group = population), position=position_dodge(0.4), col = 'grey') +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, 
                position=position_dodge(0.4))  + 
  geom_point(position=position_dodge(0.4), size = 3) + 
  ylab("Shape Score") + 
  xlab("") + 
  theme_classic() + 
  theme(text = element_text(size = 20), legend.position = "bottom")

#need to filter out the regular data
del.control <- filter(drosDel.wings, delLine == "DDC")
del.control$delLine <- "Control"
delwings.25015 <- rbind(del.control, filter(drosDel.wings, delLine == "25015"))

dd25015plot <- ggplot(dd25015, aes(x = delLine, y = estimate, shape = population)) + 
  geom_line(aes(group = population), position=position_dodge(0.4), col = 'grey') +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, 
                position=position_dodge(0.4))  + 
  geom_point(position=position_dodge(0.4), size = 3) + 
  geom_point(data = delwings.25015, aes(x = delLine, y = shapescore, shape = population), alpha = 0.2, position = position_jitterdodge(dodge.width=0.4, jitter.width = 0.2)) + 
  ylab("Shape Score") + 
  xlab("") + 
  theme_classic() + 
  theme(text = element_text(size = 20), legend.position = "none")
  
dd25693 <- interaction.est[interaction.est$delLine == "Control" |
                             interaction.est$delLine == "25693",]

delwings.25693 <- rbind(del.control, filter(drosDel.wings, delLine == "25693"))

dd25693plot <- ggplot(dd25693, aes(x = delLine, y = estimate, shape = population)) + 
  geom_line(aes(group = population), position=position_dodge(0.4), col = 'grey') +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, 
                position=position_dodge(0.4))  + 
  geom_point(position=position_dodge(0.4), size = 3) + 
  geom_point(data = delwings.25693, aes(x = delLine, y = shapescore, shape = population), alpha = 0.2, position = position_jitterdodge(dodge.width=0.4, jitter.width = 0.2)) + 
  ylab("Shape Score") + 
  xlab("") + 
  theme_classic() + 
  theme(text = element_text(size = 20), legend.position = "none")


dd25727 <- interaction.est[interaction.est$delLine == "Control" |
                             interaction.est$delLine == "25727",]

delwings.25727 <- rbind(del.control, filter(drosDel.wings, delLine == "25727"))

dd25727plot <- ggplot(dd25727, aes(x = delLine, y = estimate, shape = population)) + 
  geom_line(aes(group = population), position=position_dodge(0.4), col = 'grey') +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, 
                position=position_dodge(0.4))  + 
  geom_point(position=position_dodge(0.4), size = 3) + 
  geom_point(data = delwings.25727, aes(x = delLine, y = shapescore, shape = population), alpha = 0.2, position = position_jitterdodge(dodge.width=0.4, jitter.width = 0.2)) + 
  ylab("Shape Score") + 
  xlab("") + 
  theme_classic() + 
  theme(text = element_text(size = 20), legend.position = "none")

dd26537 <- interaction.est[interaction.est$delLine == "Control" |
                             interaction.est$delLine == "26537",]

delwings.26537 <- rbind(del.control, filter(drosDel.wings, delLine == "26537"))

dd26537plot <- ggplot(dd26537, aes(x = delLine, y = estimate, shape = population)) + 
  geom_line(aes(group = population), position=position_dodge(0.4), col = 'grey') +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, 
                position=position_dodge(0.4))  + 
  geom_point(position=position_dodge(0.4), size = 3) + 
  geom_point(data = delwings.26537, aes(x = delLine, y = shapescore, shape = population), alpha = 0.2, position = position_jitterdodge(dodge.width=0.4, jitter.width = 0.2)) + 
  ylab("Shape Score") + 
  xlab("") + 
  theme_classic() + 
  theme(text = element_text(size = 20), legend.position = "none")


dd27375 <- interaction.est[interaction.est$delLine == "Control" |
                             interaction.est$delLine == "27375",]

delwings.27375 <- rbind(del.control, filter(drosDel.wings, delLine == "27375"))



dd27375plot <- ggplot(dd27375, aes(x = delLine, y = estimate, shape = population)) + 
  geom_line(aes(group = population), position=position_dodge(0.4), col = 'grey') +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, 
                position=position_dodge(0.4))  + 
  geom_point(position=position_dodge(0.4), size = 3) + 
  geom_point(data = delwings.27375, aes(x = delLine, y = shapescore, shape = population), alpha = 0.2, position = position_jitterdodge(dodge.width=0.4, jitter.width = 0.2)) + 
  ylab("Shape Score") + 
  xlab("") + 
  theme_classic() + 
  theme(text = element_text(size = 20), legend.position = "none")


legend <- get_legend(dd25015plot_legend)
drosDelpan <- plot_grid(dd25015plot, dd25693plot, dd25727plot, 
                        dd26537plot, dd27375plot, nrow = 1)

#now with the exel pannel 
exel.wings$population <- ifelse(grepl("^z", exel.wings$background), "LA", "HA")

exel.wings$population <- relevel(factor(exel.wings$population), ref = "LA")
levels(exel.wings$population) #"LA" "HA"


#model fit
exel.wings$CS.c <- exel.wings$CS - mean(exel.wings$CS)

exel.mod2 <- lmer(shapescore ~ (CS.c+delLine+population)^2 + (1|background:population) + 
                       (1|delLine:background:rep), data = exel.wings)

summary(exel.mod2)
Anova(exel.mod2)

plot(emmeans(exel.mod2, ~ population | delLine))

interaction.est.exel <- tidy(emmeans(exel.mod2, ~ population | delLine), conf.int = TRUE, conf.level = 0.95)

#now i think I need to plot these all out individually because I can't think of another way to do it. 

#first relevel to control and change DDC to control 
interaction.est.exel$delLine <- gsub("ExC", "Control", interaction.est.exel$delLine )
interaction.est.exel$delLine <- relevel(as.factor(interaction.est.exel$delLine), ref = "Control")
levels(interaction.est.exel$delLine)

interaction.est.exel$population <- relevel(as.factor(interaction.est.exel$population), 
                                           ref = "LA")

exel7670 <- interaction.est.exel[interaction.est.exel$delLine == "Control" |
                             interaction.est.exel$delLine == "7670",]

#need to filter out the regular data
exel.control <- filter(exel.wings, delLine == "ExC")
exel.control$delLine <- "Control"
exel.7670 <- rbind(exel.control, filter(exel.wings, delLine == "7670"))

exel7670plot <- ggplot(exel7670, aes(x = delLine, y = estimate, shape = population)) + 
  geom_line(aes(group = population), position=position_dodge(0.4), col = 'grey') +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, 
                position=position_dodge(0.4))  + 
  geom_point(position=position_dodge(0.4), size = 3) + 
  geom_point(data = exel.7670, aes(x = delLine, y = shapescore, shape = population), alpha = 0.2, position = position_jitterdodge(dodge.width=0.4, jitter.width = 0.2)) + 
  ylab("Shape Score") + 
  xlab("") + 
  theme_classic() + 
  theme(text = element_text(size = 20), legend.position = "none")


exel7671 <- interaction.est.exel[interaction.est.exel$delLine == "Control" |
                              interaction.est.exel$delLine == "7671",]

exel.7671 <- rbind(exel.control, filter(exel.wings, delLine == "7671"))


exel7671plot <- ggplot(exel7671, aes(x = delLine, y = estimate, shape = population)) + 
  geom_line(aes(group = population), position=position_dodge(0.4), col = 'grey') +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, 
                position=position_dodge(0.4))  + 
  geom_point(position=position_dodge(0.4), size = 3) + 
  geom_point(data = exel.7671, aes(x = delLine, y = shapescore, shape = population), alpha = 0.2, position = position_jitterdodge(dodge.width=0.4, jitter.width = 0.2)) + 
  ylab("Shape Score") + 
  xlab("") + 
  theme_classic() + 
  theme(text = element_text(size = 20), legend.position = "none")

exel7672 <- interaction.est.exel[interaction.est.exel$delLine == "Control" |
                              interaction.est.exel$delLine == "7672",]

exel.7672 <- rbind(exel.control, filter(exel.wings, delLine == "7672"))


exel7672plot <- ggplot(exel7672, aes(x = delLine, y = estimate, shape = population)) + 
  geom_line(aes(group = population), position=position_dodge(0.4), col = 'grey') +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, 
                position=position_dodge(0.4))  + 
  geom_point(position=position_dodge(0.4), size = 3) + 
  geom_point(data = exel.7672, aes(x = delLine, y = shapescore, shape = population), alpha = 0.2, position = position_jitterdodge(dodge.width=0.4, jitter.width = 0.2)) + 
  ylab("Shape Score") + 
  xlab("") + 
  theme_classic() + 
  theme(text = element_text(size = 20), legend.position = "none")



exel7740 <- interaction.est.exel[interaction.est.exel$delLine == "Control" |
                              interaction.est.exel$delLine == "7740",]

exel.7740 <- rbind(exel.control, filter(exel.wings, delLine == "7740"))


exel7740plot <- ggplot(exel7740, aes(x = delLine, y = estimate, shape = population)) + 
  geom_line(aes(group = population), position=position_dodge(0.4), col = 'grey') +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, 
                position=position_dodge(0.4))  + 
  geom_point(position=position_dodge(0.4), size = 3) +
  geom_point(data = exel.7740, aes(x = delLine, y = shapescore, shape = population), alpha = 0.2, position = position_jitterdodge(dodge.width=0.4, jitter.width = 0.2)) + 
  ylab("Shape Score") + 
  xlab("") + 
  theme_classic() + 
  theme(text = element_text(size = 20), legend.position = "none")


Exelpan <- plot_grid(exel7670plot, exel7671plot, exel7672plot, exel7740plot, nrow = 1)

legend <- get_legend(dd25015plot_legend)


png("../Figures/DeleationLine_allPts.png", height = 3000, width = 5000, units = "px", res = 300)
plot_grid(drosDelpan, Exelpan, legend, nrow = 3, 
          labels = c("DrosDel", "Exelis", ""), rel_heights = c(1,1,0.1))
dev.off()


