###Trying to look at the relationship between my selection vec and the whole population using MP wings from her codition paper###### 
#This is the data she gave me as the compleate set in 2022 (from her dropbox). 

library(tidyverse)
library(car)
library(lme4)
library(geomorph)

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

##### Edited from the WRP_functions file to round better. 
panel.cor <- function(x, y) {
  r <- round(cor(x, y), digits = 2)
  par( usr =c(0, 1, 0, 1))
  Cor <- formatC(c(r, 0.123456789), format = "f", digits = 2)[1]
  text(0.5, 0.5, paste(Cor), cex=1.5)
}

`%!in%` <- Negate(`%in%`)

source( "~/Dropbox/DworkinLabSharedMaterial/scripts/WINGPLOTSOURCE.R" )
source( "~/Dropbox/DworkinLabSharedMaterial/scripts/WRP_FUNCTIONS.R" )
source("~/Dropbox/KatiePelletier/KP_geomorphwingfunctions.R")

######

#Total superimposed wings including my purelines 
wings.raw <- read.delim("../Data/new_2022_totalSuperImposition.dat")

#the quickest check 
#This looks good to me. 
#plot is stacked. 
ggplot(wings.raw, aes(x = CS, fill = Perp)) + 
  geom_histogram()


#about 4500 wings 
mp.wings <- wings.raw[grepl("MP_new2022", wings.raw$File),]

#about 500 wings. 
kp.wings <- wings.raw[grepl("KP_pl_stuff", wings.raw$File ),]

#cleaning up the files to get info out. 
mp.clean <- (mp.wings %>%
               separate(File, into = c(paste(rep("crap", 5), seq(5,1)), "name"), 
                        sep = "\\\\") %>%
               separate(name, into = c("perp", "line", "condition", "rep", "sex", "fly"), 
                        sep= "_") %>% 
               dplyr::select(!c(paste(rep("crap", 5), seq(5,1)), "fly"))
               )

#Now I want to get only the "normal" condition flies. It looks like most are in t24
with(mp.clean, table(line, condition))

mp.normal <- filter(mp.clean, condition == "t24")

#some of these are a little low. But can use for the whole G estimation still I think?
with(mp.normal, table(line, sex))

#Adding population 
mp.pop_dummy <- substr(mp.normal$line, 1,1)
mp.normal$population <- ifelse(mp.pop_dummy == "z", "Zambia", "Ethiopia")
mp.normal$population <- factor(mp.normal$population)

with(mp.normal, table(population, line))

#Is this enough?
#14 zi lines
length(levels(as.factor(filter(mp.normal, population == "Zambia")$line)))
#16 zi lines
length(levels(as.factor(filter(mp.normal, population == "Ethiopia")$line)))

kp.clean <- (kp.wings %>%
               separate(File, into = c(paste(rep("crap", 4), seq(4,1)), "name"), 
                        sep = "\\\\") %>%
               separate(name, into = c("perp", "line", "slide", "sex", "fly"), 
                        sep= "_") %>% 
               dplyr::select(!c(paste(rep("crap", 4), seq(4,1)), "fly", "slide"))
)

#I need to clean up my flies so that the tags match MP (I have zi and ef, she has z and e)
kp.clean$line <- sub("zi", "z", kp.clean$line)
kp.clean$line <- sub("ef", "e", kp.clean$line)
kp.clean$line 

#Adding population 
kp.pop_dummy <- substr(kp.clean$line, 1,1)
kp.clean$population <- ifelse(kp.pop_dummy == "z", "Zambia", "Ethiopia")
kp.clean$population <- factor(kp.clean$population)


#####Before I do anything else, I should check that these look like wings. 

WingBlur(kp.clean[sample(nrow(kp.clean), 50) ,15:110])
#This one has way more variation, but is way more lines. 
WingBlur(mp.normal[sample(nrow(mp.normal), 50) ,17:112])

#########################Finding MP Outliers to remove###########################
#Doing this in geomorph to look at outliers 
cord <- as.matrix(mp.normal[,17:112])
shape <- arrayspecs(cord, 48, 2)
mp.geo <- geomorph.data.frame(shape = shape,
                            CS = mp.normal$CS,
                            line = mp.normal$line,
                            condition = mp.normal$condition,
                            sex = mp.normal$sex)

#This gives the location of every point in this plot in the data (in order)
#What is a good outlier cuttoff
outliers <- unique(plotOutliers(mp.geo$shape)[1:100])
#why is this doubled now? but it was 100 observations still 
outliers

m <- mshape(mp.geo$shape)

#For now I will take the top 5 outliers. For some reason this is doubled. 
plotRefToTarget(m, mp.geo$shape[,,156], method = "points", links = wing.links)


###########################
#removing ef217 because we don't know what it is and any lines that have <10 observations in any line and those top 100 outliers from geomorph. 
with(mp.normal, table(sex, line))

problem.lines <- c("e217", "e117", "e54", "ef119", "z186", "z360", "z366", "e112")


mp.final <- (mp.normal %>%
               slice(-outliers) %>%
               filter(line %!in% problem.lines)
               
             ) 

#now I have to add e112 to the lines above (I did this, now everything has >10 individuals)
with(mp.final, table(sex, line))
with(mp.final, table(population, line)) #11 ef lines, 11 zi lines

#Now to add everything together. Need to make sure everything is the same

#"condition" "rep" 
names(mp.final)[(!(names(mp.final) %in% names(kp.clean)))]
#none
names(kp.clean)[(!(names(kp.clean) %in% names(mp.final)))]

#dropping condition and rep 
mp.final2 <- dplyr::select(mp.final, -c("condition", "rep"))

#Now it all matches. 
names(mp.final2)[(!(names(mp.final2) %in% names(kp.clean)))]

#Sticking the wings together 
wings <- rbind(kp.clean, mp.final2)

str(wings)
#zi251 is shared. 
with(wings, table(Perp, line))

with(wings, table(sex, line))

#shared line is bigger in MP data. But this is totally diffrent food and not suprising. 
(wings %>%
    group_by(line, Perp) %>%
    summarise(mean(CS)) %>%
    filter(line == "z251"))

#This distribution doesn't look totally crazy (the ef are still bigger overall ect) 
ggplot(wings, aes(x = line, y = CS, col = Perp)) +
  geom_boxplot()

#Ok now I think I am good to go with the projection 

wings$logCS <- log2(wings$CS)
####################################################################################################
#going to consider the sexes speperaltey like I have been 
wings.f <- filter(wings, sex == "f")
wings.m <- filter(wings, sex == "m")


#I want to actually do this with line means all the way down (helps to wash out environmental effects)

female.linemeans <- (wings.f %>%
                      mutate(perp.line = paste(perp, line, sep = "_")) %>%
                      group_by(perp.line) %>%
                      summarise_at(names(.)[c(15:111, 113)],
                                   .funs = c(mean="mean")) %>% 
                        separate(perp.line, into = c("perp", "line"), sep = "_")
)

female.linemeans.pop_dummy <- substr(female.linemeans$line, 1,1)
female.linemeans$population <- ifelse(female.linemeans.pop_dummy == "z", "Zambia", "Ethiopia")
female.linemeans$population <- factor(female.linemeans$population)

str(female.linemeans)

#Removing CS and spliner effect 

mod.perp.size.linemeans <-  lm(as.matrix(female.linemeans[,3:98]) ~ perp + logCS_mean,
                       data = female.linemeans)

f.linemeans.perpsize.resid <- data.frame(female.linemeans, mod.perp.size.linemeans$residuals)

#subsetting for mine and totals
ef.f.linemeans <- f.linemeans.perpsize.resid[f.linemeans.perpsize.resid$population == "Ethiopia",]

ef.f.linemeans.kp <- f.linemeans.perpsize.resid[f.linemeans.perpsize.resid$perp == "kp" &
                                                  f.linemeans.perpsize.resid$population == "Ethiopia",]

zi.f.linemeans <- f.linemeans.perpsize.resid[f.linemeans.perpsize.resid$population == "Zambia",]

zi.f.linemeans.kp <- f.linemeans.perpsize.resid[f.linemeans.perpsize.resid$perp == "kp" &
                                                  f.linemeans.perpsize.resid$population == "Zambia",]
#getting the selection vectors
ef.f.linemeans.feans <- colMeans(ef.f.linemeans[,102:197])
ef.f.linemeans.kp.feans <- colMeans(ef.f.linemeans.kp[,102:197])

zi.f.linemeans.feans <- colMeans(zi.f.linemeans[,102:197])
zi.f.linemeans.kp.feans <- colMeans(zi.f.linemeans.kp[,102:197])

all.diff.f <- ef.f.linemeans.feans - zi.f.linemeans.feans
kp.diff.f <- ef.f.linemeans.kp.feans - zi.f.linemeans.kp.feans

cor(all.diff.f, kp.diff.f) #well this is really low -0.1961824. same as before I did line means. 


######What if I leave size in there?#######

perp.fod.f.linemeans <-  lm(as.matrix(female.linemeans[,3:98]) ~ perp,
                            data = female.linemeans)


linemeans.perp.resid.f <- data.frame(female.linemeans, perp.fod.f.linemeans$residuals)

#subsetting for mine and totals
ef.f.withSize <- linemeans.perp.resid.f[linemeans.perp.resid.f$population == "Ethiopia",]

ef.f.kp.withSize <- linemeans.perp.resid.f[linemeans.perp.resid.f$perp == "kp" &
                                             linemeans.perp.resid.f$population == "Ethiopia",]


zi.f.withSize <- linemeans.perp.resid.f[linemeans.perp.resid.f$population == "Zambia",]

zi.f.kp.withSize <- linemeans.perp.resid.f[linemeans.perp.resid.f$perp == "kp" &
                                             linemeans.perp.resid.f$population == "Zambia",]


#getting the selection vectors
ef.all.f.feans.withSize <- colMeans(ef.f.withSize[,102:197])
ef.kp.f.feans.withSize <- colMeans(ef.f.kp.withSize[,102:197])

zi.all.f.feans.withSize <- colMeans(zi.f.withSize[,102:197])
zi.kp.f.feans.withSize <- colMeans(zi.f.kp.withSize[,102:197])

all.diff.f.withSize <- ef.all.f.feans.withSize - zi.all.f.feans.withSize
kp.diff.f.withSize <- ef.kp.f.feans.withSize - zi.kp.f.feans.withSize

cor(all.diff.f.withSize, kp.diff.f.withSize) # 0.4947965


#########Using just the landmarks ###################

#subsetting for mine and totals
ef.f.lm <- female.linemeans[female.linemeans$population == "Ethiopia",]

ef.f.kp.lm <- female.linemeans[female.linemeans$perp == "kp" &
                                 female.linemeans$population == "Ethiopia",]


zi.f.lm <- female.linemeans[female.linemeans$population == "Zambia",]

zi.f.kp.lm <- female.linemeans[female.linemeans$perp == "kp" &
                                 female.linemeans$population == "Zambia",]


#getting the selection vectors
ef.all.f.feans.lm <- colMeans(ef.f.lm[,3:98])
ef.kp.f.feans.lm <- colMeans(ef.f.kp.lm[,3:98])

zi.all.f.feans.lm <- colMeans(zi.f.lm[,3:98])
zi.kp.f.feans.lm <- colMeans(zi.f.kp.lm[,3:98])

all.diff.f.lm <- ef.all.f.feans.lm - zi.all.f.feans.lm
kp.diff.f.lm <- ef.kp.f.feans.lm - zi.kp.f.feans.lm

#almost exactly the same if I leave the spliner effect in and just use the landmarks
cor(all.diff.f.lm, kp.diff.f.lm)  # 0.4877295

#Also could  with unique allometries (that might be what we are seeing messing up the first one)
#Actually just going to allow for this to have a unique componet for each line.  
unique.allo.fod <- lm(as.matrix(female.linemeans[,3:98]) ~ logCS_mean + logCS_mean:line + perp,
                      data = female.linemeans)


unique.allo.linemeans <- data.frame(female.linemeans, unique.allo.fod$residuals)

#subsetting for mine and totals
ef.f.linemeans.uniqueallo <- unique.allo.linemeans[unique.allo.linemeans$population == "Ethiopia",]

ef.f.linemeans.uniqueallo.kp <- unique.allo.linemeans[unique.allo.linemeans$perp == "kp" &
                                                        unique.allo.linemeans$population == "Ethiopia",]

zi.f.linemeans.uniqueallo <- unique.allo.linemeans[unique.allo.linemeans$population == "Zambia",]

zi.f.linemeans.uniqueallo.kp <- unique.allo.linemeans[unique.allo.linemeans$perp == "kp" &
                                                        unique.allo.linemeans$population == "Zambia",]
#getting the selection vectors
zi.f.linemeans.uniqueallo.feans <- colMeans(zi.f.linemeans.uniqueallo[,102:197])
zi.f.linemeans.uniqueallo.kp.feans <- colMeans(zi.f.linemeans.uniqueallo.kp[,102:197])

zi.f.linemeans.uniqueallo.feans <- colMeans(zi.f.linemeans.uniqueallo[,102:197])
zi.f.linemeans.uniqueallo.kp.feans <- colMeans(zi.f.linemeans.uniqueallo.kp[,102:197])

uniqueallo.all.diff.f <- zi.f.linemeans.uniqueallo.feans - zi.f.linemeans.uniqueallo.feans
uniqueallo.kp.diff.f <- zi.f.linemeans.uniqueallo.kp.feans - zi.f.linemeans.uniqueallo.kp.feans

cor(all.diff.f, kp.diff.f) # -0.1961824 this did not help. 

#just using the resid going forward.  
#####As a summary: spliner effect doesn't mess with the cor here. I think that is good since there is a also a block effect wrapped up in this (my wings are from the Pool lab and MP are from her experiment)
#Modeling out logCS has a big effect on the cor. To me this means that there is a really different non-allometric shape change captured in my wings.

#########################Projection########################

#So now I want to create a quick and dirty G matrix using the line means for all of these. I am going to do this both with and without allometric shape?


af.linemeans.lm.pc <- prcomp(female.linemeans[,3:98])

summary(af.linemeans.lm.pc)

#taking them all because why not.
linemeans.f.lm.pc <- data.frame(female.linemeans, af.linemeans.lm.pc$x[,1:28])

linemeans.f.lm.pc$kp.score.lm <- projFunction(as.matrix(linemeans.f.lm.pc[,3:98]),
                                              (as.matrix(kp.diff.f.lm)))


linemeans.f.lm.pc$all.score.lm <- projFunction(as.matrix(linemeans.f.lm.pc[,3:98]),
                                               (as.matrix(all.diff.f.lm)))


linemeans.f.lm.pc$kp.score.noAllo <- projFunction(as.matrix(linemeans.f.lm.pc[,3:98]),
                                                  (as.matrix(kp.diff.f)))


linemeans.f.lm.pc$all.score.noAllo <- projFunction(as.matrix(linemeans.f.lm.pc[,3:98]),
                                                   (as.matrix(all.diff.f)))


#Yeah. My lines/selection vector is doing something different that isn't awesome for capturing the shape change between populations.
png("../Figures/females_af_cor_allo.png", width = 1500, height = 1500, units = "px", res = 300)
pairs(linemeans.f.lm.pc[,c(132:133, 102:105)], lower.panel = panel.cor)
dev.off()


#Now I want to plot just the Zi 

zi.f.lm.pc <- filter(linemeans.f.lm.pc, population == "Zambia")

zi.f.lm.pca <-  data.frame(zi.f.lm.pc[,c(1:101, 130:133)], prcomp(zi.f.lm.pc[,3:98])$x[,1:10])


png("../Figures/females_ziOnly_cor_allo.png", width = 1500, height = 1500, units = "px", res = 300)
pairs(zi.f.lm.pca[,104:109], lower.panel = panel.cor)
dev.off()

####################################Males########################################################

male.linemeans <- (wings.m %>%
                     mutate(perp.line = paste(perp, line, sep = "_")) %>%
                     group_by(perp.line) %>%
                     summarise_at(names(.)[c(15:111, 113)],
                                  .funs = c(mean="mean")) %>% 
                     separate(perp.line, into = c("perp", "line"), sep = "_")
)

male.linemeans.pop_dummy <- substr(male.linemeans$line, 1,1)
male.linemeans$population <- ifelse(male.linemeans.pop_dummy == "z", "Zambia", "Ethiopia")
male.linemeans$population <- factor(male.linemeans$population)

str(male.linemeans)

#Removing CS and spliner effect 

mod.perp.size.linemeans.m <-  lm(as.matrix(male.linemeans[,3:98]) ~ perp + logCS_mean,
                                 data = male.linemeans)

m.linemeans.perpsize.resid <- data.frame(male.linemeans, mod.perp.size.linemeans.m$residuals)

#subsetting for mine and totals
ef.m.linemeans <- m.linemeans.perpsize.resid[m.linemeans.perpsize.resid$population == "Ethiopia",]

ef.m.linemeans.kp <- m.linemeans.perpsize.resid[m.linemeans.perpsize.resid$perp == "kp" &
                                                  m.linemeans.perpsize.resid$population == "Ethiopia",]

zi.m.linemeans <- m.linemeans.perpsize.resid[m.linemeans.perpsize.resid$population == "Zambia",]

zi.m.linemeans.kp <- m.linemeans.perpsize.resid[m.linemeans.perpsize.resid$perp == "kp" &
                                                  m.linemeans.perpsize.resid$population == "Zambia",]
#getting the selection vectors
ef.m.linemeans.means <- colMeans(ef.m.linemeans[,102:197])
ef.m.linemeans.kp.means <- colMeans(ef.m.linemeans.kp[,102:197])

zi.m.linemeans.means <- colMeans(zi.m.linemeans[,102:197])
zi.m.linemeans.kp.means <- colMeans(zi.m.linemeans.kp[,102:197])

all.diff.m <- ef.m.linemeans.means - zi.m.linemeans.means
kp.diff.m <- ef.m.linemeans.kp.means - zi.m.linemeans.kp.means

cor(all.diff.m, kp.diff.m) #more similar to females than before. -0.2687458


######What if I leave size in there?#######

perp.mod.m.linemeans <-  lm(as.matrix(male.linemeans[,3:98]) ~ perp,
                            data = male.linemeans)


linemeans.perp.resid.m <- data.frame(male.linemeans, perp.mod.m.linemeans$residuals)

#subsetting for mine and totals
ef.m.withSize <- linemeans.perp.resid.m[linemeans.perp.resid.m$population == "Ethiopia",]

ef.m.kp.withSize <- linemeans.perp.resid.m[linemeans.perp.resid.m$perp == "kp" &
                                             linemeans.perp.resid.m$population == "Ethiopia",]


zi.m.withSize <- linemeans.perp.resid.m[linemeans.perp.resid.m$population == "Zambia",]

zi.m.kp.withSize <- linemeans.perp.resid.m[linemeans.perp.resid.m$perp == "kp" &
                                             linemeans.perp.resid.m$population == "Zambia",]


#getting the selection vectors
ef.all.m.means.withSize <- colMeans(ef.m.withSize[,102:197])
ef.kp.m.means.withSize <- colMeans(ef.m.kp.withSize[,102:197])

zi.all.m.means.withSize <- colMeans(zi.m.withSize[,102:197])
zi.kp.m.means.withSize <- colMeans(zi.m.kp.withSize[,102:197])

all.diff.m.withSize <- ef.all.m.means.withSize - zi.all.m.means.withSize
kp.diff.m.withSize <- ef.kp.m.means.withSize - zi.kp.m.means.withSize

cor(all.diff.m.withSize, kp.diff.m.withSize) # 0.5281066


#########Using just the landmarks ###################

#subsetting for mine and totals
ef.m.lm <- male.linemeans[male.linemeans$population == "Ethiopia",]

ef.m.kp.lm <- male.linemeans[male.linemeans$perp == "kp" &
                               male.linemeans$population == "Ethiopia",]


zi.m.lm <- male.linemeans[male.linemeans$population == "Zambia",]

zi.m.kp.lm <- male.linemeans[male.linemeans$perp == "kp" &
                               male.linemeans$population == "Zambia",]


#getting the selection vectors
ef.all.m.means.lm <- colMeans(ef.m.lm[,3:98])
ef.kp.m.means.lm <- colMeans(ef.m.kp.lm[,3:98])

zi.all.m.means.lm <- colMeans(zi.m.lm[,3:98])
zi.kp.m.means.lm <- colMeans(zi.m.kp.lm[,3:98])

all.diff.m.lm <- ef.all.m.means.lm - zi.all.m.means.lm
kp.diff.m.lm <- ef.kp.m.means.lm - zi.kp.m.means.lm

#almost exactly the same if I leave the spliner effect in and just use the landmarks
cor(all.diff.m.lm, kp.diff.m.lm)  #  0.4583184

#Also could  with unique allometries (that might be what we are seeing messing up the first one)
#Actually just going to allow for this to have a unique componet for each line.  
unique.allo.mod.m <- lm(as.matrix(male.linemeans[,3:98]) ~ logCS_mean + logCS_mean:line + perp,
                        data = male.linemeans)


unique.allo.linemeans.m <- data.frame(male.linemeans, unique.allo.mod.m$residuals)

#subsetting for mine and totals
ef.m.linemeans.uniqueallo <- unique.allo.linemeans.m[unique.allo.linemeans$population == "Ethiopia",]

ef.m.linemeans.uniqueallo.kp <- unique.allo.linemeans[unique.allo.linemeans$perp == "kp" &
                                                        unique.allo.linemeans$population == "Ethiopia",]

zi.m.linemeans.uniqueallo <- unique.allo.linemeans[unique.allo.linemeans$population == "Zambia",]

zi.m.linemeans.uniqueallo.kp <- unique.allo.linemeans[unique.allo.linemeans$perp == "kp" &
                                                        unique.allo.linemeans$population == "Zambia",]
#getting the selection vectors
zi.m.linemeans.uniqueallo.means <- colMeans(zi.m.linemeans.uniqueallo[,102:197])
zi.m.linemeans.uniqueallo.kp.means <- colMeans(zi.m.linemeans.uniqueallo.kp[,102:197])

zi.m.linemeans.uniqueallo.means <- colMeans(zi.m.linemeans.uniqueallo[,102:197])
zi.m.linemeans.uniqueallo.kp.means <- colMeans(zi.m.linemeans.uniqueallo.kp[,102:197])

uniqueallo.all.diff.m <- zi.m.linemeans.uniqueallo.means - zi.m.linemeans.uniqueallo.means
uniqueallo.kp.diff.m <- zi.m.linemeans.uniqueallo.kp.means - zi.m.linemeans.uniqueallo.kp.means

cor(all.diff.m, kp.diff.m) # --0.2657031 this did not help. 

#just using the resid going forward.  
#####As a summary: spliner effect doesn't mess with the cor here. I think that is good since there is a also a block effect wrapped up in this (my wings are from the Pool lab and MP are from her experiment)
#Modeling out logCS has a big effect on the cor. To me this means that there is a really different non-allometric shape change captured in my wings.

#########################Projection########################

#So now I want to create a quick and dirty G matrix using the line means for all of these. I am going to do this both with and without allometric shape?


af.linemeans.lm.pc <- prcomp(male.linemeans[,3:98])

summary(af.linemeans.lm.pc)

#taking them all because why not.
linemeans.m.lm.pc <- data.frame(female.linemeans, af.linemeans.lm.pc$x[,1:28])

linemeans.m.lm.pc$kp.score.lm <- projFunction(as.matrix(linemeans.m.lm.pc[,3:98]),
                                              (as.matrix(kp.diff.m.lm)))


linemeans.m.lm.pc$all.score.lm <- projFunction(as.matrix(linemeans.m.lm.pc[,3:98]),
                                               (as.matrix(all.diff.m.lm)))


linemeans.m.lm.pc$kp.score <- projFunction(as.matrix(linemeans.m.lm.pc[,3:98]),
                                                  (as.matrix(kp.diff.m)))


linemeans.m.lm.pc$all.score <- projFunction(as.matrix(linemeans.m.lm.pc[,3:98]),
                                                   (as.matrix(all.diff.m)))


#Yeah. My lines/selection vector is doing something different that isn't awesome for capturing the shape change between populations.
png("../Figures/males_af_cor_allo.png", width = 1500, height = 1500, units = "px", res = 300)
pairs(linemeans.m.lm.pc[,c(132:133, 102:105)], lower.panel = panel.cor)
dev.off()



#looking at just the zambians. 





