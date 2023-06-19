#I wan to ask if I can separate the big cell size and more cells in the big wings. 
library(tidyverse)
library(lme4)
library(effects)
library(emmeans)
library(car)

#First I want to recap Maria's results with my own from the parents. Also using some old code ID and I wrote for an undergrad. 

###############################################################
#read in clean and checks
###############################################################

#for the Tricomes, NA are coded as 0s. 
parents <- read.csv("../Data/pureLines_cellCount.csv")

parents.shape <- read.csv("../Data/pureWings_15lmk_superimposed_withFly.csv")

#Cleaning up the tricome stuff first. following Maria's stuff (although I think this filter is too harsh)
#I dont know why that is having trouble with the NAs. 
str(parents)

sum(is.na(parents$TricomeCount))

parents$TricomeCount[parents$TricomeCount == 0] <- NA
#good. 
sum(is.na(parents$TricomeCount))

#I also need to add a factor for which wings are a 2x and 4x because there are diff #mm/px in each 
#This is a dictionary I made with the sex and line of populations. 
wingSizeDict <- read.delim("pureWingImageDict.txt")



#breaking the ID col apart. Will make the fly ID 
parents <- (parents %>%
      separate(ID, into = c("perp", "line", "slide", "sex", "fly"), sep = "_") %>% 
        separate(fly, into = c("fly", "crap")))

#This matches the line and sex and adds the correct mag
parents <- left_join(parents, wingSizeDict, by = c("line", "sex"))

parents$flyID <- paste0(parents$line, parents$sex, parents$fly)

#need to make these match. 
str(parents.shape)
parents.shape$flyID <- paste0(parents.shape$line, parents.shape$sex, parents.shape$fly)

#Want to look at these. 
parents[which(parents$TricomeCount > 120),]
parents[which(parents$TricomeCount < 50),] #everything is really close to this. I don't think I am going to cut any out right away. but this is a note for myself. 

#not terrible 
hist(as.numeric(parents$TricomeCount))

#warning is for the NAs for counts. 
#this looks pretty alright to me. 
ggplot(parents, aes(y = TricomeCount, x = line, colour = sex)) +
  geom_boxplot()

#sqishing stuff together to make it so ugly but the variation is there
ggplot(parents, aes(y = TricomeCount, x = as.factor(Area), colour = line)) +
  geom_violin()

parents.shape$slide <- as.character(parents.shape$slide)

parents.shape$fly <- as.character(parents.shape$fly)

###############################################################
#CS to mm2 conversions. 
###############################################################
#I measured about 10 wings from the 4x and 2x set to get these. Used the lab macro without the crop and with the mm2 output with the correct px/mm set **needs to be changed for 2x from the one on dropbox**

#all the crap is extra stuff I should just delete (observation number and stats abotu the colour of the area mesured. its all black)
macroSizeHeader <- c("crap", "image", "size", "crap2", "crap3", "crap4")

macroSize2x <- read.csv("../Data/sizeMacroConversion2x.csv", header = F)
names(macroSize2x) <- macroSizeHeader
macroSize4x <- read.csv("../Data/sizeMacroConversion4x.csv", header = F)
names(macroSize4x) <- macroSizeHeader

#Splitting the image name apart. 
macroSize2x <- (macroSize2x %>%
              separate(image, into = c("perp", "line", "slide", "sex", "fly"), sep = "_") %>% 
              separate(fly, into = c("fly", "crap")))

hist(macroSize2x$size)

macroSize4x <- (macroSize4x %>%
                  separate(image, into = c("perp", "line", "slide", "sex", "fly"), sep = "_") %>% 
                  separate(fly, into = c("fly", "crap")))
#these numbers look about right to me. 
hist(macroSize4x$size)

#Now I want to match these with the right wings from the shape data do get the regression. 

#I need to do more because I grabbed dropped ones. 
sizeCS2x <- left_join(macroSize2x, parents.shape, by = c("line", "slide", "sex", "fly"))
#11 is good enough 
sizeCS2x <- sizeCS2x[complete.cases(sizeCS2x),]

sizeCS4x <- left_join(macroSize4x, parents.shape, by = c("line", "slide", "sex", "fly"))
#20 is also fine. 
sizeCS4x <- sizeCS4x[complete.cases(sizeCS4x),]


all.size <- rbind(sizeCS2x, sizeCS4x)
plot(all.size$CS, all.size$size)
cor(all.size$CS, all.size$size)
#I also want to mean centre size so I don't have to worry about the intercept. 
all.size$CS.c <- all.size$CS - mean(all.size$CS)
all.size$size.c <- sqrt(all.size$size) - mean(sqrt(all.size$size))



all.size.mod  <- lm(CS ~ size.c, data = all.size)
summary(all.size.mod )
#slope: 2.073490

#I also want to mean centre size so I don't have to worry about the intercept. 
# sizeCS4x$CS.c <- sizeCS4x$CS - mean(sizeCS4x$CS)
# sizeCS2x$CS.c <- sizeCS2x$CS - mean(sizeCS2x$CS)
# 
# sizeCS4x$size.c <- sqrt(sizeCS4x$size) - mean(sqrt(sizeCS4x$size))
# sizeCS2x$size.c <- sqrt(sizeCS2x$size) - mean(sqrt(sizeCS2x$size))

#take square root. 
# #mean center 
# mod2x <- lm(CS ~ size.c, data = sizeCS2x)
# summary(mod2x)
# #slope: 2.07081
# 
# mod4x <- lm(CS ~ size.c, data = sizeCS4x)
# summary(mod4x)
# #slope: 2.291337 

#I think I can just multiply CS by the slope to get mm^2. They really will let anyone get a PhD 

# parents.shape <- left_join(parents.shape, wingSizeDict, by = c("line", "sex"))
# 
# conversion.mat <- data.frame(c("2x", "4x"), c(2.07081, 2.291337))
# names(conversion.mat) <- c("mag", "converstionFact")
# 
# parents.shape <- left_join(parents.shape, conversion.mat, by = c("mag"))

parents.shape$wingSize <- parents.shape$CS*2.073490

hist(parents.shape$CS)
#why so bimodal? is this just pulling out the bigest and smallest more because of the different multiplication? Making the big so much bigger? I will have to ask Ian this question. 
hist(parents.shape$wingSize)

# parents.shape$crap <- parents.shape$CS*2.24383
ggplot(parents.shape, aes(x = line, y = wingSize, col = sex)) +
  geom_violin()

###############################################################
#Now to actually be able to do analysis. 
###############################################################
#joining the datasets together. 
#I need to print the parent shape data with the slide identifier. 
#this looks right?

#This looks like I loose a lot but I did tricomes for every wing and I did not do shape for every wing. this is done in a way so I only get macthing cases. 
wings <- inner_join(parents, parents.shape, by = c("slide", "line", "sex", "flyID"))

#need to rename my Area to wingRegion 
wings$wingRegion <- wings$Area

wings <- wings[,2:48]
#about 30 wings/sex/region. Most of these are compleate wings because I only had 78 NAs in the entire dataset (and I think a lot of those would be impossible to spline wings) 

#almost identical. Seems like with dropping some, this still turns out fine. 
(parents.shape %>% group_by(sex, line) %>% summarise(meanSize = mean(CS), sdSize = sd(CS)))

(wings %>% group_by(sex, line) %>% summarise(meanSize = mean(CS), sdSize = sd(CS)))

#and with actual size of wings. in mm^2
#This generally makes sense with what I think 
(wings %>% group_by(sex, line) %>% summarise(meanSize = mean(wingSize), sdSize = sd(wingSize)))

ggplot(wings, aes(x = line, y = CS, col = sex)) + 
  geom_violin()

#This is still wrong. but less wrong. I think that considering them together might be the way?  
ggplot(wings, aes(x = line, y = wingSize, col = sex)) + 
  geom_violin()

###############################################################
#average density across all regions. Obviously imperfect assumption
#Also averaging shape so I have it?
###############################################################
indiv_cells <-  aggregate(wings[,c(7, 10:40, 47)], 
                          c(wings["sex"],
                            wings["line"],
                            wings["slide"],
                            wings["fly.x"], 
                            wings["mag.x"], 
                            wings["population"]), 
                          mean, na.rm=TRUE )


ggplot(indiv_cells, aes(x = line, y = TricomeCount, col = sex)) + 
  geom_violin()

ggplot(indiv_cells, aes(x = CS, y = TricomeCount, col = sex, shape = line)) + 
  geom_point() + 
  geom_smooth(method = "lm")

indiv_cells$pop <- ifelse(grepl("^Z", indiv_cells$line), "LA", "HA")

#I did this with just the lines and thought that this might be easier to understand (plus I really care about the alt effect, not the specific line contrasts)
tricome.lm <- lmer(TricomeCount ~ CS*sex*pop + (1|pop:line), data = indiv_cells)

#Why does wing size have such a small effect? What did I do wrong?????
#Is it because a lot of the change is in the interactions (as in the effects depend a lot on the population)
#I think it is the weird fucking zambian males that suck. 
summary(tricome.lm)
Anova(tricome.lm)

################################################################
#yeah. not crazy. I like when things work so well. 
with(wings, table(sex, line, wingRegion))
with(indiv_cells, table(sex, line))

#Right on, which it should be. 
(wings %>% group_by(sex, line)  %>% summarise(meanSize = mean(CS), sdSize = sd(CS)))
(indiv_cells %>% group_by(sex, line) %>% summarise(meanSize = mean(CS), sdSize = sd(CS)))

(wings %>% group_by(sex, line)  %>% summarise(meanSize = mean(wingSize), sdSize = sd(wingSize)))

(indiv_cells %>% group_by(sex, line) %>% summarise(meanSize = mean(wingSize), sdSize = sd(wingSize)))

#lower density in EF than ZI. Overall this is good.  
(indiv_cells %>% group_by(sex, line) %>% summarise(meanDensity = mean(TricomeCount), sdDensity = sd(TricomeCount)))

#Now I need to get the number of tricomes/mm^2
#According to ID in the explainer he posted on basecamp, this is actually a 150x150 px square. 

#I give up on trying to understand the math and I'm just going to use ID function here. 
# trichomes_in_wings <- function(trichome_count, wing_area_mm, px_mm = 1860, box_px = 150 ) {
#   
#   sq_mm_in_px <- 1/(px_mm^2)
#   area_in_px <- box_px^2
#   area_mm_sq <- area_in_px*sq_mm_in_px
#   trichome_density_mm <- trichome_count/area_mm_sq
#   trichome_count <- wing_area_mm*trichome_density_mm 
#   
#   return(trichome_count)
# }

#calclulated using 1/(mm/px) measurement off the github for the 413 scope 
#4x: px_mm = 1860
#2x: px_mm = 925

#sq_mm_in_px = 1/px_mm^2
#4x= 2.890508e-07
#2x= 1.168736e-06

#area_in_px =  box_px^2 = 150^2 = 22500

#area_mm_sq <- area_in_px*sq_mm_in_px
#4x = 2.890508e-07*22500 = 0.006503643
#2x = 1.168736e-06*22500 = 0.02629656

#tricome denisty mm = trichome_count/area_mm_sq

#total_trichome_count <- wing_area_mm*trichome_density_mm


wings$wingSize_mm2 <- wings$wingSize^2

px.convMat <- data.frame(c("2x", "4x"), c(0.02629656, 0.006503643))
names(px.convMat) <- c("mag", "pxConversion")
  
indiv_cells$mag <- indiv_cells$mag.x
indiv_cells <- left_join(indiv_cells, px.convMat, by = "mag")

indiv_cells$tricome.dens_mm2 <- indiv_cells$TricomeCount/indiv_cells$pxConversion


#multiply the number of cells/mm^2 by wing area to get total number of cells in the wing 
indiv_cells$totalCells <- indiv_cells$wingSize * indiv_cells$tricome.dens_mm2

#Nope something is super wrong with these. The ones at 2x are fucked. 
#Must return to this problem at a later date. 
ggplot(indiv_cells, aes(x = line, y = totalCells, col = sex)) + 
  geom_violin()

#divide the total area of wing by number of cells to get average cell size. 
indiv_cells$cellsize <- indiv_cells$wingSize/indiv_cells$totalCells

#maybe some suspect stuff here. Need to look at the wings above.
#I love this for me. 
ggplot(indiv_cells, aes(x = line, y = cellsize, col = sex)) + 
  geom_violin()





