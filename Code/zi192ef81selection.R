#This is the selection for zi192xef81 wings for sequencing, following the same protocol as zi192xef96 cross submitted summer 2020

library(tidyverse)
library(MASS)

######### Loading functions #########
source('~/Dropbox/DworkinLabSharedMaterial/scripts/WRP_FUNCTIONS.R', chdir = TRUE)

source('~/Dropbox/DworkinLabSharedMaterial/scripts/WINGPLOTSOURCE.R', chdir = TRUE)

####################################

# #read in dat
# wings_raw <- read_tsv("../Data/zi192ef81allwings.dat")
# 
# #need to seperate out PL and cross wings
# crosswings <- filter(wings_raw, grepl("plate", wings_raw$CPFile))
# PLwings <- filter(wings_raw, !grepl("plate", wings_raw$CPFile))
# #602+2638 = 3240
# 
# Crap_Col <- c("File", "O1x", "O1y", "O2x", "O2y", "Date", "Time", "Tags", "Sex", "Perp")
# More_crap <- c("Wings", "KP", "AfPop", "PL", "folder", "folder2")
# 
# 
# 
# PL_fixed <- (PLwings
#                 %>% separate("CPFile", into= c("perp", "line", "slide", "Real_Sex", "Fly"), sep = "_")
#                 %>% select(-Crap_Col)
#                 %>% unite(Fly_ID, slide, Fly)
#                )
# 
# PL_fixed$pop <- ifelse(grepl("ZI", PL_fixed$line), "LA", "HA")
# levels(as.factor(PL_fixed$pop))
# 
# cross_fixed <- (crosswings
#                 %>%separate("CPFile", into= c("perp", "line", "colection", "slide", "Real_Sex", "Fly"), sep = "_")
#                 %>% select(-Crap_Col)
#                 %>% select(-colection)
#                 %>% unite(Fly_ID, slide, Fly)
#                 )
# cross_fixed$pop <- "cross"
# 
# 
# wings <- rbind(cross_fixed, PL_fixed)
# colnames(wings)[4] <- "sex"
# 
# wings$sex <- toupper(wings$sex)
# 
# wings[,c(1:4, 103)] <- lapply(wings[,c(1:4, 103)], as.factor)
# 
# save(wings, file = "../Data/zi192ef81clean.rda")

load("../Data/zi192ef81clean.rda")

wings$logCS <- log(wings$CS)

#Looks like wings. Nothing weird. 
WingBlur(wings[sample(nrow(wings), 50) ,6:101])
#WingBlur(wings[,6:101])

#calculating the PCs

PCs <- prcomp(wings[, c(6:101)])
wings_combined <- data.frame(wings, PCs$x[,1:57]) 

#ruh roh looks like something is wrong with the biggest wings. dammit. This is the same math I used with Sarah's cvl wings? 
#Went back and checked. They are just big
ggplot(wings_combined, aes(x = logCS, y = PC1, col = pop)) + geom_point()

#Calculating PL means with PCs

wingsHAF <- filter(wings_combined, sex == "F" & pop == "HA")
wingsLAF <- filter(wings_combined, sex == "F" & pop == "LA")
wingsHAM <- filter(wings_combined, sex == "M" & pop == "HA")
wingsLAM <- filter(wings_combined, sex == "M" & pop == "LA")

HAmeanF <- colMeans(wingsHAF[,105:161])
LAmeanF <- colMeans(wingsLAF[,105:161])
HAmeanM <- colMeans(wingsHAM[,105:161])
LAmeanM <- colMeans(wingsLAM[,105:161])

Fdiff <- HAmeanF - LAmeanF
Mdiff <- HAmeanM - LAmeanM
cor(Fdiff, Mdiff) #0.91

#now with landmarks 

lmwingsHAF <- filter(wings, sex == "F" & pop == "HA")
lmwingsLAF <- filter(wings, sex == "F" & pop == "LA")
lmwingsHAM <- filter(wings, sex == "M" & pop == "HA")
lmwingsLAM <- filter(wings, sex == "M" & pop == "LA")

lmHAmeanF <- colMeans(lmwingsHAF[,6:101])
lmLAmeanF <- colMeans(lmwingsLAF[,6:101])
lmHAmeanM <- colMeans(lmwingsHAM[,6:101])
lmLAmeanM <- colMeans(lmwingsLAM[,6:101])

lmFdiff <- lmHAmeanF - lmLAmeanF
lmMdiff <- lmHAmeanM - lmLAmeanM
#same answer. Of course. 
cor(lmFdiff, lmMdiff) #0.91

#In the other popualtion, I considered males and females seperatley, I want to do that again. But first I want to double check that these are sexed accuratley 

# #regressing size on shape 
# 
# res <- lm(as.matrix(wings[,6:101]) ~ logCS, data = wings)$residuals
# 
# sizewings <- data.frame(wings, res)
# 
# #Using the purelines as a training set
# 
 train <- filter(wings, pop != "cross")
 train <- train[,c(4, 6:101)]
#
 linDiscrim <- lda(formula = sex ~ ., data = train, tol = 1e-04, CV = FALSE)
#
 train_table <- table(actual = train$sex, predicted = predict(linDiscrim, newdata = train)$class)
#
#Only 97% accutate? 
 (281+304)/nrow(train)
#
 testall <- filter(wings, pop == "cross")
 test <- testall[,c(4, 6:101)]
#
 test_table <- table(actual = test$sex, predicted = predict(linDiscrim, newdata = test)$class)
#
# #this seems...a little low still. 91%
 (1248+1178)/nrow(test)
 
 
 #I want to look at the pedicted vs real sex by indiv. 
 
 testall$pred <-  predict(linDiscrim, newdata = test)$class
check <- filter(testall, testall$sex != testall$pred)
checkFly <- check[,c("Fly_ID", "sex", "pred")]
 
# #going to include size because that is so tied to dimorphism.
# 
# train <- filter(wings, pop != "cross")
# train <- train[,c(4, 6:101, 104)]
# 
# linDiscrim <- lda(formula = sex ~ ., data = train, tol = 1e-04, CV = FALSE)
# 
# train_table <- table(actual = train$sex, predicted = predict(linDiscrim, newdata = train)$class)
# 
# #99% now. Seems better. Makes sense because size matters.  
# (291+306)/nrow(train)
# 
# test <- wings[,c(4, 6:101, 104)]
# 
# test_table <- table(actual = test$sex, predicted = predict(linDiscrim, newdata = test)$class)
# 
# #95% is good enough for me. Although not perfect. 
# (1516+1593)/nrow(test)


#for the previous popualtion, we chose to look at the males and females seperatley. So I will do that again. 

females <- filter(wings, sex == "F")
#Using the landmarks plus logCS
PCAf <- prcomp(females[,c(6:101, 104)])
#So much of this variance is tied up in size/allometric diffrences! 85%!
summary(PCAf)

PCsF <- PCAf$x[,1:56]
rotationF <- PCAf$rotation

femalesPC <- data.frame(females, PCsF)

#that is a very nice line. 
ggplot(femalesPC, aes(x = logCS, y = PC1, col = pop)) + geom_point()

#This is not a line. Way more clustered than previous population I looked at. This makes sence because there is so much more variance explained by PC1 here. 
ggplot(femalesPC, aes(x = PC2, y = PC3, col = pop)) + 
  geom_point() + 
  geom_text(aes(label=as.character(Fly_ID),hjust=0,vjust=0))

#plate25c_12 and plate 23g_07 look a little suspect here. Should double check splines.
#double checked. They are fine.Keeping this line so I know I saw those and took care of it. 
WingPlot(as.matrix(filter(wings, Fly_ID == "plate25c_12"))[,6:101])


#not going to use PC1 because it captures the allometric componet. (might try to do this again below with the landmarks and modeling out size to see what happens) 


HA_mean_F <- colMeans(femalesPC[femalesPC$pop == "HA", 106:160])
LA_mean_F <- colMeans(femalesPC[femalesPC$pop == "LA", 106:160])
diff_mean_F <- HA_mean_F - LA_mean_F

femalesPC$diff <-  as.matrix(femalesPC[,106:160]) %*% diff_mean_F

female_plot <- ggplot(femalesPC, aes(x = CS, y = diff, col = pop)) +
  geom_point(alpha = 0.6, size = 2 )

ggplot(femalesPC, aes(x = CS, y = diff, col = line)) +
  geom_point(alpha = 0.6, size = 2 )

#########Males#################

males <- filter(wings, sex == "M")
#Using the landmarks plus logCS
PCAm <- prcomp(males[,c(6:101, 104)])
#So much of this variance is tied up in size/allometric diffrences! 86%!
summary(PCAm)

PCsM <- PCAm$x[,1:56]
rotationM <- PCAm$rotation

malesPC <- data.frame(males, PCsM)

#that is a very nice line. 
ggplot(malesPC, aes(x = logCS, y = PC1, col = pop)) + geom_point()

#This is not a line. Way more clustered than previous population I looked at. This makes sence because there is so much more variance explained by PC1 here. 
ggplot(malesPC, aes(x = PC2, y = PC3, col = pop)) + 
  geom_point() + 
  geom_text(aes(label=as.character(Fly_ID),hjust=0,vjust=0))

#plate9h_07 and plate9b_08 look suspect.
#checked actual splines. They are fine. 
WingPlot(as.matrix(filter(wings, Fly_ID == "plate9b_08"))[,6:101])

#not going to use PC1 because it captures the allometric componet. (might try to do this again below with the landmarks and modeling out size to see what happens) 
HA_mean_M <- colMeans(malesPC[malesPC$pop == "HA", 106:160])
LA_mean_M <- colMeans(malesPC[malesPC$pop == "LA", 106:160])
diff_mean_M <- HA_mean_M - LA_mean_M

malesPC$diff <-  as.matrix(malesPC[,106:160]) %*% diff_mean_M

#That one male is GIANT. How much do I actually trust that? I should select 51 and sex them as I go?
#I am actually really worred about mis-sexed flies now. 
male_plot <- ggplot(malesPC, aes(x = CS, y = diff, col = pop)) +
  geom_point(alpha = 0.6, size = 2 )

ggplot(malesPC, aes(x = CS, y = diff, col = line)) +
  geom_point(alpha = 0.6, size = 2 )

###Now I want to comare those two vectors I selected on. 
#rotation matrix looks pretty similar for the first couple PCs then falls apart

cor(rotationF[,1:10], rotationM[,1:10])
fvsm <- cor(rotationF[,2:56], rotationM[,2:56])
diff <- diag(fvsm)


#0.71. Not identical, but pretty similar
#Can I multiply these vectors by the diffrence vector above to compare them in the same supsace? If so, corr goes up but is not perfect. 
cor(diff_mean_F, diff_mean_M)

####Selecting actual individuals. 

################F SELECTION###########
crossF <- filter(femalesPC, pop == "cross")

#Code from Will that is cleaner and faster than mine. 
ShapeFtop <-  crossF$Fly_ID[ order( crossF$diff) ][ 1:50 ]
ShapeFbottom <-crossF$Fly_ID[ order( crossF$diff ) ][(nrow(crossF)-49): nrow(crossF)]  

length(ShapeFtop)
length(ShapeFbottom)

SizeFtop <-  crossF$Fly_ID[ order( crossF$CS) ][ 1:50 ]
SizeFbottom <-crossF$Fly_ID[ order( crossF$CS ) ][(nrow(crossF)-49): nrow(crossF)]  

length(SizeFtop)
length(SizeFbottom)

HighShapeHighSize <- intersect(ShapeFtop, SizeFbottom)
HighShapeLowSize <- intersect(ShapeFtop, SizeFtop)
LowShapeHighSize <- intersect(ShapeFbottom, SizeFbottom)
LowShapeLowSize <- intersect(ShapeFbottom, SizeFtop)

highshape <- ShapeFtop[! ShapeFtop %in% c(HighShapeHighSize,HighShapeLowSize) ]

length(highshape)

lowshape <- ShapeFbottom[!ShapeFbottom %in% c(LowShapeHighSize,LowShapeLowSize)]

highsize <- SizeFbottom[!SizeFbottom %in% c(HighShapeHighSize,LowShapeHighSize) ]

lowsize <- SizeFtop[!SizeFtop %in% c(HighShapeLowSize, LowShapeLowSize)]

crossF$class <- factor(ifelse(crossF$Fly_ID %in% HighShapeLowSize,
                              "HighShapeLowSize", 
                    ifelse(crossF$Fly_ID %in% LowShapeHighSize,
                               "LowShapeHighSize", 
                    ifelse(crossF$Fly_ID %in% HighShapeHighSize,
                               "HighShapeSize", 
                    ifelse(crossF$Fly_ID %in% LowShapeLowSize, 
                                "LowShapeSize",
                    ifelse(crossF$Fly_ID %in% highshape,
                                "HighShape", 
                    ifelse(crossF$Fly_ID %in% lowshape, 
                                "LowShape", 
                    ifelse(crossF$Fly_ID %in% lowsize,
                                "LowSize", 
                    ifelse(crossF$Fly_ID %in% highsize, 
                                "HighSize","NA"
                               )))))))))
levels(crossF$class)

with(crossF, table(class))

ggplot(crossF, aes(CS, diff, col = class)) + 
  geom_point()



crossM <- filter(malesPC, pop == "cross")

#Code from Will that is cleaner and faster than mine. 
ShapeMtop <-  crossM$Fly_ID[ order( crossM$diff) ][ 1:50 ]
ShapeMbottom <-crossM$Fly_ID[ order( crossM$diff ) ][(nrow(crossM)-49): nrow(crossM)]  

length(ShapeMtop)
length(ShapeMbottom)

SizeMtop <-  crossM$Fly_ID[ order( crossM$CS) ][ 1:50 ]
SizeMbottom <-crossM$Fly_ID[ order( crossM$CS ) ][(nrow(crossM)-49): nrow(crossM)]  

length(SizeMtop)
length(SizeMbottom)

#none
HighShapeHighSizeM <- intersect(ShapeMtop, SizeMtop)
HighShapeLowSizeM <- intersect(ShapeMtop, SizeMbottom)
LowShapeHighSizeM <- intersect(ShapeMbottom, SizeMtop)
LowShapeLowSizeM <- intersect(ShapeMbottom, SizeMbottom)


highshapeM <- ShapeMtop[! ShapeMtop %in% c(HighShapeLowSizeM, HighShapeHighSizeM ) ]

length(highshapeM)

lowshapeM <- ShapeMbottom[!ShapeMbottom %in% c(LowShapeHighSizeM, LowShapeLowSizeM )]

highsizeM <- SizeMtop[!SizeMtop %in% c(LowShapeHighSizeM,HighShapeHighSizeM)  ]

lowsizeM <- SizeMbottom[!SizeMbottom %in% c(HighShapeLowSizeM,LowShapeLowSizeM) ]

crossM$class <- factor(ifelse(crossM$Fly_ID %in% HighShapeLowSizeM,
                              "HighShapeLowSize", 
                      ifelse(crossM$Fly_ID %in% LowShapeHighSizeM,
                                     "LowShapeHighSize", 
                       ifelse(crossM$Fly_ID %in%HighShapeHighSizeM,
                                            "HighShapeSize", 
                              ifelse(crossF$Fly_ID %in% LowShapeLowSizeM, 
                                     "LowShapeSize",
                       ifelse(crossM$Fly_ID %in% highshapeM,
                       
                                                   "HighShape", 
                                                   ifelse(crossM$Fly_ID %in% lowshapeM, 
                                                          "LowShape", 
                                                          ifelse(crossM$Fly_ID %in% lowsizeM,
                                                                 "LowSize", 
                                                                 ifelse(crossM$Fly_ID %in% highsizeM, 
                                                                        "HighSize","NA"
                                                                 )))))))))

levels(crossF$class)
                       
with(crossM, table(class))
                       
ggplot(crossM, aes(CS, diff, col = class)) + 
                         geom_point()


#####Now I want to ask about the largest wings in the M data set that look suspect and also about the misclassififed wings by the LDA 

malestoobig <- filter(crossM, CS > 9)[,c("Fly_ID", "CS")]

#So most of them. 
malestoobig$Fly_ID %in% checkFlyM$Fly_ID


write.csv(checkFly, file = "suspectflies_zi192ef81.csv")

highshapeFoutlier <- ShapeFtop[  ShapeFtop %in% checkFly$Fly_ID]
#ok only 1 
length(highshapeFoutlier)
lowshapeFoutlier <- ShapeFbottom[ ShapeFbottom %in% checkFly$Fly_ID]
#6. 
length(lowshapeFoutlier)
lowsizeFoutlier <- SizeFbottom[ SizeFbottom %in% checkFly$Fly_ID]
#3
length(lowsizeFoutlier)
highsizeFoutlier <- SizeFtop[ SizeFtop %in% checkFlyF$Fly_ID]
#12
length(highsizeFoutlier)
highshapeMoutlier <- ShapeMtop[ ShapeMtop %in% checkFlyM$Fly_ID]
#5 
length(highshapeMoutlier)

lowshapeMoutlier <- ShapeMbottom[ShapeMbottom %in% checkFly$Fly_ID]
#6
length(lowshapeMoutlier)
lowsizeMoutlier <- SizeMbottom[ SizeMbottom %in% checkFly$Fly_ID]
#20!!!! 
length(lowsizeMoutlier)
highsizeMoutlier <- SizeMtop[ SizeMtop %in% checkFly$Fly_ID]
#all but 1 
length(highsizeMoutlier)

checkMe <- c(lowsizeMoutlier, lowshapeMoutlier, highshapeMoutlier, highsizeMoutlier, highsizeFoutlier, highsizeFoutlier, lowshapeFoutlier, lowsizeFoutlier )


###########Going to also try this with just landmarks and regressing out size, which is what I wish I had done before. 


#fitting a unique slope for each sex

resid <- lm(as.matrix(wings[,6:101]) ~ CS + CS:sex, 
            dat = wings)$residuals

sizewings <- data.frame(wings, resid)


lmfemale <- filter(sizewings, sex =="F")

HA_meanlm_F <- colMeans(lmfemale[lmfemale$pop == "HA", 105:200])
LA_meanlm_F <- colMeans(lmfemale[lmfemale$pop == "LA", 105:200])
diff_meanlm_F <- HA_meanlm_F - LA_meanlm_F

lmfemale$diff <-  as.matrix(lmfemale[,105:200]) %*% diff_meanlm_F

#I think this looks messier but whatever. maybe more transgressive.
female_plotlm <- ggplot(lmfemale, aes(x = CS, y = diff, col = pop)) +
  geom_point(alpha = 0.6, size = 2 )

ggplot(lmfemale, aes(x = CS, y = diff, col = line)) +
  geom_point(alpha = 0.6, size = 2 )

#Selection 
lmcrossF <- filter(lmfemale, pop == "cross")

#Code from Will that is cleaner and faster than mine. 
lmShapeFtop <-  lmcrossF$Fly_ID[ order( lmcrossF$diff) ][ 1:50 ]
lmShapeFbottom <-lmcrossF$Fly_ID[ order( lmcrossF$diff ) ][(nrow(lmcrossF)-49): nrow(lmcrossF)]  

length(lmShapeFtop)
length(lmShapeFbottom)

lmSizeFbottom <-  lmcrossF$Fly_ID[ order( lmcrossF$CS) ][ 1:50 ]
lmSizeFtop <-lmcrossF$Fly_ID[ order( lmcrossF$CS ) ][(nrow(lmcrossF)-49): nrow(lmcrossF)]  

length(lmSizeFtop)
length(lmSizeFbottom)

lmHighShapeHighSize <- intersect(lmShapeFtop, lmSizeFtop)
lmHighShapeLowSize <- intersect(lmShapeFtop, lmSizeFbottom)
lmLowShapeHighSize <- intersect(lmShapeFbottom, lmSizeFtop)
#No lowshape lowsize
intersect(lmShapeFbottom, lmSizeFbottom)

lmhighshape <- lmShapeFtop[! lmShapeFtop %in% c(lmHighShapeHighSize,lmHighShapeLowSize) ]

length(lmhighshape)

lmlowshape <- lmShapeFbottom[!lmShapeFbottom %in% lmLowShapeHighSize]

lmhighsize <- lmSizeFtop[!lmSizeFtop %in% c(lmHighShapeHighSize,lmLowShapeHighSize) ]

lmlowsize <- lmSizeFbottom[!lmSizeFbottom %in% lmHighShapeLowSize]

lmcrossF$class <- factor(ifelse(lmcrossF$Fly_ID %in% lmHighShapeLowSize,
                              "HighShapeLowSize", 
                              ifelse(lmcrossF$Fly_ID %in% lmLowShapeHighSize,
                                     "LowShapeHighSize", 
                                     ifelse(lmcrossF$Fly_ID %in% lmHighShapeHighSize,
                                            "HighShapeSize", 
                                            ifelse(crossF$Fly_ID %in% lmhighshape,
                                                   "HighShape", 
                                                   ifelse(crossF$Fly_ID %in% lmlowshape, 
                                                          "LowShape", 
                                                          ifelse(crossF$Fly_ID %in% lmlowsize,
                                                                 "LowSize", 
                                                                 ifelse(crossF$Fly_ID %in% lmhighsize, 
                                                                        "HighSize","NA"
                                                                 ))))))))
levels(crossF$class)

with(crossF, table(class))

ggplot(crossF, aes(CS, diff, col = class)) + 
  geom_point()

#going to look at the intersect between the 50 selected individuals and not worry about pooling 

#not identical but pretty similar. 
#88%
length(intersect(ShapeFtop, lmShapeFtop))/50

#94%
length(intersect(ShapeFbottom, lmShapeFbottom))/50

#Should be the same, just to check. Yep. identical. 
length(intersect(SizeFbottom, lmSizeFtop))/50


#######Boys!####################

lmmale <- filter(sizewings, sex =="M")

HA_meanlm_M <- colMeans(lmmale[lmmale$pop == "HA", 105:200])
LA_meanlm_M <- colMeans(lmmale[lmmale$pop == "LA", 105:200])
diff_meanlm_M <- HA_meanlm_M - LA_meanlm_M

lmmale$diff <-  as.matrix(lmmale[,105:200]) %*% diff_meanlm_M

#I think this looks messier but whatever. maybe more transgressive.
male_plotlm <- ggplot(lmmale, aes(x = CS, y = diff, col = pop)) +
  geom_point(alpha = 0.6, size = 2 )

ggplot(lmmale, aes(x = CS, y = diff, col = line)) +
  geom_point(alpha = 0.6, size = 2 )


lmcrossM <- filter(lmmale, pop == "cross")

#Code from Will that is cleaner and faster than mine. 
lmShapeMtop <-  lmcrossM$Fly_ID[ order( lmcrossM$diff) ][ 1:50 ]
lmShapeMbottom <-lmcrossM$Fly_ID[ order( lmcrossM$diff ) ][(nrow(lmcrossM)-49): nrow(lmcrossM)]  

length(lmShapeMtop)
length(lmShapeMbottom)

lmSizeMtop <-  lmcrossM$Fly_ID[ order( lmcrossM$CS) ][ 1:50 ]
lmSizeMbottom <-lmcrossM$Fly_ID[ order( lmcrossM$CS ) ][(nrow(lmcrossM)-49): nrow(lmcrossM)]  

length(lmSizeMtop)
length(lmSizeMbottom)

lmHighShapeHighSizeM <- intersect(lmShapeMtop, lmSizeMtop)
lmHighShapeLowSizeM <- intersect(lmShapeMtop, lmSizeMbottom)
lmLowShapeHighSizeM <- intersect(lmShapeMbottom, lmSizeMtop)
#No lowshape lowsize
intersect(lmShapeFbottom, lmSizeFbottom)

lmhighshapeM <- lmShapeMtop[! lmShapeMtop %in% c(lmHighShapeHighSizeM,lmHighShapeLowSizeM) ]

length(lmhighshapeM)

lmlowshapeM <- lmShapeMbottom[!lmShapeMbottom %in% lmLowShapeHighSizeM]

lmhighsizeM <- lmSizeMtop[!lmSizeMtop %in% c(lmHighShapeHighSizeM,lmLowShapeHighSizeM) ]

lmlowsizeM <- lmSizeMbottom[!lmSizeMbottom %in% lmHighShapeLowSizeM]

lmcrossM$class <- factor(ifelse(lmcrossM$Fly_ID %in% lmHighShapeLowSizeM,
                                "HighShapeLowSize", 
                                ifelse(lmcrossM$Fly_ID %in% lmLowShapeHighSizeM,
                                       "LowShapeHighSize", 
                                       ifelse(lmcrossM$Fly_ID %in% lmHighShapeHighSizeM,
                                              "HighShapeSize", 
                                              ifelse(crossM$Fly_ID %in% lmhighshapeM,
                                                     "HighShape", 
                                                     ifelse(crossM$Fly_ID %in% lmlowshapeM, 
                                                            "LowShape", 
                                                            ifelse(crossM$Fly_ID %in% lmlowsizeM,
                                                                   "LowSize", 
                                                                   ifelse(crossM$Fly_ID %in% lmhighsizeM, 
                                                                          "HighSize","NA"
                                                                   ))))))))
levels(crossM$class)

with(crossM, table(class))

ggplot(crossM, aes(CS, diff, col = class)) + 
  geom_point()

#going to look at the intersect between the 50 selected individuals and not worry about pooling 

#not identical but pretty similar. 
#78%
length(intersect(ShapeMtop, lmShapeMtop))/50

#90%
length(intersect(ShapeMbottom, lmShapeMbottom))/50

#Should be the same, just to check. Yep. identical. 
length(intersect(SizeMbottom, lmSizeMbottom))/50


###comparing seleciton vectors here 
#Not super similar. 
cor(diff_meanlm_F, diff_meanlm_M)


#Going to go ahead with the PCs because 80 - 90% is pretty close and it is what I did with the first population. 
#But for selecting, I need to remove some of the specimens to make things easier becasuse of problems with the samples. 
#Removing all plate 23 and plate 24 because there are 2 plate 23 in the freezer and no plate 24. Plus one was left in the fridge for an unknown period of time. 
#removing plate9f because these are not the right wings. No idea what was imaged but it wasn't this coverslip 

#I am also going to print the top 75 and just go down the list and figure it out because some flies are mis sexed. 

################F SELECTION###########

sel1 <- crossF[ grep("plate23", crossF$Fly_ID, invert = TRUE),]
selF <- sel1[ grep("plate24", sel1$Fly_ID, invert = TRUE),]
#1230 + 45 + 48 = 1323 Looks right 

 
#Code from Will that is cleaner and faster than mine. 
ShapeFtop <-  selF$Fly_ID[ order( selF$diff) ][ 1:50 ]
ShapeFbottom <-selF$Fly_ID[ order( selF$diff ) ][(nrow(selF)-49): nrow(selF)]  

length(ShapeFtop)
length(ShapeFbottom)

SizeFtop <-  selF$Fly_ID[ order( selF$CS) ][ 1:50 ]
SizeFbottom <-selF$Fly_ID[ order( selF$CS ) ][(nrow(selF)-49): nrow(selF)]  

length(SizeFtop)
length(SizeFbottom)

HighShapeHighSize <- intersect(ShapeFtop, SizeFbottom)
HighShapeLowSize <- intersect(ShapeFtop, SizeFtop)
LowShapeHighSize <- intersect(ShapeFbottom, SizeFbottom)
LowShapeLowSize <- intersect(ShapeFbottom, SizeFtop)

highshape <- ShapeFtop[! ShapeFtop %in% c(HighShapeHighSize,HighShapeLowSize) ]

length(highshape)

lowshape <- ShapeFbottom[!ShapeFbottom %in% c(LowShapeHighSize,LowShapeLowSize)]

highsize <- SizeFbottom[!SizeFbottom %in% c(HighShapeHighSize,LowShapeHighSize) ]

lowsize <- SizeFtop[!SizeFtop %in% c(HighShapeLowSize, LowShapeLowSize)]

selF$class <- factor(ifelse(selF$Fly_ID %in% HighShapeLowSize,
                              "HighShapeLowSize", 
                              ifelse(selF$Fly_ID %in% LowShapeHighSize,
                                     "LowShapeHighSize", 
                                     ifelse(selF$Fly_ID %in% HighShapeHighSize,
                                            "HighShapeSize", 
                                            ifelse(selF$Fly_ID %in% LowShapeLowSize, 
                                                   "LowShapeSize",
                                                   ifelse(selF$Fly_ID %in% highshape,
                                                          "HighShape", 
                                                          ifelse(selF$Fly_ID %in% lowshape, 
                                                                 "LowShape", 
                                                                 ifelse(selF$Fly_ID %in% lowsize,
                                                                        "LowSize", 
                                                                        ifelse(selF$Fly_ID %in% highsize, 
                                                                               "HighSize","NA"
                                                                        )))))))))
levels(selF$class)

with(selF, table(class))

ggplot(selF, aes(CS, diff, col = class)) + 
  geom_point()



#########Males################

sel1 <- crossM[ grep("plate23", crossM$Fly_ID, invert = T),]
sel2 <- sel1[ grep("plate24", sel1$Fly_ID, invert = T),]
selM <- sel2[ grep("plate9f", sel2$Fly_ID, invert = T),]
#1213 + 45 + 46 +11 = 1315 Looks right 


 
ShapeMtop <-  selM$Fly_ID[ order( selM$diff) ][ 1:50 ]
ShapeMbottom <-selM$Fly_ID[ order( selM$diff ) ][(nrow(selM)-49): nrow(selM)]  

length(ShapeMtop)
length(ShapeMbottom)

SizeMbottom <-  selM$Fly_ID[ order( selM$CS) ][ 1:50 ]
crap <-  selM$Fly_ID[ order( selM$CS) ]
SizeMtop <-selM$Fly_ID[ order( selM$CS ) ][(nrow(selM)-49): nrow(selM)]  

length(SizeMtop)
length(SizeMbottom)

#none
HighShapeHighSizeM <- intersect(ShapeMtop, SizeMtop)
HighShapeLowSizeM <- intersect(ShapeMtop, SizeMbottom)
LowShapeHighSizeM <- intersect(ShapeMbottom, SizeMtop)
LowShapeLowSizeM <- intersect(ShapeMbottom, SizeMbottom)


highshapeM <- ShapeMtop[! ShapeMtop %in% c(HighShapeLowSizeM, HighShapeHighSizeM ) ]

length(highshapeM)

lowshapeM <- ShapeMbottom[!ShapeMbottom %in% c(LowShapeHighSizeM, LowShapeLowSizeM )]

highsizeM <- SizeMtop[!SizeMtop %in% c(LowShapeHighSizeM,HighShapeHighSizeM)  ]

lowsizeM <- SizeMbottom[!SizeMbottom %in% c(HighShapeLowSizeM,LowShapeLowSizeM) ]

selM$class <- factor(ifelse(selM$Fly_ID %in% HighShapeLowSizeM,
                              "HighShapeLowSize", 
                              ifelse(selM$Fly_ID %in% LowShapeHighSizeM,
                                     "LowShapeHighSize", 
                                     ifelse(selM$Fly_ID %in%HighShapeHighSizeM,
                                            "HighShapeSize", 
                                            ifelse(selM$Fly_ID %in% LowShapeLowSizeM, 
                                                   "LowShapeSize",
                                                   ifelse(selM$Fly_ID %in% highshapeM,
                                                          
                                                          "HighShape", 
                                                          ifelse(selM$Fly_ID %in% lowshapeM, 
                                                                 "LowShape", 
                                                                 ifelse(selM$Fly_ID %in% lowsizeM,
                                                                        "LowSize", 
                                                                        ifelse(selM$Fly_ID %in% highsizeM, 
                                                                               "HighSize","NA"
                                                                        )))))))))

levels(selM$class)

with(selM, table(class))

ggplot(selM, aes(CS, diff, col = class)) + 
  geom_point()

#Now I need to rank the big and small lists. 


pullM <- (selM %>%
            filter(class != "NA") %>%
            dplyr::select(Fly_ID, CS, diff, class)
          )

write.csv(pullM, "../Outputs/zi192ef81_males.csv")


pullF <- (selF %>%
            filter(class != "NA") %>%
            dplyr::select(Fly_ID, CS, diff, class)
)

write.csv(pullF, "../Outputs/zi192ef81_females.csv")




#######plotting for paper#######


m.sel.plot <- ggplot(crossM, aes(x = CS, y = diff, col = class)) +
  geom_point(alpha = 0.5, size = 0.75) + 
  ylab("Shape Score") + 
  xlab("Centroid Size") + 
  scale_colour_manual(values=c('red', 'red', 'red', 'red', 'red', 'red', 'darkgrey')) +
  geom_vline(aes(xintercept = max(crossM[crossM$Fly_ID %in% SizeMbottom,]$CS))) + 
  geom_vline(aes(xintercept = min(crossM[crossM$Fly_ID %in% SizeMtop,]$CS))) +
  geom_hline(aes(yintercept = min(crossM[crossM$Fly_ID %in% ShapeMbottom,]$diff))) +
  geom_hline(aes(yintercept = max(crossM[crossM$Fly_ID %in% ShapeMtop,]$diff))) +
  theme(panel.background = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        legend.key = element_blank(),
        text = element_text(size=15),
        axis.text.x= element_text(size=12), 
        axis.text.y= element_text(size=12))


f.sel.plot <-ggplot(crossF, aes(x = CS, y = diff, col = class)) +
  geom_point(alpha = 0.5, size = 0.75) + 
  ylab("Shape Score") + 
  xlab("Centroid Size") + 
  scale_colour_manual(values=c('red', 'red', 'red', 'red', 'red', 
                               'red', 'red','darkgrey')) +
  geom_vline(aes(xintercept = min(crossF[crossF$Fly_ID %in% SizeFbottom,]$CS))) + 
  geom_vline(aes(xintercept = max(crossF[crossF$Fly_ID %in% SizeFtop,]$CS))) +
  geom_hline(aes(yintercept = min(crossF[crossF$Fly_ID %in% ShapeFbottom,]$diff))) +
  geom_hline(aes(yintercept = max(crossF[crossF$Fly_ID %in% ShapeFtop,]$diff))) +
  theme(panel.background = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        legend.key = element_blank(),
        text = element_text(size=15),
        axis.text.x= element_text(size=12), 
        axis.text.y= element_text(size=12))



png("../Figures/zi192ef81_bothSex_selection.png", 
    width =1700, height = 750, 
    units = "px", res = 300)

plot_grid(f.sel.plot, m.sel.plot, labels = c("Female", "Male"))

dev.off()
