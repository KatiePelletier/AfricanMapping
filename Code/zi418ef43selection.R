#This is the selection for zi418xef43 wings for sequencing, following the same protocol as zi192xef96 cross submitted summer 2020

library(tidyverse)

library(MASS)

######### Loading functions #########
source('~/Dropbox/DworkinLabSharedMaterial/scripts/WRP_FUNCTIONS.R', chdir = TRUE)

source('~/Dropbox/DworkinLabSharedMaterial/scripts/WINGPLOTSOURCE.R', chdir = TRUE)

####################################

# #read in dat
# wings_raw <- read_tsv("zi418ef43_Output_Landmarks+Veins+Outline.dat")
# 
# #need to seperate out PL and cross wings
# crosswings <- filter(wings_raw, grepl("plate", wings_raw$CPFile))
# PLwings <- filter(wings_raw, !grepl("plate", wings_raw$CPFile))
# #602+2415 = 3017
# 
# Crap_Col <- c("File", "O1x", "O1y", "O2x", "O2y", "Date", "Time", "Tags", "Sex", "Perp")
# More_crap <- c("Wings", "KP", "AfPop", "PL", "folder", "folder2")
# 
# 
# 
# PL_fixed <- (PLwings
#                 %>% separate("CPFile", into= c("perp", "line", "slide", "Real_Sex", "Fly"), sep = "_")
#                 %>% dplyr::select(-Crap_Col)
#                 %>% unite(Fly_ID, slide, Fly)
#                )
# 
# PL_fixed$pop <- ifelse(grepl("ZI", PL_fixed$line), "LA", "HA")
# levels(as.factor(PL_fixed$pop))
# 
# cross_fixed <- (crosswings
#                 %>%separate("CPFile", into= c("perp", "line", "colection", "slide", "Real_Sex", "Fly"), sep = "_")
#                 %>% dplyr::select(-Crap_Col)
#                 %>% dplyr::select(-colection)
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
# save(wings, file = "zi418ef43clean.rda")

load("../Data/zi418ef43clean.rda")

wings$logCS <- log(wings$CS)

#Looks like wings. Nothing weird. 
WingBlur(wings[sample(nrow(wings), 50) ,6:101])
#WingBlur(wings[,6:101])

#calculating the PCs

PCs <- prcomp(wings[, c(6:101)])
wings_combined <- data.frame(wings, PCs$x[,1:57]) 

#Went back and checked. They are just big
ggplot(wings_combined, aes(x = logCS, y = PC1, col = pop)) + 
  geom_point()

#taking a look at size alone 
cross <- filter(wings_combined, line == "Zi418Ef43")

#females are maybe a little bimodal. I think what we saw earlier (Archana data) was the plotting because ggplot his stacks by defalt
#Will ask about the really small females below. 
ggplot(cross, aes( x = CS, col = sex)) + 
  geom_freqpoly(binwidth = 0.05)

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
cor(Fdiff, Mdiff) #0.93

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
#almost same answer. Of course. 
cor(lmFdiff, lmMdiff) #0.92

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
# #this seems...Still pretty low. 90%
 #Should double check these wings when I make the pools.
 (1115+1048)/nrow(test)
 
 
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

#I want to correct for the LDA predicted sex in the 
wings$sex2 <- NA


for (i in 1:nrow(wings)) {
  ID <- droplevels(wings$Fly_ID[i])
  if (wings$Fly_ID[i] %in% checkFly$Fly_ID) {
    wings$sex2[i] <- checkFly[ID == checkFly$Fly_ID,]$pred
  } else {
    wings$sex2[i] <- wings$sex[i]
  }
}
 
#this codes as 2 and 1 rather than M and F. but F == 1
wings$sex2

hist(wings$sex2)

#some are probably still wrong (although more males just may have been removed because of bad dissections) Better than before. 
wings$sex.corr <- as.factor(ifelse(wings$sex2 == 1, "F", "M"))
plot(wings$sex.corr)


#for the previous popualtion, we chose to look at the males and females seperatley. So I will do that again. 

females <- filter(wings, sex.corr == "F")
#Using the landmarks plus logCS
PCAf <- prcomp(females[,c(6:101, 104)])

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

#plate9c_08 spline is messed up. Needs to be fixed. 
#check plate9c_08 as well. 
#Other outliers have been checked at least twice.
WingPlot(as.matrix(filter(wings, Fly_ID == "plate9c_09"))[,6:101])


#for now I am going to remove this because the PC is in use.
dim(wings)
wings <- wings[wings$Fly_ID != "plate9c_09",]
dim(wings)

#not going to use PC1 because it captures the allometric componet. (might try to do this again below with the landmarks and modeling out size to see what happens) 


HA_mean_F <- colMeans(femalesPC[femalesPC$pop == "HA", 108:162])
LA_mean_F <- colMeans(femalesPC[femalesPC$pop == "LA", 108:162])
diff_mean_F <- HA_mean_F - LA_mean_F

femalesPC$diff <-  as.matrix(femalesPC[,108:162]) %*% diff_mean_F

female_plot <- ggplot(femalesPC, aes(x = CS, y = diff, col = pop)) +
  geom_point(alpha = 0.6, size = 2 )

ggplot(femalesPC, aes(x = CS, y = diff, col = line)) +
  geom_point(alpha = 0.6, size = 2 ) +
  geom_text(aes(label=as.character(Fly_ID),hjust=0,vjust=0))



#########Males#################

males <- filter(wings, sex.corr == "M")
#Using the landmarks plus logCS
PCAm <- prcomp(males[,c(6:101, 104)])

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


#Outliers look fine here and as splines. 
WingPlot(as.matrix(filter(wings, Fly_ID == "plate24h_03"))[,6:101])

#not going to use PC1 because it captures the allometric componet. (might try to do this again below with the landmarks and modeling out size to see what happens) 
HA_mean_M <- colMeans(malesPC[malesPC$pop == "HA", 108:162])
LA_mean_M <- colMeans(malesPC[malesPC$pop == "LA", 108:162])
diff_mean_M <- HA_mean_M - LA_mean_M

malesPC$diff <-  as.matrix(malesPC[,108:162]) %*% diff_mean_M

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


#This one is low. but I think they are in oppisite directions. 
#Can I multiply these vectors by the diffrence vector above to compare them in the same supsace? If so, corr goes up but is not perfect. 
cor(diff_mean_F, diff_mean_M)

####Selecting actual individuals. 

# ################F SELECTION###########
crossF <- filter(femalesPC, pop == "cross")

#Code from Will that is cleaner and faster than mine.
ShapeFtop <-  crossF$Fly_ID[ order( crossF$diff) ][ 1:50 ]
ShapeFbottom <-crossF$Fly_ID[ order( crossF$diff ) ][(nrow(crossF)-49): nrow(crossF)]

#need 11 more because they were males. 
ShapeFbottom2 <-crossF$Fly_ID[ order( crossF$diff ) ][(nrow(crossF)-60): nrow(crossF)]

newShapeFBottom <- ShapeFbottom2[!(ShapeFbottom2 %in% ShapeFbottom)] 
length(newShapeBottom)

length(ShapeFtop)
length(ShapeFbottom)

#and Again. How fun. 
ShapeFbottom3 <-crossF$Fly_ID[ order( crossF$diff ) ][(nrow(crossF)-65): nrow(crossF)]

fuckthis <- ShapeFbottom3[!(ShapeFbottom3 %in% ShapeFbottom2)] 
length(fuckthis)
fuckthis

length(ShapeFtop)
length(ShapeFbottom)

#actually the small ones.
SizeFtop <-  crossF$Fly_ID[ order( crossF$CS) ][ 1:50 ]
#actually the big ones.
SizeFbottom <-crossF$Fly_ID[ order( crossF$CS ) ][(nrow(crossF)-49): nrow(crossF)]

#need 38 more because they were males. 
SizeFbottom2 <-crossF$Fly_ID[ order( crossF$CS ) ][(nrow(crossF)-91): nrow(crossF)]

newSizeFBottom <- SizeFbottom2[!(SizeFbottom2 %in% SizeFbottom)] 
length(newSizeBottom)

length(SizeFtop)
length(SizeFbottom)

HighShapeHighSize <- intersect(ShapeFtop, SizeFbottom)
HighShapeLowSize <- intersect(ShapeFtop, SizeFtop)
LowShapeHighSize <- intersect(ShapeFbottom, SizeFbottom)
LowShapeLowSize <- intersect(ShapeFbottom, SizeFtop)

#I don't want these to intersect with things I already pulled. 
#3 of these wont work. 
#added 3 and 1 of those till wont work.  
# exclude these ones from new size OH WELL "plate26c_06" "plate1a_08"  "plate26a_08" "plate20c_04"
newHighShapeHighSize <- intersect(ShapeFtop, newSizeFBottom)
#newHighShapeLowSize <- intersect(ShapeFtop, SizeFtop)
newLowShapeHighSize <- intersect(newShapeFBottom, newSizeFBottom)
newLowShapeLowSize <- intersect(newShapeFBottom, SizeFtop)

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

pullF <-(crossF %>%
           filter(class != "NA") %>%
           dplyr::select(Fly_ID, diff, CS, class))


flyF <- c(as.character(newShapeFBottom), as.character(newSizeFBottom))
classF <- c(rep("LowShape", length(newShapeFBottom)), rep("LowSize", length(newSizeFBottom)) )

pullsFnew <- cbind(flyF, classF)
write.csv(pullsFnew, "zi418ef43_extras_pullsF.csv")

#Writing it out
write.csv(pullF, "zi418ef43_withSexCorrection_pullsF.csv")

######################################################

crossM <- filter(malesPC, pop == "cross")

#Code from Will that is cleaner and faster than mine.
ShapeMtop <-  crossM$Fly_ID[ order( crossM$diff) ][ 1:50 ]
ShapeMbottom <-crossM$Fly_ID[ order( crossM$diff ) ][(nrow(crossM)-49): nrow(crossM)]

MshapeMax <- max(crossM$diff[ 1:50 ])
MshapeMin <- min(crossM$diff[nrow(crossM)-49: nrow(crossM)])

length(ShapeMtop)
length(ShapeMbottom)


ShapeMtop2 <-  crossM$Fly_ID[ order( crossM$diff) ][ 1:53 ]
ShapeMbottom2 <-crossM$Fly_ID[ order( crossM$diff ) ][(nrow(crossM)-50): nrow(crossM)]

newShapeMbottom <- ShapeMbottom2[!(ShapeMbottom2 %in% ShapeMbottom)] 
length(newShapeMbottom)

newShapeMtop <- ShapeMtop2[!(ShapeMtop2 %in% ShapeMtop)] 
length(newShapeMtop)

SizeMtop <-  crossM$Fly_ID[ order( crossM$CS) ][ 1:50 ]
SizeMbottom <-crossM$Fly_ID[ order( crossM$CS ) ][(nrow(crossM)-49): nrow(crossM)]

MsizeMax <- max(crossM$CS[ 1:50 ])
MsizeMin <- min(crossM$CS[nrow(crossM)-49: nrow(crossM)])

length(SizeMtop)
length(SizeMbottom)

SizeMtop2 <-  crossM$Fly_ID[ order( crossM$CS) ][ 1:55 ]
SizeMbottom2 <-crossM$Fly_ID[ order( crossM$CS ) ][(nrow(crossM)-66): nrow(crossM)]

#and three more for clean up 
SizeMbottom3 <-crossM$Fly_ID[ order( crossM$CS ) ][(nrow(crossM)-70): nrow(crossM)]
crap <- SizeMbottom3[!(SizeMbottom3 %in% SizeMbottom2)] 
length(crap)
crap


newSizeMbottom <- SizeMbottom2[!(SizeMbottom2 %in% SizeMbottom)] 
length(newSizeMbottom)

#double checking these are big
#These look right
crossM[crossM$Fly_ID %in% newSizeMbottom,]$CS

newSizeMtop <- SizeMtop2[!(SizeMtop2 %in% SizeMtop)] 
length(newSizeMtop)
crossM[crossM$Fly_ID %in% newSizeMtop,]$CS

#are these in other groups already?
intersect(newShapeMtop, SizeMtop) #0
intersect(ShapeMtop, newSizeMtop) #0 
intersect(newShapeMtop, SizeMbottom) #0
intersect(ShapeMtop, newSizeMbottom) #1 "plate13f_1"
intersect(newShapeMbottom, SizeMtop) #0
intersect(ShapeMbottom, newSizeMtop) #1 "plate27d_01"
intersect(newShapeMbottom, SizeMbottom)#0
intersect(ShapeMbottom, newSizeMbottom) #0

intersect(newShapeMtop, newSizeMtop) #0
intersect(newShapeMtop, newSizeMbottom) #0
intersect(newShapeMbottom, newSizeMtop) #0
intersect(newShapeMbottom, newSizeMbottom) #0





HighShapeHighSizeM <- intersect(newShapeMtop, SizeMtop)
HighShapeHighSizeM <- intersect(newShapeMtop, SizeMtop)
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
                              ifelse(crossM$Fly_ID %in% LowShapeLowSizeM,
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

levels(crossM$class)

with(crossM, table(class))

ggplot(crossM, aes(CS, diff, col = class)) +
                         geom_point()


pullM <-(crossM %>%
           filter(class != "NA") %>%
           dplyr::select(Fly_ID, diff, CS, class))

flyM <- c(as.character(newShapeMtop), as.character(newShapeMbottom), as.character(newSizeMtop), as.character(newSizeMbottom) )
classM <- c(rep("HighShape", length(newShapeMtop)), rep("LowShape", length(newShapeMbottom)), rep("HighSize", length(newSizeMtop)), rep("LowSize", length(newSizeMbottom)) )

pullsMnew <- cbind(flyM, classM)
write.csv(pullsMnew, "zi418ef43_extras_pullsM.csv")

#Writing it out
write.csv(pullM, "zi418ef43_withSexCorrection_pullsM.csv")

#I decided to use the predicted sex because the LDA does a pretty good job and it means I can make fewer corrections later.

# ################F SELECTION###########
# crossF <- filter(femalesPC, pop == "cross")
# 
# #Code from Will that is cleaner and faster than mine. 
# ShapeFtop <-  crossF$Fly_ID[ order( crossF$diff) ][ 1:50 ]
# ShapeFbottom <-crossF$Fly_ID[ order( crossF$diff ) ][(nrow(crossF)-49): nrow(crossF)]  
# 
# length(ShapeFtop)
# length(ShapeFbottom)
# 
# #actually the small ones. 
# SizeFtop <-  crossF$Fly_ID[ order( crossF$CS) ][ 1:50 ]
# #actually the big ones. 
# SizeFbottom <-crossF$Fly_ID[ order( crossF$CS ) ][(nrow(crossF)-49): nrow(crossF)]  
# 
# length(SizeFtop)
# length(SizeFbottom)
# 
# HighShapeHighSize <- intersect(ShapeFtop, SizeFbottom)
# HighShapeLowSize <- intersect(ShapeFtop, SizeFtop)
# LowShapeHighSize <- intersect(ShapeFbottom, SizeFbottom)
# LowShapeLowSize <- intersect(ShapeFbottom, SizeFtop)
# 
# highshape <- ShapeFtop[! ShapeFtop %in% c(HighShapeHighSize,HighShapeLowSize) ]
# 
# length(highshape)
# 
# lowshape <- ShapeFbottom[!ShapeFbottom %in% c(LowShapeHighSize,LowShapeLowSize)]
# 
# highsize <- SizeFbottom[!SizeFbottom %in% c(HighShapeHighSize,LowShapeHighSize) ]
# 
# lowsize <- SizeFtop[!SizeFtop %in% c(HighShapeLowSize, LowShapeLowSize)]
# 
# crossF$class <- factor(ifelse(crossF$Fly_ID %in% HighShapeLowSize,
#                               "HighShapeLowSize", 
#                               ifelse(crossF$Fly_ID %in% LowShapeHighSize,
#                                      "LowShapeHighSize", 
#                                      ifelse(crossF$Fly_ID %in% HighShapeHighSize,
#                                             "HighShapeSize", 
#                                             ifelse(crossF$Fly_ID %in% LowShapeLowSize, 
#                                                    "LowShapeSize",
#                                                    ifelse(crossF$Fly_ID %in% highshape,
#                                                           "HighShape", 
#                                                           ifelse(crossF$Fly_ID %in% lowshape, 
#                                                                  "LowShape", 
#                                                                  ifelse(crossF$Fly_ID %in% lowsize,
#                                                                         "LowSize", 
#                                                                         ifelse(crossF$Fly_ID %in% highsize, 
#                                                                                "HighSize","NA"
#                                                                         )))))))))
# levels(crossF$class)
# 
# with(crossF, table(class))
# 
# ggplot(crossF, aes(CS, diff, col = class)) + 
#   geom_point()
# 
# ######################################################
# 
# crossM <- filter(malesPC, pop == "cross")
# 
# #Code from Will that is cleaner and faster than mine. 
# ShapeMtop <-  crossM$Fly_ID[ order( crossM$diff) ][ 1:50 ]
# ShapeMbottom <-crossM$Fly_ID[ order( crossM$diff ) ][(nrow(crossM)-49): nrow(crossM)]  
# 
# length(ShapeMtop)
# length(ShapeMbottom)
# 
# SizeMtop <-  crossM$Fly_ID[ order( crossM$CS) ][ 1:50 ]
# SizeMbottom <-crossM$Fly_ID[ order( crossM$CS ) ][(nrow(crossM)-49): nrow(crossM)]  
# 
# length(SizeMtop)
# length(SizeMbottom)
# 
# #none
# HighShapeHighSizeM <- intersect(ShapeMtop, SizeMtop)
# HighShapeLowSizeM <- intersect(ShapeMtop, SizeMbottom)
# LowShapeHighSizeM <- intersect(ShapeMbottom, SizeMtop)
# LowShapeLowSizeM <- intersect(ShapeMbottom, SizeMbottom)
# 
# 
# highshapeM <- ShapeMtop[! ShapeMtop %in% c(HighShapeLowSizeM, HighShapeHighSizeM ) ]
# 
# length(highshapeM)
# 
# lowshapeM <- ShapeMbottom[!ShapeMbottom %in% c(LowShapeHighSizeM, LowShapeLowSizeM )]
# 
# highsizeM <- SizeMtop[!SizeMtop %in% c(LowShapeHighSizeM,HighShapeHighSizeM)  ]
# 
# lowsizeM <- SizeMbottom[!SizeMbottom %in% c(HighShapeLowSizeM,LowShapeLowSizeM) ]
# 
# crossM$class <- factor(ifelse(crossM$Fly_ID %in% HighShapeLowSizeM,
#                               "HighShapeLowSize", 
#                               ifelse(crossM$Fly_ID %in% LowShapeHighSizeM,
#                                      "LowShapeHighSize", 
#                                      ifelse(crossM$Fly_ID %in%HighShapeHighSizeM,
#                                             "HighShapeSize", 
#                                             ifelse(crossM$Fly_ID %in% LowShapeLowSizeM, 
#                                                    "LowShapeSize",
#                                                    ifelse(crossM$Fly_ID %in% highshapeM,
#                                                           
#                                                           "HighShape", 
#                                                           ifelse(crossM$Fly_ID %in% lowshapeM, 
#                                                                  "LowShape", 
#                                                                  ifelse(crossM$Fly_ID %in% lowsizeM,
#                                                                         "LowSize", 
#                                                                         ifelse(crossM$Fly_ID %in% highsizeM, 
#                                                                                "HighSize","NA"
#                                                                         )))))))))
# 
# levels(crossM$class)
# 
# with(crossM, table(class))
# 
# ggplot(crossM, aes(CS, diff, col = class)) + 
#   geom_point()


#####Now I want to ask about the largest wings in the M data set that look suspect and also about the misclassififed wings by the LDA 

malestoobig <- filter(crossM, CS > 9.5)[,c("Fly_ID", "CS")]

#13/23
#At least the largest ones are sexed incorretly 
nrow(malestoobig)
sum(malestoobig$Fly_ID %in% checkFly$Fly_ID)


write.csv(checkFly, file = "suspectflies_zi418ef43.csv")

highshapeFoutlier <- ShapeFtop[  ShapeFtop %in% checkFly$Fly_ID]

length(highshapeFoutlier)
lowshapeFoutlier <- ShapeFbottom[ ShapeFbottom %in% checkFly$Fly_ID]
#7. 
length(lowshapeFoutlier)
lowsizeFoutlier <- SizeFbottom[ SizeFbottom %in% checkFly$Fly_ID]
#3
length(lowsizeFoutlier)
highsizeFoutlier <- SizeFtop[ SizeFtop %in% checkFly$Fly_ID]
#21
length(highsizeFoutlier)
highshapeMoutlier <- ShapeMtop[ ShapeMtop %in% checkFly$Fly_ID]
#23
length(highshapeMoutlier)

lowshapeMoutlier <- ShapeMbottom[ShapeMbottom %in% checkFly$Fly_ID]
#3
length(lowshapeMoutlier)
lowsizeMoutlier <- SizeMbottom[ SizeMbottom %in% checkFly$Fly_ID]
#26
length(lowsizeMoutlier)
highsizeMoutlier <- SizeMtop[ SizeMtop %in% checkFly$Fly_ID]
#3
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

#92%
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
#74%
length(intersect(ShapeMtop, lmShapeMtop))/50

#80%
length(intersect(ShapeMbottom, lmShapeMbottom))/50

#Should be the same, just to check. Yep. identical. 
length(intersect(SizeMbottom, lmSizeMbottom))/50


###comparing seleciton vectors here 
#Not super similar. 
cor(diff_meanlm_F, diff_meanlm_M)

###############I also want to ask if I would get the same answers if I use size data from AT thesis (not the CS size) 

#this is imperfect data becasue there are still duplicate values in here. 
ATsize <- read_csv("sizedat_forKP_zi418ed43.csv")

hist(ATsize$Area)

#I need to drop out the FlyIDs for unsplinable wings 

size.dat <- ATsize[ATsize$FlyID %in% cross$Fly_ID,]

#these are not the same at all. About 500 wings are missing??????
#Her data set has about 100 fewer wings than mine to start. She must have dropped out some wings that I splined. 
nrow(size.dat)
nrow(cross)
nrow(ATsize)

crap <- cross[ATsize$FlyID %in% ATsize$FlyID,]
nrow(crap)
crap2 <- cross[ATsize$FlyID %in% cross$Fly_ID,]

plot(crap$Area, crap2$CS)
cor(crap$Area, crap2$CS)

size.f <- filter(size.dat, sex == "f")
size.f.top <-  size.f$FlyID[order(size.f$Area) ][ 1:50 ]
size.f.bottom <-size.f$FlyID[ order( size.f$Area ) ][(nrow(size.f)-49): nrow(size.f)]  





#not awesome. 
length(intersect(size.f.top, SizeFtop))
length(intersect(size.f.bottom, SizeFbottom))


size.m <- filter(size.dat, sex == "m")

size.m.top <-  size.m$FlyID[order(size.m$Area) ][ 1:50 ]
size.m.bottom <-size.m$FlyID[ order( size.m$Area ) ][(nrow(size.m)-49): nrow(size.m)]  

#Also not awesome 
length(intersect(size.m.top, SizeMtop))
length(intersect(size.m.bottom, SizeMbottom))


#Going to just take the biggest and smallest from the AT data
AT_f <- filter(size.dat, sex == "f")
AT_m <- filter(size.dat, sex == "m")

ATbigF <-  AT_f$FlyID[ order(AT_f$Area) ][ 1:50 ]
ATsmallF <-AT_f$FlyID[ order(AT_f$Area) ][(nrow(AT_f)-49): nrow(AT_f)] 

#how many of these are not splined?



# Keeping this code from the last selection to edit what I did for this cross but I want to make other edits to this data (obviously)
# 
# #I am also going to print the top 75 and just go down the list and figure it out because some flies are mis sexed. 
# 
# ################F SELECTION###########
# 
# sel1 <- crossF[ grep("plate23", crossF$Fly_ID, invert = TRUE),]
# selF <- sel1[ grep("plate24", sel1$Fly_ID, invert = TRUE),]
# #1230 + 45 + 48 = 1323 Looks right 
# 
#  
# #Code from Will that is cleaner and faster than mine. 
# ShapeFtop <-  selF$Fly_ID[ order( selF$diff) ][ 1:75 ]
# ShapeFbottom <-selF$Fly_ID[ order( selF$diff ) ][(nrow(selF)-74): nrow(selF)]  
# 
# length(ShapeFtop)
# length(ShapeFbottom)
# 
# SizeFtop <-  selF$Fly_ID[ order( selF$CS) ][ 1:75 ]
# SizeFbottom <-selF$Fly_ID[ order( selF$CS ) ][(nrow(selF)-74): nrow(selF)]  
# 
# length(SizeFtop)
# length(SizeFbottom)
# 
# HighShapeHighSize <- intersect(ShapeFtop, SizeFbottom)
# HighShapeLowSize <- intersect(ShapeFtop, SizeFtop)
# LowShapeHighSize <- intersect(ShapeFbottom, SizeFbottom)
# LowShapeLowSize <- intersect(ShapeFbottom, SizeFtop)
# 
# highshape <- ShapeFtop[! ShapeFtop %in% c(HighShapeHighSize,HighShapeLowSize) ]
# 
# length(highshape)
# 
# lowshape <- ShapeFbottom[!ShapeFbottom %in% c(LowShapeHighSize,LowShapeLowSize)]
# 
# highsize <- SizeFbottom[!SizeFbottom %in% c(HighShapeHighSize,LowShapeHighSize) ]
# 
# lowsize <- SizeFtop[!SizeFtop %in% c(HighShapeLowSize, LowShapeLowSize)]
# 
# selF$class <- factor(ifelse(selF$Fly_ID %in% HighShapeLowSize,
#                               "HighShapeLowSize", 
#                               ifelse(selF$Fly_ID %in% LowShapeHighSize,
#                                      "LowShapeHighSize", 
#                                      ifelse(selF$Fly_ID %in% HighShapeHighSize,
#                                             "HighShapeSize", 
#                                             ifelse(selF$Fly_ID %in% LowShapeLowSize, 
#                                                    "LowShapeSize",
#                                                    ifelse(selF$Fly_ID %in% highshape,
#                                                           "HighShape", 
#                                                           ifelse(selF$Fly_ID %in% lowshape, 
#                                                                  "LowShape", 
#                                                                  ifelse(selF$Fly_ID %in% lowsize,
#                                                                         "LowSize", 
#                                                                         ifelse(selF$Fly_ID %in% highsize, 
#                                                                                "HighSize","NA"
#                                                                         )))))))))
# levels(selF$class)
# 
# with(selF, table(class))
# 
# ggplot(selF, aes(CS, diff, col = class)) + 
#   geom_point()
# 
# 
# 
# #########Males################
# 
# sel1 <- crossM[ grep("plate23", crossM$Fly_ID, invert = T),]
# sel2 <- sel1[ grep("plate24", sel1$Fly_ID, invert = T),]
# selM <- sel2[ grep("plate9f", sel2$Fly_ID, invert = T),]
# #1213 + 45 + 46 +11 = 1315 Looks right 
# 
# 
#  
# ShapeMtop <-  selM$Fly_ID[ order( selM$diff) ][ 1:75 ]
# ShapeMbottom <-selM$Fly_ID[ order( selM$diff ) ][(nrow(selM)-74): nrow(selM)]  
# 
# length(ShapeMtop)
# length(ShapeMbottom)
# 
# SizeMbottom <-  selM$Fly_ID[ order( selM$CS) ][ 1:75 ]
# crap <-  selM$Fly_ID[ order( selM$CS) ]
# SizeMtop <-selM$Fly_ID[ order( selM$CS ) ][(nrow(selM)-74): nrow(selM)]  
# 
# length(SizeMtop)
# length(SizeMbottom)
# 
# #none
# HighShapeHighSizeM <- intersect(ShapeMtop, SizeMtop)
# HighShapeLowSizeM <- intersect(ShapeMtop, SizeMbottom)
# LowShapeHighSizeM <- intersect(ShapeMbottom, SizeMtop)
# LowShapeLowSizeM <- intersect(ShapeMbottom, SizeMbottom)
# 
# 
# highshapeM <- ShapeMtop[! ShapeMtop %in% c(HighShapeLowSizeM, HighShapeHighSizeM ) ]
# 
# length(highshapeM)
# 
# lowshapeM <- ShapeMbottom[!ShapeMbottom %in% c(LowShapeHighSizeM, LowShapeLowSizeM )]
# 
# highsizeM <- SizeMtop[!SizeMtop %in% c(LowShapeHighSizeM,HighShapeHighSizeM)  ]
# 
# lowsizeM <- SizeMbottom[!SizeMbottom %in% c(HighShapeLowSizeM,LowShapeLowSizeM) ]
# 
# selM$class <- factor(ifelse(selM$Fly_ID %in% HighShapeLowSizeM,
#                               "HighShapeLowSize", 
#                               ifelse(selM$Fly_ID %in% LowShapeHighSizeM,
#                                      "LowShapeHighSize", 
#                                      ifelse(selM$Fly_ID %in%HighShapeHighSizeM,
#                                             "HighShapeSize", 
#                                             ifelse(selM$Fly_ID %in% LowShapeLowSizeM, 
#                                                    "LowShapeSize",
#                                                    ifelse(selM$Fly_ID %in% highshapeM,
#                                                           
#                                                           "HighShape", 
#                                                           ifelse(selM$Fly_ID %in% lowshapeM, 
#                                                                  "LowShape", 
#                                                                  ifelse(selM$Fly_ID %in% lowsizeM,
#                                                                         "LowSize", 
#                                                                         ifelse(selM$Fly_ID %in% highsizeM, 
#                                                                                "HighSize","NA"
#                                                                         )))))))))
# 
# levels(selM$class)
# 
# with(selM, table(class))
# 
# ggplot(selM, aes(CS, diff, col = class)) + 
#   geom_point()
# 
# #Now I need to rank the big and small lists. 
# 
# 
# pullM <- (selM %>%
#             filter(class != "NA") %>%
#             dplyr::select(Fly_ID, CS, diff, class)
#           )
# 
# write.csv(pullM, "../Outputs/zi192ef81_males.csv")
# 
# 
# pullF <- (selF %>%
#             filter(class != "NA") %>%
#             dplyr::select(Fly_ID, CS, diff, class)
# )
# 
# write.csv(pullF, "../Outputs/zi192ef81_females.csv")


###Because I have seen this before, I think it will just be easier to use the predicted sex rather than throw out all those flies (and of course check those) 



####################plotting out the distribution of selection##############


m.sel.plot <-ggplot(crossM, aes(x = CS, y = diff, col = class)) +
  geom_point(alpha = 0.5, size = 0.75) + 
  ylab("Shape Score") + 
  xlab("Centroid Size") + 
  scale_colour_manual(values=c('red', 'red', 'red', 'red', 'red', 'red', 'red', 'darkgrey')) +
  geom_vline(aes(xintercept = min(crossM[crossM$Fly_ID %in% SizeMbottom,]$CS))) + 
  geom_vline(aes(xintercept = max(crossM[crossM$Fly_ID %in% SizeMtop,]$CS))) +
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
  scale_colour_manual(values=c('red', 'red', 'red', 'red', 'red', 'red', 
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



png("../Figures/zi418ef43_bothSex_selection.png", 
    width =1700, height = 750, 
    units = "px", res = 300)

plot_grid(f.sel.plot, m.sel.plot, labels = c("Female", "Male"))

dev.off()









