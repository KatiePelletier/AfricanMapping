#a quick script to put all the populations into the same shape space. Parents and zi192ef81 are currently elsewhere because if my code wasn't chaotic, no one would know I wrote it. 

library(tidyverse)
library(geomorph)
library(abind) 

source("~/Dropbox/KatiePelletier/KP_geomorphwingfunctions.R")


pt.names <- c("X1.X", "X1.Y", "X2.X", "X2.Y", "X3.X", "X3.Y", "X4.X", "X4.Y","X5.X", "X5.Y", "X6.X", "X6.Y", "X7.X", "X7.Y", "X8.X", "X8.Y", "X9.X", "X9.Y", "X10.X", "X10.Y", "X11.X", "X11.Y", "X12.X", "X12.Y", "X13.X", "X13.Y", "X14.X", "X14.Y", "X15.X", "X15.Y")

#need the common mean too. 
common.mean <- read.csv("../Data/afMap_commonMean.csv")

############################

zi192ef43.15pt <- read.delim("../Data/zi192ef43_allwings_15pt_working.txt", header = FALSE)

#going directly into geomorph. Going to do the cleanup the names after. 
cord <- as.matrix(zi192ef43.15pt[,2:31])
shape <- arrayspecs(cord, 15, 2)

shape2 <- abind(common.mean, shape, along = 3)

#superimposition
cross.shape <- gpagen(shape2, ProcD = T, max.iter = 1)

#normal enough.  
hist(cross.shape$Csize)

plot(cross.shape, links = wing15.links)

plotOutliers(cross.shape$coords)

#this is not a wing. 
plot(cross.shape[["coords"]][,,829])

#2278 1234 1218, 829 are notreal 
crap.wings <- c(2278, 1234, 1218, 829)

#this looks much better
#plotOutliers(cross.shape[["coords"]][,,-crap.wings], inspect.outliers = TRUE)

crap <- as.data.frame(two.d.array(cross.shape$coords))
crap2 <- crap[-1,]
crap2$ID <- zi192ef43.15pt$V1

crap2[crap.wings,]$ID
# "KP_zi192ef43_c1208_plate25g_f_08"  "KP_Zi192Ef43_c1206_plate13e_f_04"
# "KP_Zi192Ef43_c1206_plate13c_f_08"  "KP_zi192xef43_c1204_plate8e_f_11"

crap2$CS <- cross.shape$Csize[-1]

hist(crap2$CS)

cross.15lmk <- crap2[-crap.wings,]

nrow(filter(cross.15lmk, CS < 2)) #75 wings. fuck this. 

#this is a problem. and I'm just getting rid of them at least for now. 
hist(cross.15lmk$CS)
cross.15lmk <- filter(cross.15lmk, CS > 2)
hist(cross.15lmk$CS)

#cleaning up the ID col 

garbage <- (cross.15lmk %>%
              separate(ID, into = c("perp", "line", "col", "plate", "sex", "fly")))

cross.15lmk <- (cross.15lmk %>%
                  separate(ID, into = c("perp", "line", "col", "plate", "sex", "fly")) %>%
                  dplyr::select(-c("perp", "col", "plate", "fly")))



#adding cross and log size to make data go together nice 
cross.15lmk$population <- "cross"
cross.15lmk$logCS <- log2(cross.15lmk$CS)
cross.15lmk$logCS_c <- cross.15lmk$logCS - mean(cross.15lmk$logCS)

str(cross.15lmk)

#not really needed becasue no stats but should match the other data. 
cross.15lmk[,c(31:32, 34)] <- lapply(cross.15lmk[,c(31:32, 34)], as.factor)
str(cross.15lmk)

#there is also a problem with the col names of the landmarks. 
colnames(cross.15lmk) <- c(pt.names, names(cross.15lmk[31:36]))

write.csv(cross.15lmk, "../Data/zi192ef43_15pt_superimposed.csv", row.names = FALSE, quote = FALSE)


write.csv(garbage, "../Data/zi192ef43_15pt_superimposed_withFly.csv", row.names = FALSE, quote = FALSE)

###############################################
zi192ef96.15pt <- read.delim("../Data/tempzi192ef96.txt", header = FALSE)

#problems
goodCrap <- zi192ef96.15pt


goodIDs <- zi192ef96.15pt[c(1:1003, 1082:1853),]

#going directly into geomorph. Going to do the cleanup the names after. 
cord <- as.matrix(goodCrap[,2:31])
shape <- arrayspecs(cord, 15, 2)

shape2 <- abind(common.mean, shape, along = 3)
#dont know why that wrote into the file. I'll remove. 
#1004-1082 all in set 5? Will redo. 
shape2[,,1081]
goodshape <- shape2[,,-c(1004:1081)]

dim(shape2)
dim(goodshape)

#superimposition
cross.shape <- gpagen(goodshape, ProcD = T, max.iter = 1)

#normal enough.  
hist(cross.shape$Csize)

#some problems still. 
plot(cross.shape, links = wing15.links)

plotOutliers(cross.shape$coords)

#this is not a wing. 
plot(cross.shape[["coords"]][,,255])

#944 are notreal 
crap.wings <- c(944, 255)

#this looks much better
#plotOutliers(cross.shape[["coords"]][,,-crap.wings], inspect.outliers = TRUE)

crap <- as.data.frame(two.d.array(cross.shape$coords))
crap2 <- crap[-1,]
crap2$ID <- goodIDs$V1

crap2[crap.wings,]$ID
#"KP_zi192xef96_c1114_plate17h_m_04"
#"KP_Zi192xEf96_c1114_plate3h_m_06"

crap2$CS <- cross.shape$Csize[-1]

hist(crap2$CS)

crap.wings

cross.15lmk <- crap2[-c(945, 256),]
hist(cross.15lmk$CS)

nrow(filter(cross.15lmk, CS < 2)) #75 wings. fuck this. 

#this is a problem. and I'm just getting rid of them at least for now. 
hist(cross.15lmk$CS)
cross.15lmk <- filter(cross.15lmk, CS > 2)
hist(cross.15lmk$CS)

#cleaning up the ID col 

garbage <- (cross.15lmk %>%
              separate(ID, into = c("perp", "line", "col", "plate", "sex", "fly")))

cross.15lmk <- (cross.15lmk %>%
                  separate(ID, into = c("perp", "line", "col", "plate", "sex", "fly")) )



#adding cross and log size to make data go together nice 
cross.15lmk$population <- "cross"
cross.15lmk$logCS <- log2(cross.15lmk$CS)
cross.15lmk$logCS_c <- cross.15lmk$logCS - mean(cross.15lmk$logCS)

str(cross.15lmk)

#not really needed becasue no stats but should match the other data. 
cross.15lmk[,c(31:32, 34)] <- lapply(cross.15lmk[,c(31:32, 34)], as.factor)
str(cross.15lmk)

#there is also a problem with the col names of the landmarks. 
colnames(cross.15lmk) <- c(pt.names, names(cross.15lmk[31:40]))

write.csv(cross.15lmk, "../Data/temp_zi192ef96_15pt_superimposed.csv", row.names = FALSE, quote = FALSE)


write.csv(garbage, "../Data/temp_zi192ef96_15pt_superimposed_withFly.csv", row.names = FALSE, quote = FALSE)





