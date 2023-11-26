#comparing the 15 point and full spline data to ask a couple questions about if this was actually an ok thing to do. Do I select the same flies? Do I think the vectors are similar. Using one set of data that I have done both ways. 

library(tidyverse)
library(cowplot)
library(geomorph)
library(abind) 

source("~/Dropbox/KatiePelletier/KP_geomorphwingfunctions.R")


#for easy projections and so I don't forget to divide. 
projFunction <- function(x, y) {
  scalarProj <- (x %*% y) / norm(y, type = "2")
  return(scalarProj)
}

#first reading in the 15lmk data. this is superimposed to the 15pt common mean elsewhere (I really hope that code is doing what I think it is otherwise I'm FUCKED)
pure.15pt <- read.csv("../Data/pureWings_15lmk_superimposed.csv")

str(pure.15pt)
pure.15pt[,32:34] <- lapply(pure.15pt[,32:34], as.factor)
str(pure.15pt)

#making the high and low vectors from (allometry corrected) landmarks and pcs. Want the other ones for fun, but really the PCs are what I used. 

#males first. 
pure.males.15pt <- filter(pure.15pt, sex == "M")

pure.males.fit.15pt <- data.frame(pure.males.15pt, lm(as.matrix(pure.males.15pt[,1:30])~ logCS_c, data = pure.males.15pt)$resid)

males.ha.15pt <- colMeans(pure.males.fit.15pt[pure.males.fit.15pt$population == "HA", 37:66])
males.la.15pt <- colMeans(pure.males.fit.15pt[pure.males.fit.15pt$population == "LA", 37:66])

males.diff.15pt <- males.ha.15pt - males.la.15pt

#That is unexpected. 
sqrt(sum(males.diff.15pt^2)) #0.005161079

###Bossman recomended also doing this without the size correction. 
males.ha.15ptlm <- colMeans(pure.males.fit.15pt[pure.males.fit.15pt$population == "HA", 1:30])
males.la.15ptlm <- colMeans(pure.males.fit.15pt[pure.males.fit.15pt$population == "LA", 1:30])

males.diff.15ptlm <- males.ha.15ptlm - males.la.15ptlm

#That is unexpected. 
sqrt(sum(males.diff.15ptlm^2)) #0.01929885

#I want to use ID plotting to look at the variation here between high and low. I think he might have made these functions by now... somewhere.
#commented out so I remember to swap this for the real function. 

line_colours <- c(rgb(0, 0, 1, 0.1),
                  rgb(1, 0, 1, 0.1))

point_colours <- c(rgb(0,0,0.8, 0.1),
                   rgb(0.8,0,0.8, 0.1))

wing.links.id <- c(1,7,12,13,14,15,11,10,9,8,
                6, 12, NA,
                6, 2, NA,
                8, 13, NA,
                10, 14, NA,
                9,3, NA,
                4,5, NA,
                5, 3, NA,
                5, 11)

index_x <- seq.int(from = 1, to = 29, by = 2)
index_y <- seq.int(from = 2, to = 30, by = 2)
line_index_x <- (wing.links.id*2 - 1) # landmark number, indexing column number in matrix for x coordinate
line_index_y <- wing.links.id*2

individuals <- as.matrix(pure.males.15pt[,1:30])

#png("../Figures/15pt_highLowMales_blur.png")
matplot(x = t(individuals[,line_index_x]),
        y = t(individuals[,line_index_y]),
        lty = 1, type = "l", lwd = 0.4,
        xlab = "", ylab = "", asp = 1,
        frame.plot = F, ann = F, axes = F,
        col = line_colours[pure.males.15pt$population])

matpoints(x = t(individuals[,index_x]),
          y = t(individuals[,index_y]),
          pch = 20, cex = 0.25,
          col = point_colours[pure.males.15pt$population])

wing_labels = c("highland", "lowland")

legend("bottom", ncol = 2, inset = -0.05,
       x.intersp = 0.5, y.intersp = 0.9,
       legend = wing_labels,
       text.col =  c(rgb(0,0,1, 1),
                     rgb(1,0,1, 1)),
       bty = "n")

#dev.off()

#Now with females. 
pure.females.15pt <- filter(pure.15pt, sex == "F")

pure.females.fit.15pt <- data.frame(pure.females.15pt, lm(as.matrix(pure.females.15pt[,1:30])~ logCS_c, data = pure.females.15pt)$resid)

females.ha.15pt <- colMeans(pure.females.fit.15pt[pure.females.fit.15pt$population == "HA", 37:66])
females.la.15pt <- colMeans(pure.females.fit.15pt[pure.females.fit.15pt$population == "LA", 37:66])

females.diff.15pt <- females.ha.15pt - females.la.15pt

#That is unexpected. 
sqrt(sum(females.diff.15pt^2)) #0.005106257


females.ha.15ptlm <- colMeans(pure.females.fit.15pt[pure.females.fit.15pt$population == "HA", 1:30])
females.la.15ptlm <- colMeans(pure.females.fit.15pt[pure.females.fit.15pt$population == "LA", 1:30])

females.diff.15ptlm <- females.ha.15ptlm - females.la.15ptlm

#That is unexpected. 
sqrt(sum(females.diff.15ptlm^2)) #0.01521056


#Now to compare males and female effects because those seem to be diff. 
cor(females.diff.15pt, males.diff.15pt) #0.9203572
cor(females.diff.15ptlm, males.diff.15ptlm) #0.9599264

#now for the PC based approach. 

pure.pca.m.15pt <- prcomp(pure.males.15pt[c(1:30, 35)])

summary(pure.pca.m.15pt)

pure.wings.pca.m.15pt <- data.frame(pure.males.15pt, pure.pca.m.15pt$x)

# all good. 
ggplot(pure.wings.pca.m.15pt, aes(x = logCS, y = PC1, col = line)) + 
  geom_point(alpha = 0.5)

#More interesting plots? 
#I've seen this before. Later I can fit the actual models and get the effect sizes back. Its not the best and I remember this now. 

#double check that the red and black works with other figures I've made. 
pc23.15pt.m <- ggplot(pure.wings.pca.m.15pt, aes(x = PC2, y = PC3, col = line)) + 
  geom_point(alpha = 0.6) + 
  scale_color_manual(values = c("gray80", "gray38", "black", "lightcoral", "red", "darkred")) + 
  theme_classic() + 
  theme(text = element_text(size=12),
        axis.text.x= element_text(size=10), 
        axis.text.y= element_text(size=10), 
        legend.position = "none")

pc45.15pt.m <- ggplot(pure.wings.pca.m.15pt, aes(x = PC4, y = PC5, col = line)) + 
  geom_point(alpha = 0.6) + 
  scale_color_manual(values = c("gray80", "gray38", "black", "lightcoral", "red", "darkred")) + 
  theme_classic() + 
  theme(text = element_text(size=12),
        axis.text.x= element_text(size=10), 
        axis.text.y= element_text(size=10))

pc45.15pt.m.noleg <- ggplot(pure.wings.pca.m.15pt, aes(x = PC4, y = PC5, col = line)) + 
  geom_point(alpha = 0.6) + 
  scale_color_manual(values = c("gray80", "gray38", "black", "lightcoral", "red", "darkred")) + 
  theme_classic() + 
  theme(text = element_text(size=12),
        axis.text.x= element_text(size=10), 
        axis.text.y= element_text(size=10), 
        legend.position = "none")

# png("../Figures/pca15pt_m.png", width = 3000, height = 1200, units = "px", res = 300)
# plot_grid(pc23.15pt.m, pc45.15pt.m, rel_widths = c(0.9, 1))
# dev.off()

ha.m.pcs.15pt <- colMeans(pure.wings.pca.m.15pt[pure.wings.pca.m.15pt$population == "HA",38:64])
la.m.pcs.15pt <- colMeans(pure.wings.pca.m.15pt[pure.wings.pca.m.15pt$population == "LA",38:64])

male.diff.pcs.15pt <- ha.m.pcs.15pt - la.m.pcs.15pt

sqrt(sum(male.diff.pcs.15pt^2)) #0.005152562

#and the females 
pure.pca.f.15pt <- prcomp(pure.females.15pt[c(1:30, 35)])

summary(pure.pca.f.15pt)

pure.wings.pca.f.15pt <- data.frame(pure.females.15pt, pure.pca.f.15pt$x)

# all good. 
ggplot(pure.wings.pca.f.15pt, aes(x = logCS, y = PC1, col = line)) + 
  geom_point(alpha = 0.5)

#More interesting plots? 
#Not exactly the same at the males. 
#double check that the red and black works with other figures I've made. 
pc23.15pt.f <- ggplot(pure.wings.pca.f.15pt, aes(x = PC2, y = PC3, col = line)) + 
  geom_point(alpha = 0.6) + 
  scale_color_manual(values = c("gray80", "gray38", "black", "lightcoral", "red", "darkred")) + 
  theme_classic() + 
  theme(text = element_text(size=12),
        axis.text.x= element_text(size=10), 
        axis.text.y= element_text(size=10), 
        legend.position = "none")

pc45.15pt.f <- ggplot(pure.wings.pca.f.15pt, aes(x = PC4, y = PC5, col = line)) + 
  geom_point(alpha = 0.6) + 
  scale_color_manual(values = c("gray80", "gray38", "black", "lightcoral", "red", "darkred")) + 
  theme_classic() + 
  theme(text = element_text(size=12),
        axis.text.x= element_text(size=10), 
        axis.text.y= element_text(size=10))

# png("../Figures/pca15pt_f.png", width = 3000, height = 1200, units = "px", res = 300)
# plot_grid(pc23.15pt.f, pc45.15pt.f, rel_widths = c(0.9, 1))
# dev.off()

ha.f.pcs.15pt <- colMeans(pure.wings.pca.f.15pt[pure.wings.pca.f.15pt$population == "HA",38:64])
la.f.pcs.15pt <- colMeans(pure.wings.pca.f.15pt[pure.wings.pca.f.15pt$population == "LA",38:64])

female.diff.pcs.15pt <- ha.f.pcs.15pt - la.f.pcs.15pt

sqrt(sum(female.diff.pcs.15pt^2)) #0.005152562

#################a quick break to compare the males and female matrices. For fun?####################################

#first 5 because it prints nice
#these actually look sort of similar. No obvious "flipping"  
cor(pure.pca.f.15pt$rotation[,1:5], pure.pca.m.15pt$rotation[,1:5])


################back the the original question#################
#I want to look at the distances between the populations 
#doing this one sex at a time for ease and also to match everything else. 
cord <- as.matrix(pure.males.15pt[,1:30])
shape <- arrayspecs(cord, 15, 2)

gdf.pure.males15 <- geomorph.data.frame(shape = shape,
                                 CS = pure.males.15pt$CS, 
                                line = pure.males.15pt$line, 
                                population = pure.males.15pt$population, 
                                logCS = pure.males.15pt$logCS, 
                                logCS_c = pure.males.15pt$logCS_c)

#Now I want to ask how similar the populations are to one another 
#nesting lines within the popualtion to account for relatedness
pure15.mod.males <- procD.lm(shape ~ 1 + logCS_c*population + 
                        population:line + logCS_c:population:line,
                       data = gdf.pure.males15 )

summary(pure15.mod.males)

pure15.null.males <- procD.lm(shape ~ 1 + logCS_c*population
                              + logCS_c:population:line, 
                        data = gdf.pure.males15)

#Does this actually work with line being nested or am I asking the computer to do something crazy here. Biologically what I did makes senese to me? But I don't think you can get these numbers from the model?
pure15.pair.males <- pairwise(pure15.mod.males, pure15.null.males, pure.males.15pt$line)

#distances by default
summary(pure15.pair.males)

read.out.pure15.pair <- data.frame(summary(pure15.pair.males)[["summary.table"]])

#write.csv(read.out.pure15.pair, "../Data/15pointLines.pairdistances.csv")

#do I recover a similar distance doing it this was as I was above?
pure15.null.Mpop <- procD.lm(shape ~ 1 + logCS_c+ 
                               logCS_c:population +  
                               population:line + 
                               logCS_c:population:line,
                             data = gdf.pure.males15 )


pure15.males.parirPops <-pairwise(pure15.mod.males, pure15.null.Mpop, pure.males.15pt$population)

#bigger by the estimate but with big CI 
summary(pure15.males.parirPops)

#this doesn't explain it. 
summary(pure15.males.parirPops, test.type = "var")

###################################################################

#I need to take the code and actually do the selection with it. 
#reading in the 15 point data and superimposing to the common mean so its in the same shapespace as the parent data. 
zi192ef81.15pt <- read.delim("../Data/zi192ef81_15point.txt", header = FALSE)

#need the common mean too. 
common.mean <- read.csv("../Data/afMap_commonMean.csv")

#going directly into geomorph. Going to do the cleanup the names after. 

#I checked all these to see if there are any crazy numbers. but I don't see it here. these all look like not crazy numbers. 
#hist(zi192ef81.15pt[,31])

cord <- as.matrix(zi192ef81.15pt[,2:31])
shape <- arrayspecs(cord, 15, 2)

#problematic 
shape[,,2193]
better.shape <- shape[,,-2193]

shape2 <- abind(common.mean, better.shape, along = 3)

#superimposition
cross.shape <- gpagen(shape2, ProcD = T, max.iter = 1)

#normal enough.  
hist(cross.shape$Csize)

plot(cross.shape, links = wing15.links)

plotOutliers(cross.shape$coords)

#this is not a wing. 
plot(cross.shape[["coords"]][,,563])

#607, 38 and 539 have some pretty funky positions as well. I will drop those out. 
#plotOutliers(cross.shape[["coords"]][,,-563], inspect.outliers = TRUE)

crap.wings <- c(607, 38, 539, 563)

crap <- as.data.frame(two.d.array(cross.shape$coords))
crap2 <- crap[-1,]
crap2$ID <- zi192ef81.15pt$V1[-2193]

crap2[crap.wings,]$ID
# "IED_Zi192Ef81_c1217_plate4h_m_02"  "IED_Zi192Ef81_c1217_plate14f_m_04"
# "IED_Zi192Ef81_c1217_plate4a_f_07"  "IED_Zi192Ef81_c1217_plate4c_f_09"

crap2$CS <- cross.shape$Csize[-1]

cross.15lmk <- crap2[-crap.wings,]

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
colnames(cross.15lmk) <- c(names(pure.15pt[1:30]), names(cross.15lmk[31:36]))


#I want to write out the superimposed cross
#write.csv(cross.15lmk, "../Data/zi192ef81_15pt_superimposed.csv", row.names = FALSE, quote = FALSE)

#write.csv(garbage, "../Data/zi192ef81_15pt_superimposed_withFly.csv", row.names = FALSE, quote = FALSE)

#Now to do the projections.
#males first. 
cross.15lmk.males <- filter(cross.15lmk, sex == "m")

#need to do the PCA all together.
allMales.15lmk <- rbind(cross.15lmk.males, pure.males.15pt)

allMales.15lmk.pcs <- data.frame(allMales.15lmk, prcomp(allMales.15lmk[,c(1:30, 35)])$x)

#checking
ggplot(allMales.15lmk.pcs, aes(x = logCS, y = PC1)) + 
  geom_point()

#now to make the projection vectors. 
pure.forProj.males15 <- filter(allMales.15lmk.pcs, population != "cross")

ha.m.proj.15pt <- colMeans(pure.forProj.males15[pure.forProj.males15$population == "HA",38:64])
la.m.proj.15pt <- colMeans(pure.forProj.males15[pure.forProj.males15$population == "LA",38:64])

male.diff.proj.15pt <- ha.m.proj.15pt - la.m.proj.15pt

cross.forProj.males15 <- filter(allMales.15lmk.pcs, population == "cross")

#now for the actual projection. 
#I want to add the fly ID tag back here because I want to compare these lists later This is a wild way to do this but I am working fast to get something to look at later. 
cross.forProj.males15$Fly_ID <- paste(garbage[garbage$sex == "m",]$plate, garbage[garbage$sex == "m",]$fly, sep = "_")


cross.forProj.males15$shapeScore <- projFunction(as.matrix(cross.forProj.males15[,38:64]), male.diff.proj.15pt)

#that is a distribution. One weird guy out there but thats a future probelm. I have the vec now. 
hist(cross.forProj.males15$shapeScore)

top50sizeMales.15lmk  <- (cross.forProj.males15 %>% top_n(50, CS))
bottom50sizeMales.15lmk  <- (cross.forProj.males15 %>% top_n(-50, CS))


top50shapeMales.15lmk  <- (cross.forProj.males15 %>% top_n(50, shapeScore))
bottom50shapeMales.15lmk  <- (cross.forProj.males15 %>% top_n(-50, shapeScore))

#now for the full splines. I am going to take this data that has already been cleaned up elsewhere and use that becasue I don't want to reinvent the wheel. s

load("../Data/zi192ef81clean.rda")

#renaming for clarity. 
wings.spline <- wings
str(wings.spline)

#This is where the important info is stored. 
levels(wings.spline$line)

pure.spline <- filter(wings.spline, line != "Zi192Ef81")
cross.spline <- filter(wings.spline, line == "Zi192Ef81")

#now grabbing the effect vectors both ways to compare to the 15point ones. 
#adding population first 
pure.spline$population <- ifelse(grepl("Z", pure.spline$line), "LA", "HA")
pure.spline$logCS <- log(pure.spline$CS)
pure.spline$logCS_c <- pure.spline$logCS - mean(pure.spline$logCS)

pure.males.spline<- filter(pure.spline, sex == "M")

pure.males.fit.spline <- data.frame(pure.males.spline, lm(as.matrix(pure.males.spline[,6:101]) ~ logCS_c, data = pure.males.spline)$resid)

males.ha.spline <- colMeans(pure.males.fit.spline[pure.males.fit.spline$population == "HA", 107:202])
males.la.spline <- colMeans(pure.males.fit.spline[pure.males.fit.spline$population == "LA", 107:202])

males.diff.spline <- males.ha.spline - males.la.spline

#with allometry in there. 
males.ha.splinelm <- colMeans(pure.males.fit.spline[pure.males.fit.spline$population == "HA", 6:101])
males.la.splinelm <- colMeans(pure.males.fit.spline[pure.males.fit.spline$population == "LA", 6:101])

males.diff.splinelm <- males.ha.splinelm - males.la.splinelm

#That is unexpected. 
#15 point: 0.005161079
#How did we not catch this before? baby grad student Katie was really bad at her job sometimes. 
sqrt(sum(males.diff.spline^2)) #0.004746787

#Now with females. 
pure.females.spline <- filter(pure.spline, sex == "F")

pure.females.fit.spline <- data.frame(pure.females.spline, lm(as.matrix(pure.females.spline[,6:101])~ logCS_c, data = pure.females.spline)$resid)

females.ha.spline <- colMeans(pure.females.fit.spline[pure.females.fit.spline$population == "HA", 107:202])
females.la.spline <- colMeans(pure.females.fit.spline[pure.females.fit.spline$population == "LA", 107:202])

females.diff.spline <- females.ha.spline - females.la.spline

#landmarks
females.ha.splinelm <- colMeans(pure.females.fit.spline[pure.females.fit.spline$population == "HA", 6:101])
females.la.splinelm <- colMeans(pure.females.fit.spline[pure.females.fit.spline$population == "LA", 6:101])

females.diff.splinelm <- females.ha.splinelm - females.la.splinelm

#15point:0.005106257
sqrt(sum(females.diff.spline^2)) #0.005328175

#Little less then before (with 15). but still. 
cor(females.diff.spline, males.diff.spline) #0.8573106

#now for the PC based approach. 
pure.pca.m.spline <- prcomp(pure.males.spline[c(7:101, 105)])

summary(pure.pca.m.spline)

pure.wings.pca.m.spline <- data.frame(pure.males.spline, pure.pca.m.spline$x[,1:59])

# all good. 
ggplot(pure.wings.pca.m.spline, aes(x = logCS, y = PC1, col = line)) + 
  geom_point(alpha = 0.5)

#More interesting plots? 
#looks exactly the same as the other plot (there is a flip in the direction of PC2 but nothing exciting)
pc23.spline.m <- ggplot(pure.wings.pca.m.spline, aes(x = PC2, y = PC3, col = line)) + 
  geom_point(alpha = 0.6) + 
  scale_color_manual(values = c("gray80", "gray38", "black", "lightcoral", "red", "darkred")) + 
  theme_classic() + 
  theme(text = element_text(size=12),
        axis.text.x= element_text(size=10), 
        axis.text.y= element_text(size=10), 
        legend.position = "none")

pc45.spline.m <- ggplot(pure.wings.pca.m.spline, aes(x = PC4, y = PC5, col = line)) + 
  geom_point(alpha = 0.6) + 
  scale_color_manual(values = c("gray80", "gray38", "black", "lightcoral", "red", "darkred")) + 
  theme_classic() + 
  theme(text = element_text(size=12),
        axis.text.x= element_text(size=10), 
        axis.text.y= element_text(size=10))

pc45.spline.m.noleg <- ggplot(pure.wings.pca.m.spline, aes(x = PC4, y = PC5, col = line)) + 
  geom_point(alpha = 0.6) + 
  scale_color_manual(values = c("gray80", "gray38", "black", "lightcoral", "red", "darkred")) + 
  theme_classic() + 
  theme(text = element_text(size=12),
        axis.text.x= element_text(size=10), 
        axis.text.y= element_text(size=10), 
        legend.position = "none")

# png("../Figures/pcaSpline_m.png", width = 3000, height = 1200, units = "px", res = 300)
# plot_grid(pc23.spline.m, pc45.spline.m, rel_widths = c(0.9, 1))
# dev.off()

ha.m.pcs.spline <- colMeans(pure.wings.pca.m.spline[pure.wings.pca.m.spline$population == "HA",108:164])
la.m.pcs.spline <- colMeans(pure.wings.pca.m.spline[pure.wings.pca.m.spline$population == "LA",108:164])

male.diff.pcs.spline <- ha.m.pcs.spline - la.m.pcs.spline

#15pt:  0.005152562
sqrt(sum(male.diff.pcs.spline^2)) #0.004719508

#and the females 
pure.pca.f.spline <- prcomp(pure.females.spline[c(7:101, 105)])

summary(pure.pca.f.spline)

pure.wings.pca.f.spline <- data.frame(pure.females.spline, pure.pca.f.spline$x[,1:59])

# all good. 
ggplot(pure.wings.pca.f.spline, aes(x = logCS, y = PC1, col = line)) + 
  geom_point(alpha = 0.5)

#More interesting plots? 
#more diff sort of here? I think because that one wing is here and gets dropped in the 15point. 
pc23.spline.f <- ggplot(pure.wings.pca.f.spline, aes(x = PC2, y = PC3, col = line)) + 
  geom_point(alpha = 0.6) + 
  scale_color_manual(values = c("gray80", "gray38", "black", "lightcoral", "red", "darkred")) + 
  theme_classic() + 
  theme(text = element_text(size=12),
        axis.text.x= element_text(size=10), 
        axis.text.y= element_text(size=10), 
        legend.position = "none")

pc45.spline.f <- ggplot(pure.wings.pca.f.spline, aes(x = PC4, y = PC5, col = line)) + 
  geom_point(alpha = 0.6) + 
  scale_color_manual(values = c("gray80", "gray38", "black", "lightcoral", "red", "darkred")) + 
  theme_classic() + 
  theme(text = element_text(size=12),
        axis.text.x= element_text(size=10), 
        axis.text.y= element_text(size=10))

pc45.spline.f.noleg <- ggplot(pure.wings.pca.f.spline, aes(x = PC4, y = PC5, col = line)) + 
  geom_point(alpha = 0.6) + 
  scale_color_manual(values = c("gray80", "gray38", "black", "lightcoral", "red", "darkred")) + 
  theme_classic() + 
  theme(text = element_text(size=12),
        axis.text.x= element_text(size=10), 
        axis.text.y= element_text(size=10), 
        legend.position = "none")

# png("../Figures/pcaSpline_f.png", width = 3000, height = 1200, units = "px", res = 300)
#plot_grid(pc23.spline.f, pc45.spline.f, rel_widths = c(0.9, 1))
# dev.off()

pcas.leg <- get_legend(pc45.spline.m)

bothSexPCAplot <- plot_grid(pc23.15pt.m, pc45.15pt.m.noleg, pc23.spline.f,
                             pc45.spline.f.noleg, nrow = 2, ncol = 2, 
                             labels = c("male", "", "female", ""))

bothWaysPCAplot <- plot_grid(pc23.15pt.m, pc45.15pt.m.noleg, pc23.spline.m,
                            pc45.spline.m.noleg, nrow = 2, ncol = 2, 
                            labels = c("15 point", "", "spline", ""), 
                            label_size = 12,
                            label_x = 0.2, label_y = 1)

#png("../Figures/bothSexPCAplot.png", width =2500, height = 2200, units = "px",res = 300)
#plot_grid(bothSexPCAplot, pcas.leg, rel_widths = c(1,0.25))
#dev.off()

png("../Figures/bothWaysMalesPCAplot.png", width =2500, height = 2500, units = "px",res = 300)
plot_grid(bothWaysPCAplot, pcas.leg, rel_widths = c(1,0.25))
dev.off()

ha.f.pcs.spline <- colMeans(pure.wings.pca.f.spline[pure.wings.pca.f.spline$population == "HA",108:164])
la.f.pcs.spline <- colMeans(pure.wings.pca.f.spline[pure.wings.pca.f.spline$population == "LA",108:164])

female.diff.pcs.spline <- ha.f.pcs.spline - la.f.pcs.spline

#15point: 0.005152562
sqrt(sum(female.diff.pcs.spline^2)) #0.005288473

#################a quick break to compare the males and female matrices. For fun?####################################

#first 5 because it prints nice
#look pretty similar to before. Not really that interesting to look at tho. 
cor(pure.pca.f.spline$rotation[,1:5], pure.pca.m.spline$rotation[,1:5])

################back the the original question#################
#I want to look at the distances between the populations 
#doing this one sex at a time for ease and also to match everything else. 
cord <- as.matrix(pure.males.spline[,6:101])
shape <- arrayspecs(cord, 48, 2)

gdf.pure.males.spline <- geomorph.data.frame(shape = shape,
                                        CS = pure.males.spline$CS, 
                                        line = pure.males.spline$line, 
                                        population = pure.males.spline$population, 
                                        logCS = pure.males.spline$logCS, 
                                        logCS_c = pure.males.spline$logCS_c)

#Now I want to ask how similar the populations are to one another 
#nesting lines within the popualtion to account for relatedness
spline.mod.males <- procD.lm(shape ~ 1 + logCS_c*population + 
                               population:line + logCS_c:population:line,
                             data = gdf.pure.males.spline )

summary(spline.mod.males)

spline.null.males <- procD.lm(shape ~ 1 + logCS_c*population
                              + logCS_c:population:line, 
                              data = gdf.pure.males.spline)

#Does this actually work with line being nested or am I asking the computer to do something crazy here. Biologically what I did makes senese to me? But I don't think you can get these numbers from the model?
spline.pair.males <- pairwise(spline.mod.males, spline.null.males, pure.males.spline$line)

#distances by default
summary(spline.pair.males)

read.out.spline.pair <- data.frame(summary(spline.pair.males)[["summary.table"]])

#write.csv(read.out.spline.pair, "../Data/splineLines.pairdistances_male.csv")

#I want to project the cross wings, and compare shapescores with the 15 and spline generated vectors. 

#Now to do the projections.
#need to do the PCA all together.
allMales.spline <- filter(wings.spline, sex == "M")
allMales.spline$logCS <- log(allMales.spline$CS)

allMales.spline.pcs <- data.frame(allMales.spline, prcomp(allMales.spline[,c(6:101, 104)])$x[,1:58])

#checking
ggplot(allMales.spline.pcs, aes(x = logCS, y = PC1)) + 
  geom_point()

#now to make the projection vectors. 
pure.forProj.malesSpline <- filter(allMales.spline.pcs, line != "Zi192Ef81")

pure.forProj.malesSpline$population <- ifelse(grepl("Z", pure.forProj.malesSpline$line), "LA", "HA")

ha.m.proj.Spline <- colMeans(pure.forProj.malesSpline[pure.forProj.malesSpline$population == "HA",106:162])
la.m.proj.Spline <- colMeans(pure.forProj.malesSpline[pure.forProj.malesSpline$population == "LA",106:162])

male.diff.proj.Spline <- ha.m.proj.Spline - la.m.proj.Spline

cross.forProj.malesSpline <- filter(allMales.spline.pcs, line == "Zi192Ef81")

#now for the actual projection. 

cross.forProj.malesSpline$shapeScore <- projFunction(as.matrix(cross.forProj.malesSpline[,106:162]), male.diff.proj.Spline)

#this is a diffrent distribution but that means nothing really
hist(cross.forProj.malesSpline$shapeScore)

top50sizeMales.spline  <- (cross.forProj.malesSpline %>% top_n(50, CS))
bottom50sizeMales.spline  <- (cross.forProj.malesSpline %>% top_n(-50, CS))


top50shapeMales.spline  <- (cross.forProj.malesSpline %>% top_n(50, shapeScore))
bottom50shapeMales.spline  <- (cross.forProj.malesSpline %>% top_n(-50, shapeScore))

#Now similar? 
#thankfully similar. Some diffs may be because some flies are not in both. 
sum(top50sizeMales.spline$Fly_ID %in% top50sizeMales.15lmk$Fly_ID)
sum(bottom50sizeMales.spline$Fly_ID %in% bottom50sizeMales.15lmk$Fly_ID)

#less awesome actually. fuck. fuck. 
sum(top50shapeMales.spline$Fly_ID %in% top50shapeMales.15lmk$Fly_ID)
sum(bottom50shapeMales.spline$Fly_ID %in% bottom50shapeMales.15lmk$Fly_ID)

#I want to just look at all the shape scores. First I need to get the IDs in both. 
inboth <- cross.forProj.males15[cross.forProj.males15$Fly_ID %in% cross.forProj.malesSpline$Fly_ID,]

inboth$splineScore <- cross.forProj.malesSpline$shapeScore[c(which(cross.forProj.males15$Fly_ID %in% cross.forProj.malesSpline$Fly_ID))]

#well that sucks. 
ggplot(inboth, aes(x = shapeScore, y = splineScore)) + 
  geom_point() + 
  geom_smooth(col = "red") + 
  xlab("Shape Score from 15pt") + 
  ylab("Shape Score from spline")


######################
#now using the landmarks alone for the projection (not corrected for allometry)

cross.forProj.males15$shapeScore2 <- projFunction(as.matrix(cross.forProj.males15[,1:30]), males.diff.15ptlm)

cross.forProj.malesSpline$shapeScore2 <- projFunction(as.matrix(cross.forProj.malesSpline[,6:101]), males.diff.splinelm)


#I want to just look at all the shape scores. First I need to get the IDs in both. 
inboth <- cross.forProj.males15[cross.forProj.males15$Fly_ID %in% cross.forProj.malesSpline$Fly_ID,]

inboth$splineScore2 <- cross.forProj.malesSpline$shapeScore2[c(which(cross.forProj.males15$Fly_ID %in% cross.forProj.malesSpline$Fly_ID))]

#Not making me super hopeful. fuck. 
png("../Figures/shapeScoreCor_usingLandmarks_noCorrection.png")
ggplot(inboth, aes(x = shapeScore2, y = splineScore2)) + 
  geom_point() + 
  geom_smooth(col = "red") + 
  xlab("Shape Score from 15pt") + 
  ylab("Shape Score from spline")

dev.off()

#Ok so do I at least select the same? 
top50shapeMales.15lmklm  <- (cross.forProj.males15 %>% top_n(50, shapeScore2))
bottom50shapeMales.15lmklm <- (cross.forProj.males15 %>% top_n(-50, shapeScore2))

top50shapeMales.splinelm  <- (cross.forProj.malesSpline %>% top_n(50, shapeScore2))
bottom50shapeMales.splinelm  <- (cross.forProj.malesSpline %>% top_n(-50, shapeScore2))


#Ok. 30 and 32. Not the worst. 
sum(top50shapeMales.splinelm$Fly_ID %in% top50shapeMales.15lmklm$Fly_ID)
sum(bottom50shapeMales.splinelm$Fly_ID %in% bottom50shapeMales.15lmklm$Fly_ID)


####drawing some pictures for ID 
#going to do this with a blur 
top50shapeMales.15lmk$group <- "top"
bottom50shapeMales.15lmk$group <- "bottom"

selected.15shape <- rbind(top50shapeMales.15lmk, bottom50shapeMales.15lmk)

individuals <- as.matrix(selected.15shape[,1:30])

#png("../Figures/upDown15ptBlur_nofit.png")
# matplot(x = t(individuals[,line_index_x]),
#         y = t(individuals[,line_index_y]),
#         lty = 1, type = "l", lwd = 0.4,
#         xlab = "", ylab = "", asp = 1,
#         frame.plot = F, ann = F, axes = F,
#         col = line_colours[selected.15shape$group])
# 
# matpoints(x = t(individuals[,index_x]),
#           y = t(individuals[,index_y]),
#           pch = 20, cex = 0.25,
#           col = point_colours[selected.15shape$group])
# 
# wing_labels = c("top", "bottom")
# 
# legend("bottom", ncol = 2, inset = -0.05,
#        x.intersp = 0.5, y.intersp = 0.9,
#        legend = wing_labels,
#        text.col =  c(rgb(0,0,1, 1),
#                      rgb(1,0,1, 1)),
#        bty = "n")

#dev.off()

#adding in 

#that wasn't working so I'm going to just play with the geomorph plots for now. Blurs can come later. 
source("~/Dropbox/KatiePelletier/KP_geomorphwingfunctions.R")

#Why is this so small? None of these capture the posterior cross vein shift. 
cord1 <- as.matrix(selected.15shape[,1:30])
shape1 <- arrayspecs(cord1, 15, 2)
selected.15shape.gdf <- geomorph.data.frame(shape = shape1,
                                  group = selected.15shape$group)


new.coords1 <- coords.subset(selected.15shape.gdf$shape, group = selected.15shape$group)
names(new.coords) # see the list levels
# group shape means

pool_mean1 <- lapply(new.coords1, mshape)

png("../Figures/ShapePools_15pt_fromPC.png")
plotRefToTarget(pool_mean1[["top"]], 
                pool_mean1[["bottom"]], 
                links = fifteen.links, method = "points", mag = 2, 
                gridPars=wing.spec )
dev.off()

top50shapeMales.15lmklm$group <- "top"
bottom50shapeMales.15lmklm$group <- "bottom"

selected.15shapelm <- rbind(top50shapeMales.15lmklm, bottom50shapeMales.15lmklm)

cord2 <- as.matrix(selected.15shapelm[,1:30])
shape2 <- arrayspecs(cord2, 15, 2)
selected.15shapelm.gdf <- geomorph.data.frame(shape = shape2,
                                            group = selected.15shapelm$group)

new.coords2 <- coords.subset(selected.15shapelm.gdf$shape, group = selected.15shapelm$group)
names(new.coords2) # see the list levels
# group shape means

pool_mean2 <- lapply(new.coords2, mshape)



png("../Figures/ShapePools_15pt_nofit.png")
plotRefToTarget(pool_mean2[["top"]], 
                pool_mean2[["bottom"]], 
                links = fifteen.links, method = "points", mag = 2, 
                gridPars=wing.spec )
dev.off()
#######

top50shapeMales.splinelm$group <- "top"
bottom50shapeMales.splinelm$group <- "bottom"

selected.splinelm <- rbind(top50shapeMales.splinelm, bottom50shapeMales.splinelm)

cord3 <- as.matrix(selected.splinelm[,6:101])
shape3 <- arrayspecs(cord3, 48, 2)
selected.splinelm.gdf <- geomorph.data.frame(shape = shape3,
                                              group = selected.splinelm$group)

new.coords3 <- coords.subset(selected.splinelm.gdf$shape, group = selected.splinelm$group)
names(new.coords3) # see the list levels
# group shape means
pool_mean3<- lapply(new.coords3, mshape)

png("../Figures/ShapePools_spline_nofit.png")
plotRefToTarget(pool_mean3[["top"]], 
                pool_mean3[["bottom"]], 
                links = wing.links, method = "points", mag = 2, 
                gridPars=wing.spec )
dev.off()

#######

top50shapeMales.spline$group <- "top"
bottom50shapeMales.spline$group <- "bottom"

selected.spline <- rbind(top50shapeMales.spline, bottom50shapeMales.spline)

cord4 <- as.matrix(selected.spline[,6:101])
shape4 <- arrayspecs(cord4, 48, 2)
selected.spline.gdf <- geomorph.data.frame(shape = shape4,
                                             group = selected.spline$group)

new.coords4 <- coords.subset(selected.spline.gdf$shape, group = selected.spline$group)
names(new.coords4) # see the list levels
# group shape means
pool_mean4<- lapply(new.coords4, mshape)

#this wing is so much smaller?
png("../Figures/ShapePools_spline_fromPC.png")
plotRefToTarget(pool_mean4[["top"]], 
                pool_mean4[["bottom"]], 
                links = wing.links, method = "points", mag = 2, 
                gridPars=wing.spec )
dev.off()

###################################################################
#Making figure for paper to compare the shape changes
##################################################################

library(cowplot)
library(magick)

nomod.15.raw <- image_read("../Figures/ShapePools_15pt_nofit.png")
nomod.15.crop <- image_trim(nomod.15.raw)
nomod.15 <- ggdraw() + draw_image(nomod.15.crop)

mod.15.raw <- image_read("../Figures/ShapePools_15pt_fromPC.png")
mod.15.crop <- image_trim(mod.15.raw)
mod.15 <- ggdraw() + draw_image(mod.15.crop)

#need to flip these as well. 
nomod.spline.raw <- image_read("../Figures/ShapePools_spline_nofit.png")
nomod.spline.crop <- image_trim(nomod.spline.raw)
nomod.spline.flop <- image_flop(nomod.spline.crop)
nomod.spline <- ggdraw() + draw_image(nomod.spline.flop)

mod.spline.raw <- image_read("../Figures/ShapePools_spline_fromPC.png")
mod.spline.crop <- image_trim(mod.spline.raw)
mod.spline.flop <- image_flop(mod.spline.crop)
mod.spline <- ggdraw() + draw_image(mod.spline.flop)

png("../Figures/15vsSplineWireframes.png", width =2000, height = 1500, units = "px",res = 300)
plot_grid(nomod.15, mod.15, nomod.spline, mod.spline, nrow = 2, ncol = 2, 
          labels = c("no allometry correction", "allometry correction"))

dev.off()





