library(tidyverse)
#library(ggbeeswarm)

########################
source( "~/Dropbox/DworkinLabSharedMaterial/scripts/WRP_FUNCTIONS.R" )
source( "~/Dropbox/DworkinLabSharedMaterial/scripts/WINGPLOTSOURCE.R" )
###########################

#Plotting an intro figure with the distribution of shapes between high and low alt populations for both size and shape. For this figure I want to use full splines. 

load("../Data/zi192ef96_allwings.rda")

wings.parents <- filter(wings, Genotype != "zi192xef96")

parent.size <- ggplot(wings.parents, aes(x = Genotype, y = CS, fill = sex)) + 
  geom_violin(alpha = 0.6) + 
  xlab("") + 
  ylab("Size") + 
  scale_fill_manual(values = c("grey", "black"), labels = c("Female", "Male")) +
  theme_classic()
  
#I also want to plot the variation between high and low altitude popualtions. 
#Right now this plots everything and is pretty low.

wings.parents2 <- filter(wings.parents, sex == "m")
wings.parents3 <- filter(wings.parents, sex == "f")


WingBlur3(as.matrix(wings.parents[,4:99]), groups = TRUE, droplevels(wings.parents$Alt), blur_transparency = 0.02)

#males only. no correction. 
png("../Figures/malesPureBlur_noCorrection.png")
WingBlur3(as.matrix(wings.parents2[,4:99]), groups = TRUE, droplevels(wings.parents2$Alt), blur_transparency = 0.02)
dev.off()

#What if I correct for allometry first?
# corrected.lm <- lm(as.matrix(wings.parents2[,4:99])~ logCS, data =wings.parents2)$resid
# 
# #need to add a mean shape back for plotting
# mean.shape <- colMeans(wings.parents2[,4:99])
# 
# plotting.lm <- corrected.lm + rep(mean.shape, each = nrow(corrected.lm))
# 
# #this looks worse. I think I should fit within group. 
# WingBlur3(plotting.lm, groups = TRUE, droplevels(wings.parents2$Alt), blur_transparency = 0.02)


Xall.res <- matrix(NA, nrow = nrow(wings.parents2), ncol = 96)

for (ln in levels(as.factor(wings.parents2$Genotype))) {
  ind <- which(wings.parents2$Genotype==ln) #In this case line is the cross or parent genotype
  Xln <- as.matrix(wings.parents2[ind,4:99]) # extracting out the right rows (line) 
  lcs <- wings.parents2$logCS[ind] # log centroid size
  reg <- lm(Xln ~ lcs) # regression model
  Xm <- colMeans(predict(reg)) #CHECK IF ITS PREDICTING THE RIGHT WAY # predicted mean configuration at mean centroid size to add back to the residuals. Predict is fitting it by default at the mean for lcs.
  Xall.res[ind,] <- reg$residuals + matrix(Xm, nr=length(ind), nc=length(Xm), byrow=TRUE) # recentered on the line mean
}

#I will go with this one. It didn't make all that much difference but I feel good about it. 

png("../Figures/HighLowParent_male_wingBlur.png")
WingBlur3(Xall.res, groups = TRUE, droplevels(wings.parents2$Alt), blur_transparency = 0.02)
dev.off()

Xall.resF <- matrix(NA, nrow = nrow(wings.parents3), ncol = 96)

for (ln in levels(as.factor(wings.parents3$Genotype))) {
  ind <- which(wings.parents3$Genotype==ln) #In this case line is the cross or parent genotype
  Xln <- as.matrix(wings.parents3[ind,4:99]) # extracting out the right rows (line) 
  lcs <- wings.parents3$logCS[ind] # log centroid size
  reg <- lm(Xln ~ lcs) # regression model
  Xm <- colMeans(predict(reg)) #CHECK IF ITS PREDICTING THE RIGHT WAY # predicted mean configuration at mean centroid size to add back to the residuals. Predict is fitting it by default at the mean for lcs.
  Xall.resF[ind,] <- reg$residuals + matrix(Xm, nr=length(ind), nc=length(Xm), byrow=TRUE) # recentered on the line mean
}


png("../Figures/HighLowParent_female_wingBlur.png")
WingBlur3(Xall.resF, groups = TRUE, droplevels(wings.parents3$Alt), blur_transparency = 0.02)
dev.off()


  
#adding a map as well 
#this will work but needs an API key and I don't have that sort of brain power right now. 

# library(ggmap)
# 
# map <- get_googlemap("Montpellier, France", zoom = 8, maptype = "terrain")
# 
# 
# ggmap(map) + 
#   theme_void() + 
#   ggtitle("terrain") + 
#   theme(
#     plot.title = element_text(colour = "orange"), 
#     panel.border = element_rect(colour = "grey", fill=NA, size=2)
#   )



parentmales.raw <- image_read("../Figures/HighLowParent_male_wingBlur.png")
parentmales.crop <- image_trim(parentmales.raw)
parentmales.crop.r <- image_flop(parentmales.crop)
parentmales.blur <- ggdraw() + draw_image(parentmales.crop.r)


# parentnoCor.raw <- image_read("../Figures/malesPureBlur_noCorrection.png")
# parentnoCor.crop <- image_trim(parentnoCor.raw)
# parentnoCor.crop.r <- image_flop(parentnoCor.crop)
# parentnoCor.blur <- ggdraw() + draw_image(parentnoCor.crop.r)

parentFemales.raw <- image_read("../Figures/HighLowParent_female_wingBlur.png")
parentFemales.crop <- image_trim(parentFemales.raw)
parentFemales.crop.r <- image_flop(parentFemales.crop)
parentFemales.blur <- ggdraw() + draw_image(parentFemales.crop.r)

wing_pannel <- plot_grid(parentmales.blur, parentFemales.blur, labels = c("Male", "Female"), ncol = 1)


png("../Figures/ParentDiffs_PureLines.png", height = 1500, width = 3000, units = "px", res = 300)
plot_grid(parent.size, wing_pannel, nrow = 1, rel_widths = c(1,0.4), greedy = FALSE, labs = c("A", ))

dev.off()


