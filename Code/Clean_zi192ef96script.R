#Analysis for selection of lines (adapted from Rmd file). This is missing some of the checks done in that file but is meant to be a more clean file that can be run quickly from top to bottom to work with this data. 
#For new data/populations this should not be all the code run, more checks are needed!! 

load("../Data/zi192ef96_allwings.rda")
source("~/Dropbox/KatiePelletier/KP_geomorphwingfunctions.R")

########
library(tidyverse)
library(MASS)
library(geomorph)
library(magick)
library(gridGraphics)
#######

######### Loading functions #########
source('~/Dropbox/DworkinLabSharedMaterial/scripts/WRP_FUNCTIONS.R', chdir = TRUE)

source('~/Dropbox/DworkinLabSharedMaterial/scripts/WINGPLOTSOURCE.R', chdir = TRUE)

####################################

#A look at the wings 
PCs <- prcomp(wings[, c(4:99, 104)])
wings_combined <- data.frame(wings, PCs$x[,1:57]) 

#this isn't particularly useful, just looking at some parental lines. 

wingsHAF <- filter(wings_combined, sex == "f" & Alt == "HA")
wingsLAF <- filter(wings_combined, sex == "f" & Alt == "LA")
wingsHAM <- filter(wings_combined, sex == "m" & Alt == "HA")
wingsLAM <- filter(wings_combined, sex == "m" & Alt == "LA")

HAmeanF <- colMeans(wingsHAF[,106:160])
LAmeanF <- colMeans(wingsLAF[,106:160])
HAmeanM <- colMeans(wingsHAM[,106:160])
LAmeanM <- colMeans(wingsLAM[,106:160])

#This is using the PCs from a PCA including all sexes together. 
Fdiff <- HAmeanF - LAmeanF
Mdiff <- HAmeanM - LAmeanM
#0.91. The diffrence between populations for males and females is really similar. (with allometry)
#0.85 when logCS is included in PCA and then PC1 is dropped. 
#There is something going on here with allometry - probably not the same between the sexes. 
cor(Fdiff, Mdiff)

ggplot(wings_combined, aes(x = PC2, y = PC3, col = Alt, shape =  sex)) +
  geom_point(alpha = 0.6, size = 2 )


ggplot(wings_combined, aes(x = PC3, y = PC4, col = Alt, shape =  sex)) +
  geom_point(alpha = 0.6, size = 2 )

##########################
#We chose to look at/select the sexes seperatley. 

wings_female <- subset(wings, sex == "f")
wings_male <- subset(wings, sex == "m")


pureHAF <- (wings %>%
              subset(sex=="f" & Alt == "HA") 
)
pureHAM <- (wings %>%
              subset(sex=="m" & Alt == "HA") 
)

pureLAF <- (wings %>%
              subset(sex=="f" & Alt == "LA") 
)
pureLAM <- (wings %>%
              subset(sex=="m" & Alt == "LA") 
)

HAmean <- colMeans(pureHAF[,4:99])
LAmean <- colMeans(pureLAF[,4:99])
Fdiff <- HAmean - LAmean

HAmeanM <- colMeans(pureHAM[,4:99])
LAmeanM <- colMeans(pureLAM[,4:99])
Mdiff <- HAmeanM - LAmeanM

#Still very high, this is using only the landmarks. 
cor(Fdiff, Mdiff)

########Female selection first 
PC_female <- prcomp(wings_female[, c(4:99, 104)]) 
#I included logCS in this, but PC1 doesn't explain as much varriance as I would expect (55%)
summary(PC_female)

wings_female <- data.frame(wings_female, PC_female$x[,1:56])

#PC1 is still the allometric componet 

#I want to print this as a supplemental figure to show that PC1 is capturing the allometric componet of shape. 

ggplot(wings_female, aes(logCS, PC1, col = Alt)) + 
  geom_point(alpha = 0.5) + 
  theme_classic()
#This looks ok. 
ggplot(wings_female, aes(PC3, PC2, col = Alt)) + geom_point(alpha = 0.5) 

#doing my selection in PC space. 
#This is PC 2-56
HA_mean_F <- colMeans(wings_female[wings_female$Alt == "HA", 106:160])
LA_mean_F <- colMeans(wings_female[wings_female$Alt == "LA", 106:160])
diff_mean_F <- HA_mean_F - LA_mean_F

wings_female$diff <-  as.matrix(wings_female[,106:160]) %*% diff_mean_F

female_plot <- ggplot(wings_female, aes(x = CS, y = diff, col = Alt)) +
  geom_point(alpha = 0.6, size = 2 )

##Now to select the top and bottom 50 for shape and size 

wings_cross_F <- (wings_female %>% filter(Genotype == "zi192xef96"))

#selecting 52 incase there is a problem when pulling one of the flies (bad disection ect. these can be ignored later), then I have a couple backups
female.top50shape <- (wings_cross_F %>%
                        top_n(52, diff) %>%
                        mutate(class.shape = "highshape"))

female.bottom50shape <- (wings_cross_F %>%
                           top_n(-52, diff) %>%
                           mutate(class.shape = "lowshape"))

female.shape <- rbind(female.bottom50shape, female.top50shape)
female.sel <- rbind(female.bottom50shape, female.top50shape)
female.shape$shaperank <- rank(female.shape$diff)

female.top50size <- (wings_cross_F %>%
                       top_n(52, CS) %>%
                       mutate(class.size = "lowsize"))


female.bottom50size <- (wings_cross_F %>%
                          top_n(-52, CS) %>%
                          mutate(class.size = "highsize"))



female.size <- rbind(female.bottom50size, female.top50size)
female.size$sizerank <- rank(female.shape$CS)


female.sel <- (wings_cross_F %>%
                 left_join(female.shape) %>%
                 left_join(female.size) %>% 
                 mutate(class = paste(class.shape, class.size, sep = ".")))

female.selection <- (female.sel %>%
                       filter(class != "NA.NA") %>%
                       dplyr::select(Fly_ID, sex, CS, diff, class))

with(female.sel, table(class))

#adding high/low labels 
library(grid)

lowshape <- grobTree(textGrob("Low Altitude Shape", x=0.4,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=10)))
highshape <- grobTree(textGrob("High Altitude Shape", x=0.4,  y=0.05, hjust=0, gp=gpar(col="black", fontsize=10)))
lowsize <- grobTree(textGrob("Low Altitude Size", x=0.05,  y=0.5, hjust=0, gp=gpar(col="black", fontsize=10)))
highsize <- grobTree(textGrob("High Altitude Size", x=0.95,  y=0.5, hjust=0, gp=gpar(col="black", fontsize=10)))

f.sel <- ggplot(female.sel, aes(x = CS, y = diff, col = class)) +
  geom_point(alpha = 0.5) + 
  ylab("High-Low Shape Diffrence") + 
  xlab("Centroid Size") + 
  scale_colour_manual(values=c('red', 'blue', 'green', 'deeppink', 'darkviolet', 'darkorange', 'grey46')) +
  theme(panel.background = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        legend.key = element_blank(),
        text = element_text(size=15),
        axis.text.x= element_text(size=12), 
        axis.text.y= element_text(size=12))

pdf("../Outputs/zi192ef96femalesforseq.pdf")
f.sel
dev.off()

#Now I want to calculate the average shape of the left and right and plot these. 

HAselShape <- colMeans(female.top50shape[,4:99])
LAselShape <- colMeans(female.bottom50shape[,4:99])

shapeeffect <- HAselShape - LAselShape

WingPlot(HAselShape, wingcol="blue")
WingPlot(LAselShape, wingcol = "black")


# pdf("../Outputs/SelectiondiffshapeF.pdf")
# WingEffect(HAselShape, shapeeffect, shapeeffect, scale.factor=1, wingframe = F, scale.display = F)
# dev.off()

##and in geomorph for less fucky wings. 

#I am going to make a df with the the left and right pool mean shapes and the mean shape because I am going to plot the selection effect while I am here and doing things. 

m <- colMeans(wings_cross_F[,4:99])

plotting.f <- as.data.frame(rbind(m, HAselShape, LAselShape))
plotting.f$crap <- c("mshape", "HA", "LA")

cord <- as.matrix(plotting.f[,1:96])
shape <- arrayspecs(cord, 48, 2)
f.plot.wings <- geomorph.data.frame(shape = shape,
                            crap = plotting.f$crap)

#Just checking because this needs to be a mshape object for plotting. 
mean.shape <- mshape(f.plot.wings$shape)

plot(mean.shape, links = wing.links)

#Need to back transform the selection vector? not sure how to handle that with the dropped PC? Just use the other one for plotting?

png("../Figures/zi192ef96_RvL_meanshape_2x.png")
plotRefToTarget(f.plot.wings$shape[,,2], 
                f.plot.wings$shape[,,3], 
                links = wing.links, method = "points", mag = 2, 
                gridPars=wing.spec )
dev.off()


######Males

PC_male <- prcomp(wings_male[, c(4:99, 104)])
#this is also true for males, only 53% explained by PC1
summary(PC_male)
wings_male <- data.frame(wings_male, PC_male$x[,1:56])

HA_mean_M <- colMeans(wings_male[wings_male$Alt == "HA", 106:160])
LA_mean_M <- colMeans(wings_male[wings_male$Alt == "LA", 106:160])
diff_mean_M <- HA_mean_M - LA_mean_M

#This is diffrent, but they are two diffrent PCAs so of course. Right?
cor(diff_mean_M, diff_mean_F)

wings_male$diff <-  as.matrix(wings_male[,106:160]) %*% diff_mean_M

male_plot <- ggplot(wings_male, aes(x = CS, y = diff, col = Alt)) +
  geom_point(alpha = 0.6, size = 2 )


wings_cross_m <- (wings_male %>% filter(Genotype == "zi192xef96"))

male.top50shape <- (wings_cross_m %>%
                      top_n(52, diff) %>%
                      mutate(class.shape = "lowshape"))

male.bottom50shape <- (wings_cross_m %>%
                         top_n(-52, diff) %>%
                         mutate(class.shape = "highshape"))

male.shape <- rbind(male.bottom50shape, male.top50shape)
male.shape$shaperank <- rank(male.shape$diff)


male.top50size <- (wings_cross_m %>%
                     top_n(52, CS) %>%
                     mutate(class.size = "highsize"))

male.bottom50size <- (wings_cross_m %>%
                        top_n(-52, CS) %>%
                        mutate(class.size = "lowsize"))

male.size <- rbind(male.bottom50size, male.top50size)
male.shape$sizerank <- rank(male.shape$CS)


male.sel <- (wings_cross_m %>%
               left_join(male.shape) %>%
               left_join(male.size) %>% 
               mutate(class = paste(class.shape, class.size, sep = ".")))

male.selection <- (male.sel %>%
                     filter(class != "NA.NA") %>%
                     dplyr::select(Fly_ID, sex, CS,diff, class))

m.sel <- ggplot(male.sel, aes(x = CS, y = diff, col = class)) + geom_point()

selection <- rbind(female.selection, male.selection)

write_csv(selection, "../zi192ef96_forSequencing.csv")



###############Some paper plots###################################

library(cowplot)
#female and male cs vs PC1 to show allometry. 


f.allo <- ggplot(wings_female, aes(logCS, PC1, col = Alt)) + 
  geom_point(alpha = 0.3) + 
  theme_classic() + 
  theme(legend.title = element_blank()) + 
  xlab("Centroid Size") +
  scale_color_manual(values = c("black", "red", "darkgray"), 
                     name = "Alt", 
                     labels = c("Low Altitude", "F20 progeny", "High Altitude"))  

m.allo <- ggplot(wings_male, aes(logCS, PC1, col = Alt)) + 
  geom_point(alpha = 0.3) + 
  theme_classic() + 
  #theme(legend.position = "none") + 
  xlab("Centroid Size") +
  scale_color_manual(values = c("black", "red", "darkgray"), 
                     name = "Alt", 
                     labels = c("Low Altitude", "F20 progeny", "High Altitude"))  

allo.legend <- get_legend(f.allo)

allo.plot <- plot_grid(m.allo + 
                         theme(legend.position = "none"), 
                       f.allo + 
                         theme(legend.position = "none"), 
              labels = c("Males", "Females"), 
              scale = c(1, 1))


png("../Figures/zi192ef96_bothsex_PC1CS_alloPlot.png", 
    width =2600, height = 1000, 
    units = "px", res = 300)
plot_grid(allo.plot, allo.legend, rel_widths = c(3, .4))
dev.off()


#####Size and shape change plots. 

f.sel.plot <-ggplot(female.sel, aes(x = CS, y = diff, col = class)) +
  geom_point(alpha = 0.5, size = 0.75) + 
  ylab("Shape Score") + 
  xlab("Centroid Size") + 
  scale_colour_manual(values=c('red', 'red', 'red', 'red', 'red', 'red', 'darkgrey')) +
    geom_vline(aes(xintercept = min(female.top50size$CS))) + 
    geom_vline(aes(xintercept = max(female.bottom50size$CS))) +
    geom_hline(aes(yintercept = max(female.bottom50shape$diff))) +
    geom_hline(aes(yintercept = min(female.top50shape$diff))) +
  theme(panel.background = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        legend.key = element_blank(),
        text = element_text(size=15),
        axis.text.x= element_text(size=12), 
        axis.text.y= element_text(size=12))


f.RvL.raw <- image_read("../Figures/zi192ef96_RvL_meanshape_2x.png")
f.RvL.crop <- image_trim(f.RvL.raw)
f.RvL.flop <- image_flop(f.RvL.crop)
f.RvL <- ggdraw() + draw_image(f.RvL.flop)

f.UvD.raw <- image_read("../Figures/Fcomparesize.png")
f.UvD.crop <- image_trim(f.UvD.raw)
f.UvD <- ggdraw() + draw_image(f.UvD.crop)

shape.size.col <- plot_grid(f.RvL, f.UvD, ncol = 1,
                            labels = c("B", "C"), 
                            scale = c(0.6, 1))


png("../Figures/zi192ef96_f_selectionFig.png", 
    width =3000, height = 1800, 
    units = "px", res = 300)

plot_grid(f.sel.plot, shape.size.col, 
          labels = c("A", ""), 
          scale = c(1, 0.9))

dev.off()

m.sel.plot <-ggplot(male.sel, aes(x = CS, y = diff, col = class)) +
  geom_point(alpha = 0.5, size = 0.75) + 
  ylab("Shape Score") + 
  xlab("Centroid Size") + 
  scale_colour_manual(values=c('red', 'red', 'red', 'red', 'red', 'red', 'red', 'darkgrey')) +
  geom_vline(aes(xintercept = min(male.top50size$CS))) + 
  geom_vline(aes(xintercept = max(male.bottom50size$CS))) +
  geom_hline(aes(yintercept = max(male.bottom50shape$diff))) +
  geom_hline(aes(yintercept = min(male.top50shape$diff))) +
  theme(panel.background = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        legend.key = element_blank(),
        text = element_text(size=15),
        axis.text.x= element_text(size=12), 
        axis.text.y= element_text(size=12))

png("../Figures/zi192ef96_m_selectionFig.png", 
    width =1500, height = 1500, 
    units = "px", res = 300)

m.sel.plot

dev.off()


png("../Figures/zi192ef96_bothSex_selection.png", 
    width =1700, height = 750, 
    units = "px", res = 300)

plot_grid(f.sel.plot, m.sel.plot, labels = c("Female", "Male"))

dev.off()

