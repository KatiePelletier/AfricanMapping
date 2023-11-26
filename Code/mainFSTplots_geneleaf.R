#Remaking the plots for the paper (this uses only male data for boring and tecnical reasons. Both sexes are in the supp.)


###I think my renumbering code will still work for this one#####
library(data.table)
library(tidyverse)
library(ggnewscale)
library(cowplot)
source("KP_genomescan_source.R")

labs <- c("chrom", "start", "end", "snps", "fst", "number")

#################################################################

zi192ef81 <- fread("../Data/zi192_ef81_fst_noMerge.csv")


#I am interested in the 2.5 and 3.7 cols. they are indexed as 13, 19. I also want the first 4 cols.
zi192ef81 <- zi192ef81[,c(1:4, 13, 19)]


zi192ef81.arranged <- chrNumbering(zi192ef81)
zi192ef81.middle <- middleChr(zi192ef81.arranged)


zi192ef81.shapeM <- (zi192ef81.arranged %>% 
                       select(c(1:4, 5 ,7)) %>%
                       filter(snps >= 10) ) #this line deals with the NA in the data.

names(zi192ef81.shapeM)<- labs

#converts all negative numbers to 0
zi192ef81.shapeM$fst <- pmax(zi192ef81.shapeM$fst, 0)


zi192ef81.shapeMplot <- ggplot(data = zi192ef81.shapeM, aes(x=number, y=fst, color=chrom)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  scale_x_discrete(limits=c(zi192ef81.middle),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  theme(panel.background = element_blank())  + 
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 


zi192ef81.sizeM <- (zi192ef81.arranged %>% 
                      select(c(1:4, 6 ,7)) %>%
                      filter(snps >= 10) ) #this line deals with the NA in the data.

names(zi192ef81.sizeM)<- labs

#converts all negative numbers to 0
zi192ef81.sizeM$fst <- pmax(zi192ef81.sizeM$fst, 0)

zi192ef81.sizeMplot <- ggplot(data = zi192ef81.sizeM, aes(x=number, y=fst, color=chrom)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  scale_x_discrete(limits=c(zi192ef81.middle),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  theme(panel.background = element_blank())  + 
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 


#################################################################
zi192ef96 <- fread("../Data/zi192ef96_fst_noMerge.csv")


#I am interested in the 2.5 and 4.7 cols. they are indexed as 13 and 22. I also want the first 4 cols.
zi192ef96 <- zi192ef96[,c(1:4, 13, 22)]

head(zi192ef96)


zi192ef96.arranged <- chrNumbering(zi192ef96)
zi192ef96.middle <- middleChr(zi192ef96.arranged)


zi192ef96.shapeM <- (zi192ef96.arranged %>% 
                       select(c(1:4, 5 ,7)) %>%
                       filter(snps >= 10) ) #this line deals with the NA in the data.

names(zi192ef96.shapeM)<- labs

#converts all negative numbers to 0
zi192ef96.shapeM$fst <- pmax(zi192ef96.shapeM$fst, 0)


zi192ef96.shapeMplot <- ggplot(data = zi192ef96.shapeM, aes(x=number, y=fst, color=chrom)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  scale_x_discrete(limits=c(zi192ef96.middle),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  theme(panel.background = element_blank())  + 
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 


zi192ef96.sizeM <- (zi192ef96.arranged %>% 
                      select(c(1:4, 6 ,7)) %>%
                      filter(snps >= 10) ) #this line deals with the NA in the data.

names(zi192ef96.sizeM)<- labs

#converts all negative numbers to 0
zi192ef96.sizeM$fst <- pmax(zi192ef96.sizeM$fst, 0)


zi192ef96.sizeMplot <- ggplot(data = zi192ef96.sizeM, aes(x=number, y=fst, color=chrom)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  scale_x_discrete(limits=c(zi192ef81.middle),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0.001, 0.6) +
  theme(panel.background = element_blank())  + 
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 

#################################################################
zi418ef43 <- fread("../Data/zi418ef43_fst_noMerge.csv")


#I am interested in the 3.7 and 5.9 cols. they are indexed as 23 and 34. I also want the first 4 cols.
zi418ef43 <- zi418ef43[,c(1:4, 23, 34)]


head(zi418ef43)

zi418ef43.arranged <- chrNumbering(zi418ef43)
zi418ef43.middle <- middleChr(zi418ef43.arranged)


zi418ef43.shapeM <- (zi418ef43.arranged %>% 
                       select(c(1:4, 5 ,7)) %>%
                       filter(snps >= 10) ) #this line deals with the NA in the data.

names(zi418ef43.shapeM)<- labs

#converts all negative numbers to 0
zi418ef43.shapeM$fst <- pmax(zi418ef43.shapeM$fst, 0)


zi418ef43.shapeMplot <- ggplot(data = zi418ef43.shapeM, aes(x=number, y=fst, color=chrom)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  scale_x_discrete(limits=c(zi418ef43.middle),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  theme(panel.background = element_blank())  + 
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 


zi418ef43.sizeM <- (zi418ef43.arranged %>% 
                      select(c(1:4, 6 ,7)) %>%
                      filter(snps >= 10) ) #this line deals with the NA in the data.

names(zi418ef43.sizeM)<- labs

#converts all negative numbers to 0
zi418ef43.sizeM$fst <- pmax(zi418ef43.sizeM$fst, 0)


zi418ef43.sizeMplot <- ggplot(data = zi418ef43.sizeM, aes(x=number, y=fst, color=chrom)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  scale_x_discrete(limits=c(zi192ef81.middle),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  theme(panel.background = element_blank())  + 
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 

####################################################################################
#Now to arrange these all into the final figures for the paper. 

png("../Figures/shapeFst_malesOnly_5000windows_geneleaf.png", width =4000, height = 4000, units = "px",res = 200)
plot_grid(zi192ef96.shapeMplot, zi192ef81.shapeMplot, zi418ef43.shapeMplot, ncol = 1, 
          labels = c("Zi192 x Ef96","Zi192 x Ef81", "Zi418 x Ef43"), label_size = 25)

dev.off()

png("../Figures/sizeFst_malesOnly_5000windows_geneleaf.png", width =4000, height = 4000, units = "px",res = 200)
plot_grid(zi192ef96.sizeMplot, zi192ef81.sizeMplot, zi418ef43.sizeMplot, ncol = 1, 
          labels = c("Zi192 x Ef96","Zi192 x Ef81", "Zi418 x Ef43"), label_size = 25)

dev.off()



####################################################################################
#plotting out the 3R region for the supplemental 

#before I just pulled out all of 3R and plotted that. 

zi418ef43.shapeM.3R <- filter(zi418ef43.shapeM, chrom == "3R")
zi418ef43.shapeM.3R$cross <- "zi418ef43"
zi192ef96.shapeM.3R <- filter(zi192ef96.shapeM, chrom == "3R")
zi192ef96.shapeM.3R$cross <- "zi192ef96"
zi192ef81.shapeM.3R <- filter(zi192ef81.shapeM, chrom == "3R")
zi192ef81.shapeM.3R$cross <- "zi192ef81"

all.3R <- rbind(zi192ef96.shapeM.3R, zi192ef81.shapeM.3R, zi418ef43.shapeM.3R)

#start marks the start position of that window on 3R
png("../Figures/3Rfst_overplotting_genedalf.png", width =2000, height = 1000, units = "px",res = 200)
ggplot(all.3R, aes(x = start, y = fst, col = cross)) + 
  geom_point(size=1, alpha = 0.4) + 
  scale_colour_manual(values=c('grey25', 'grey46', 'black') ) +
  theme(panel.background = element_blank())  + 
  xlab("") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 
dev.off()

####################################################################################
#I want to make a list of genes within the top 5% of sites to compare with the CMH and give some maybe candidates?

#the first thing I need to do is get the windows with the top 5% of Fst. I did this in the first paper already. 
#this will also let me look for oultier windows. 

quantile(zi418ef43.shapeM$fst, 0.99 ,na.rm = TRUE)

#262
nrow(filter(zi418ef43.shapeM, fst >= quantile(zi418ef43.shapeM$fst, 0.99 ,na.rm = TRUE)))

#254
nrow(filter(zi192ef96.shapeM, fst >= quantile(zi192ef96.shapeM$fst, 0.99 ,na.rm = TRUE)))

#253
nrow(filter(zi192ef81.shapeM, fst >= quantile(zi192ef81.shapeM$fst, 0.99 ,na.rm = TRUE)))


zi418ef43.shapeM.outlier <- filter(zi418ef43.shapeM, fst >= quantile(zi418ef43.shapeM$fst, 
                                                              0.99 ,na.rm = TRUE))
zi418ef43.shapeM.outlier$windowID <- paste(zi418ef43.shapeM.outlier$chrom, 
                                           zi418ef43.shapeM.outlier$start, sep = ".")

zi192ef96.shapeM.outlier <- filter(zi192ef96.shapeM, fst >= quantile(zi192ef96.shapeM$fst, 
                                                                     0.99 ,na.rm = TRUE))
zi192ef96.shapeM.outlier$windowID <- paste(zi192ef96.shapeM.outlier$chrom, 
                                           zi192ef96.shapeM.outlier$start, sep = ".")


zi192ef81.shapeM.outlier <- filter(zi192ef81.shapeM, fst >= quantile(zi192ef81.shapeM$fst, 
                                                                     0.99 ,na.rm = TRUE))
zi192ef81.shapeM.outlier$windowID <- paste(zi192ef81.shapeM.outlier$chrom, 
                                           zi192ef81.shapeM.outlier$start, sep = ".")

#number shared windows. 

#0 
length(which(zi418ef43.shapeM.outlier$windowID %in%  zi192ef96.shapeM.outlier$windowID))
#same
length(which(zi192ef96.shapeM.outlier$windowID %in%  zi418ef43.shapeM.outlier$windowID))

#2
length(which(zi418ef43.shapeM.outlier$windowID %in%  zi192ef81.shapeM.outlier$windowID))
#2
length(which(zi192ef81.shapeM.outlier$windowID %in%  zi418ef43.shapeM.outlier$windowID))


#141
length(which(zi192ef96.shapeM.outlier$windowID %in%  zi192ef81.shapeM.outlier$windowID))

zi192ef96.shapeM.outlier$windowID[which(zi192ef96.shapeM.outlier$windowID %in%  zi192ef81.shapeM.outlier$windowID)]
#141
length(which(zi192ef81.shapeM.outlier$windowID %in%  zi192ef96.shapeM.outlier$windowID))

#clearly none shared between all three. 

#I don't actually think this annotation of genes is useful because the region is so fucking big and has a lot of other unrelated stuff. 


#overlap for size. 

#same as above because it is a total of the whole genome. 
#262
nrow(filter(zi418ef43.sizeM, fst >= quantile(zi418ef43.sizeM$fst, 0.99 ,na.rm = TRUE)))

#254
nrow(filter(zi192ef96.sizeM, fst >= quantile(zi192ef96.sizeM$fst, 0.99 ,na.rm = TRUE)))

#253
nrow(filter(zi192ef81.sizeM, fst >= quantile(zi192ef81.sizeM$fst, 0.99 ,na.rm = TRUE)))


zi418ef43.sizeM.outlier <- filter(zi418ef43.sizeM, fst >= quantile(zi418ef43.sizeM$fst, 
                                                                     0.99 ,na.rm = TRUE))
zi418ef43.sizeM.outlier$windowID <- paste(zi418ef43.sizeM.outlier$chrom, 
                                           zi418ef43.sizeM.outlier$start, sep = ".")

zi192ef96.sizeM.outlier <- filter(zi192ef96.shapeM, fst >= quantile(zi192ef96.sizeM$fst, 
                                                                     0.99 ,na.rm = TRUE))
zi192ef96.sizeM.outlier$windowID <- paste(zi192ef96.sizeM.outlier$chrom, 
                                           zi192ef96.sizeM.outlier$start, sep = ".")


zi192ef81.sizeM.outlier <- filter(zi192ef81.sizeM, fst >= quantile(zi192ef81.sizeM$fst, 
                                                                     0.99 ,na.rm = TRUE))
zi192ef81.sizeM.outlier$windowID <- paste(zi192ef81.sizeM.outlier$chrom, 
                                           zi192ef81.sizeM.outlier$start, sep = ".")

#number shared windows. 

#4
length(which(zi418ef43.sizeM.outlier$windowID %in%  zi192ef96.sizeM.outlier$windowID))


#3
length(which(zi418ef43.sizeM.outlier$windowID %in%  zi192ef81.sizeM.outlier$windowID))



#2
length(which(zi192ef96.sizeM.outlier$windowID %in%  zi192ef81.sizeM.outlier$windowID))

out1 <- zi192ef96.sizeM.outlier$windowID[which(zi192ef96.sizeM.outlier$windowID %in%  zi192ef81.sizeM.outlier$windowID)]
length(which(out1 %in% zi418ef43.sizeM.outlier$windowID))
