#Remaking the female plots for the supplemental


###I think my renumbering code will still work for this one#####
library(data.table)
library(tidyverse)
library(ggnewscale)
library(cowplot)
source("KP_genomescan_source.R")

labs <- c("chrom", "start", "end", "snps", "fst", "number")

#################################################################
#this one is seperate because the females got left out of the sync
zi192ef96.shape <- fread("../Data/zi192ef96_shapeOnly_fst_genedalf.csv")

#Females then males. 
zi192ef96.shape2 <- zi192ef96.shape[,c(1:4, 5,10)]


zi192ef96.shape2.arranged <- chrNumbering(zi192ef96.shape2)
zi192ef96.shape.middle <- middleChr(zi192ef96.shape2.arranged)


zi192ef96.shapeF <- (zi192ef96.shape2.arranged%>% 
                       select(c(1:4, 5,7)) %>%
                       filter(snps >= 10) ) #this line deals with the NA in the data.

names(zi192ef96.shapeF)<- labs

zi192ef96.shapeM <- (zi192ef96.shape2.arranged%>% 
                       select(c(1:4, 6,7)) %>%
                       filter(snps >= 10) ) #this line deals with the NA in the data.

names(zi192ef96.shapeM)<- labs

#converts all negative numbers to 0
zi192ef96.shapeM$fst <- pmax(zi192ef96.shapeM$fst, 0)
zi192ef96.shapeF$fst <- pmax(zi192ef96.shapeF$fst, 0)


zi192ef96shape.plot <- ggplot(data = zi192ef96.shapeM, aes(x=number, y=fst, color=chrom)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  scale_x_discrete(limits=c(zi192ef96.shape.middle ),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  new_scale_colour() +
  geom_point(data = zi192ef96.shapeF, aes(x =number, y = fst, col = chrom), 
             alpha = 0.08, show.legend = F,) +
  scale_colour_manual(values=c('tomato', 'coral3', 'tomato', 'coral3', 'tomato', 'coral3')) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  theme(panel.background = element_blank())  + 
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 


#################################################################
zi192ef81 <- fread("../Data/zi192_ef81_fst_noMerge.csv")

#I am interested in the 1.4, 2.5 and 3.7 cols. they are indexed as 7, 13, 19. I also want the first 4 cols.
zi192ef81 <- zi192ef81[,c(1:4,7, 13, 19)]


head(zi192ef81)


zi192ef81.arranged <- chrNumbering(zi192ef81)
zi192ef81.middle <- middleChr(zi192ef81.arranged)

#I think I want to pull these out individually and then drop NaN cols 

labs <- c("chrom", "start", "end", "snps", "fst", "number")

zi192ef81.shapeF <- (zi192ef81.arranged %>% 
                       select(c(1:5,8)) %>%
                       filter(snps >= 10) ) #this line deals with the NA in the data.

names(zi192ef81.shapeF)<- labs


#missing values becaue of 0 limit. Does look cleaner than popoolation. 
zi192ef81.shapeFplot <- ggplot(data = zi192ef81.shapeF, aes(x=number, y=fst, color=chrom)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  scale_x_discrete(limits=c(zi192ef81.middle),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(-0.3, 0.6) +
  theme(panel.background = element_blank())  + 
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 




zi192ef81.shapeM <- (zi192ef81.arranged %>% 
                       select(c(1:4, 6 ,8)) %>%
                       filter(snps >= 10) ) #this line deals with the NA in the data.

names(zi192ef81.shapeM)<- labs


zi192ef81.shapeMplot <- ggplot(data = zi192ef81.shapeM, aes(x=number, y=fst, color=chrom)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  scale_x_discrete(limits=c(zi192ef81.middle),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(-0.3, 0.6) +
  theme(panel.background = element_blank())  + 
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 


zi192ef81.sizeM <- (zi192ef81.arranged %>% 
                      select(c(1:4, 7 ,8)) %>%
                      filter(snps >= 10) ) #this line deals with the NA in the data.

names(zi192ef81.sizeM)<- labs


zi192ef81.sizeMplot <- ggplot(data = zi192ef81.sizeM, aes(x=number, y=fst, color=chrom)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  scale_x_discrete(limits=c(zi192ef81.middle),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(-0.3, 0.6) +
  theme(panel.background = element_blank())  + 
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 


#converts all negative numbers to 0
zi192ef81.shapeM$fst <- pmax(zi192ef81.shapeM$fst, 0)
zi192ef81.shapeF$fst <- pmax(zi192ef81.shapeF$fst, 0)

zi192ef81shape.plot <- ggplot(data = zi192ef81.shapeM, aes(x=number, y=fst, color=chrom)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  scale_x_discrete(limits=c(zi192ef96.shape.middle ),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  new_scale_colour() +
  geom_point(data = zi192ef81.shapeF, aes(x =number, y = fst, col = chrom), 
             alpha = 0.08, show.legend = F,) +
  scale_colour_manual(values=c('tomato', 'coral3', 'tomato', 'coral3', 'tomato', 'coral3')) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  theme(panel.background = element_blank())  + 
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 


#################################################################

zi418ef43 <- fread("../Data/zi418ef43_fst_noMerge.csv")

#only doing shape because I am missing size data (and its all messed up)
#these are 2:6 (F) col 16 and 3:7 (M) col 23. I also want the first 4 cols.
zi418ef43 <- zi418ef43[,c(1:4, 16, 23)]


head(zi418ef43)


zi418ef43.arranged <- chrNumbering(zi418ef43)
zi418ef43.middle <- middleChr(zi418ef43.arranged)

#I think I want to pull these out individually and then drop NaN cols 

labs <- c("chrom", "start", "end", "snps", "fst", "number")

zi418ef43.shapeF <- (zi418ef43.arranged %>% 
                       select(c(1:5,7)) %>%
                       filter(snps >= 10) ) #this line deals with the NA in the data.

names(zi418ef43.shapeF)<- labs


#missing values becaue of 0 limit. Does look cleaner than popoolation. 
zi418ef43.shapeFplot <- ggplot(data = zi418ef43.shapeF, aes(x=number, y=fst, color=chrom)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  scale_x_discrete(limits=c(zi192ef81.middle),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(-0.3, 0.6) +
  theme(panel.background = element_blank())  + 
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 




zi418ef43.shapeM <- (zi418ef43.arranged %>% 
                       select(c(1:4, 6,7)) %>%
                       filter(snps >= 10) ) #this line deals with the NA in the data.

names(zi418ef43.shapeM)<- labs


zi418ef43.shapeMplot <- ggplot(data = zi418ef43.shapeM, aes(x=number, y=fst, color=chrom)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  scale_x_discrete(limits=c(zi418ef43.middle),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(-0.3, 0.6) +
  theme(panel.background = element_blank())  + 
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 


#converts all negative numbers to 0
zi418ef43.shapeM$fst <- pmax(zi418ef43.shapeM$fst, 0)
zi418ef43.shapeF$fst <- pmax(zi418ef43.shapeF$fst, 0)

zi418ef43shape.plot <- ggplot(data = zi418ef43.shapeM, aes(x=number, y=fst, color=chrom)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  scale_x_discrete(limits=c(zi418ef43.middle),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  new_scale_colour() +
  geom_point(data = zi418ef43.shapeF, aes(x =number, y = fst, col = chrom), 
             alpha = 0.08, show.legend = F,) +
  scale_colour_manual(values=c('tomato', 'coral3', 'tomato', 'coral3', 'tomato', 'coral3')) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  theme(panel.background = element_blank())  + 
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 


png("../Figures/shapeFst_bothSex_5000windows_geneleaf.png", width =4000, height = 4000, units = "px",res = 200)
plot_grid(zi192ef96shape.plot, zi192ef81shape.plot, zi418ef43shape.plot, ncol = 1, 
          labels = c("Zi192 x Ef96","Zi192 x Ef81", "Zi418 x Ef43"), label_size = 25)

dev.off()
