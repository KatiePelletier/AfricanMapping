#because of the weirdness with the third cross I want to plot this out with the crosses on top of one anohter. 

library(data.table)
library(tidyverse)
library(ggnewscale)
library(cowplot)
source("KP_genomescan_source.R")

###################################################################
#dataReadIn
###################################################################

zi192ef96.all <- fread("../Data/zi192ef96_nomerge_5000windows.fst")
#3:6 and 4:7 are the ones I care about. 
#col 19 and 23. also want the first 5

temp <- zi192ef96.all[,c(1:5, 19, 23)]

shitshit <- zi192ef96.all[,c(1:5, 10, 14)]

#one of these is the 

zi192ef96.size <-populaionFst_cleanup(temp, x = c("female", "male"))

zi192ef96.size <- chrNumbering(zi192ef96.size)
chrlabel.19296.size <- middleChr(zi192ef96.size)

zi192ef96.size.long <- pivot_longer(zi192ef96.size, cols = c("female", "male"),
                                     names_to = "sex", values_to = "meanFst")



#Something got overwritten here. I need to go back to my notes. 
zi192ef96.all2 <- fread("../Data/zi192ef96_shapeOnly_5000windows.fst")
#1:3 and 2:4 are the comps I want. col 7 and 10
temp.shape <- zi192ef96.all[,c(1:5, 7, 10)]

zi192ef96.shape <-populaionFst_cleanup(temp.shape, x = c("female", "male"))

zi192ef96.shape <- chrNumbering(zi192ef96.shape)
chrlabel.19296.shape <- middleChr(zi192ef96.shape)

zi192ef96.shape.long <- pivot_longer(zi192ef96.shape, cols = c("female", "male"),
                                    names_to = "sex", values_to = "meanFst")


cor(zi192ef96.size$male, zi192ef96.shape$male)

zi192ef81.all <- fread("../Data/zi192ef81_nomerge_5000windows.fst")

temp19281shape <- zi192ef81.all[,c(1:5, 8, 14)]

zi192ef81.shape <-populaionFst_cleanup(temp19281shape, x = c("female", "male"))

zi192ef81.shape <- chrNumbering(zi192ef81.shape)
chrlabel.19281.shape <- middleChr(zi192ef81.shape)

zi192ef81.shape.long <- pivot_longer(zi192ef81.shape, cols = c("female", "male"),
                                     names_to = "sex", values_to = "meanFst")

#only males. 
temp19281size <- zi192ef81.all[,c(1:5, 11, 20)]

brokeMycode <-populaionFst_cleanup(temp19281size, x = c("crap", "male"))

zi192ef81.size <- dplyr::select(brokeMycode, - "crap")

zi192ef81.size <- chrNumbering(zi192ef81.size)
chrlabel.19281.size <- middleChr(zi192ef81.size)

zi418ef43.all <- fread("../Data/zi418ef43_nomerge_5000windows.fst")

temp41843shape <- zi418ef43.all[,c(1:5, 17, 24)]

zi418ef43.shape <-populaionFst_cleanup(temp41843shape, x = c("female", "male"))

zi418ef43.shape <- chrNumbering(zi418ef43.shape)
chrlabel.41843.shape <- middleChr(zi418ef43.shape)

zi418ef43.shape.long <- pivot_longer(zi418ef43.shape, cols = c("female", "male"),
                                    names_to = "sex", values_to = "meanFst")

temp41843size <- zi418ef43.all[,c(1:5, 30, 35)]

zi418ef43.size <-populaionFst_cleanup(temp41843size, x = c("female", "male"))

zi418ef43.size <- chrNumbering(zi418ef43.size)
chrlabel.41843.size <- middleChr(zi418ef43.size)

zi418ef43.size.long <- pivot_longer(zi418ef43.size, cols = c("female", "male"),
                                    names_to = "sex", values_to = "meanFst")

###################################################################
#size plotting
###################################################################

#I want to plot the sexes on top of one another. 

zi192ef96.sizePlot <- ggplot(data = zi192ef96.size, aes(x=number, y=male, color=chr)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  scale_x_discrete(limits=c(chrlabel.19296.size),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  new_scale_colour() +
  geom_point(aes(y = female, col = chr), alpha = 0.08, show.legend = F,) +
  scale_colour_manual(values=c('tomato', 'coral3', 'tomato', 'coral3', 'tomato', 'coral3')) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  theme(panel.background = element_blank())  + 
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 

zi192ef81.sizePlot <- ggplot(data = zi192ef81.size, aes(x=number, y=male, color=chr)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  scale_x_discrete(limits=c(chrlabel.19296.size),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  #new_scale_colour() +
  #geom_point(aes(y = female, col = chr), alpha = 0.08, show.legend = F,) +
  #scale_colour_manual(values=c('tomato', 'coral3', 'tomato', 'coral3', 'tomato', 'coral3')) +
  theme(panel.background = element_blank())  + 
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 

#Why am I like this? What is going on? 
zi418ef43.sizePlot <- ggplot(data = zi418ef43.size, aes(x=number, y=male, color=chr)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) +
  scale_x_discrete(limits=c(chrlabel.19296.size),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  new_scale_colour() +
  geom_point(aes(y = female, col = chr), alpha = 0.08, show.legend = F,) +
  scale_colour_manual(values=c('tomato', 'coral3', 'tomato', 'coral3', 'tomato', 'coral3')) +
  theme(panel.background = element_blank())  + 
    xlab("Chromosome") +
    ylab("meanFst") +
    ylim(0, 0.6) +
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 


png("../Figures/sizeFst_bySex_5000wingows.png", width =4000, height = 4000, units = "px", res = 200)

plot_grid(zi192ef96.sizePlot, zi192ef81.sizePlot, zi418ef43.sizePlot, ncol = 1, 
          labels = c("Zi192 x Ef96", "Zi192 x Ef81", "Zi418 x Ef43"),  label_size = 25)

dev.off()

###################################################################
#size plotting male only
###################################################################

zi192ef96.sizePlotMales <- ggplot(data = zi192ef96.size, aes(x=number, y=male, color=chr)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  scale_x_discrete(limits=c(chrlabel.19296.size),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  theme(panel.background = element_blank())  + 
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 

zi192ef81.sizePlotMales <- ggplot(data = zi192ef81.size, aes(x=number, y=male, color=chr)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  scale_x_discrete(limits=c(chrlabel.19296.size),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  theme(panel.background = element_blank())  + 
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 

#Why am I like this? What is going on? 
zi418ef43.sizePlotMales <- ggplot(data = zi418ef43.size, aes(x=number, y=male, color=chr)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) +
  scale_x_discrete(limits=c(chrlabel.19296.size),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  theme(panel.background = element_blank()) + 
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 


png("../Figures/sizeFst_malesOnly_5000wingows.png", width =4000, height = 4000, units = "px",res = 200)

plot_grid(zi192ef96.sizePlotMales, zi192ef81.sizePlotMales, zi418ef43.sizePlotMales, ncol = 1, 
          labels = c("Zi192 x Ef96", "Zi192 x Ef81", "Zi418 x Ef43"), label_size = 25)

dev.off()

###################################################################
#shape plotting
###################################################################

#I want to plot the sexes on top of one another. 

zi192ef96.shapePlot <- ggplot(data = zi192ef96.shape, aes(x=number, y=male, color=chr)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) +
  scale_x_discrete(limits=c(chrlabel.19296.shape),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  new_scale_colour() +
  geom_point(aes(y = female, col = chr), alpha = 0.08, show.legend = F,) +
  scale_colour_manual(values=c('tomato', 'coral3', 'tomato', 'coral3', 'tomato', 'coral3')) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  theme(panel.background = element_blank())  +
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20))

zi192ef81.shapePlot <- ggplot(data = zi192ef81.shape, aes(x=number, y=male, color=chr)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  scale_x_discrete(limits=c(chrlabel.19296.size),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  new_scale_colour() +
  geom_point(aes(y = female, col = chr), alpha = 0.08, show.legend = F,) +
  scale_colour_manual(values=c('tomato', 'coral3', 'tomato', 'coral3', 'tomato', 'coral3')) +
  theme(panel.background = element_blank())  + 
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 

#This looks better.  
zi418ef43.shapePlot <- ggplot(data = zi418ef43.shape, aes(x=number, y=male, color=chr)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) +
  scale_x_discrete(limits=c(chrlabel.19296.size),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  new_scale_colour() +
  geom_point(aes(y = female, col = chr), alpha = 0.08, show.legend = F,) +
  scale_colour_manual(values=c('tomato', 'coral3', 'tomato', 'coral3', 'tomato', 'coral3')) +
  theme(panel.background = element_blank())  + 
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 


png("../Figures/shapeFst_bySex_5000wingows.png", width =4000, height = 4000, units = "px",res = 200)
#zi192ef96.shapePlot
#"Zi192 x Ef96", 
plot_grid(zi192ef96.shapePlot, zi192ef81.shapePlot, zi418ef43.shapePlot, ncol = 1, 
          labels = c("Zi192 x Ef96", "Zi192 x Ef81", "Zi418 x Ef43"), label_size = 25)

dev.off()

###################################################################
#shape plotting, males only 
###################################################################

#I want to plot the sexes on top of one another. 

zi192ef96.shapePlotMales <- ggplot(data = zi192ef96.shape, aes(x=number, y=male, color=chr)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) +
  scale_x_discrete(limits=c(chrlabel.19296.size),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  theme(panel.background = element_blank())  +
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20))

zi192ef81.shapePlotMales <- ggplot(data = zi192ef81.shape, aes(x=number, y=male, color=chr)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  scale_x_discrete(limits=c(chrlabel.19296.size),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  theme(panel.background = element_blank())  + 
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 

#This looks better.  
zi418ef43.shapePlotMales <- ggplot(data = zi418ef43.shape, aes(x=number, y=male, color=chr)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) +
  scale_x_discrete(limits=c(chrlabel.19296.size),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  theme(panel.background = element_blank())  + 
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 


png("../Figures/shapeFst_malesOnly_5000windows.png", width =4000, height = 4000, units = "px",res = 200)
#zi192ef96.sizePlot
#"Zi192 x Ef96", 
plot_grid(zi192ef96.shapePlotMales, zi192ef81.shapePlotMales, zi418ef43.shapePlotMales, ncol = 1, 
          labels = c("Zi192 x Ef96","Zi192 x Ef81", "Zi418 x Ef43"), label_size = 25)

dev.off()

###################################################################
#shape plotting, females only 
###################################################################

zi192ef96.shapePlotFemales <- ggplot(data = zi192ef96.shape, aes(x=number, y=female, color=chr)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) +
  scale_x_discrete(limits=c(chrlabel.19296.size),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  theme(panel.background = element_blank())  +
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20))

zi192ef81.shapePlotFemales <- ggplot(data = zi192ef81.shape, aes(x=number, y=female, color=chr)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  scale_x_discrete(limits=c(chrlabel.19296.size),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  theme(panel.background = element_blank())  + 
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 

#This looks better.  
zi418ef43.shapePlotFemales <- ggplot(data = zi418ef43.shape, aes(x=number, y=female, color=chr)) +
  geom_point(size=1, show.legend = F, alpha = 0.5) +
  scale_x_discrete(limits=c(chrlabel.19296.size),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  theme(panel.background = element_blank())  + 
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 


png("../Figures/shapeFst_femalesOnly_5000windows.png", width =4000, height = 4000, units = "px",res = 200)
#zi192ef96.sizePlot
#"Zi192 x Ef96", 
plot_grid(zi192ef96.shapePlotFemales, zi192ef81.shapePlotFemales, zi418ef43.shapePlotFemales, ncol = 1, 
          labels = c("Zi192 x Ef96", "Zi192 x Ef81", "Zi418 x Ef43"), label_size = 25)

dev.off()

###################################################################
#extra plots to add to things (specifically for exit seminar)
###################################################################

#I want to plot what a shared and not shared archetecture would look like
#For shared, I will plot the 3R peak because that can be reused. 

zi192ef96.3R <- filter(zi192ef96.shape, chr == "3R")
zi192ef96.3R$Cross <- "zi192ef96"
zi192ef81.3R <- filter(zi192ef81.shape, chr == "3R")
zi192ef81.3R$Cross <- "zi192ef81"
zi418ef43.3R <- filter(zi418ef43.shape, chr == "3R")
zi418ef43.3R$Cross <- "zi418ef43"

all.3R <- rbind(zi192ef96.3R, zi192ef81.3R, zi418ef43.3R)
peak.3R <- rbind(zi192ef96.3R, zi192ef81.3R)

png("../Figures/sharedGeneticBasisEx.png", width =2000, height = 1000, units = "px",res = 200)
ggplot(peak.3R, aes(x = window, y = male, col = Cross)) + 
  geom_point(size=1, show.legend = F, alpha = 0.4) + 
  scale_colour_manual(values=c('red', 'black') ) +
  theme(panel.background = element_blank())  + 
  xlab("") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 
dev.off()

png("../Figures/3Rfst_overplotting.png", width =2000, height = 1000, units = "px",res = 200)
ggplot(all.3R, aes(x = window, y = male, col = Cross)) + 
  geom_point(size=1, alpha = 0.5) + 
  scale_colour_manual(values=c('grey25', 'grey46', 'black') ) +
  theme(panel.background = element_blank())  + 
  xlab("") +
  ylab("meanFst") +
  ylim(0, 0.6) +
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 
dev.off()


#zi192ef96 male and female size X is a good example of a not shared basis
no.share <- filter(zi192ef96.size.long, chr == "X")

png("../Figures/noaredBasisExample_0and1.png", width =2000, height = 1000, units = "px",res = 200)
ggplot(no.share, aes(x = window, y = meanFst, col = sex)) + 
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  scale_colour_manual(values=c('grey46', 'black') ) +
  theme(panel.background = element_blank())  + 
  xlab("") +
  ylab("meanFst") +
  ylim(0, 1) +
  geom_hline(yintercept = 1, col = "red") + 
  geom_hline(yintercept = 0, col = "red") + 
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 
dev.off()

png("../Figures/noaredBasisExample.png", width =2000, height = 1000, units = "px",res = 200)
ggplot(no.share, aes(x = window, y = meanFst, col = sex)) + 
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  scale_colour_manual(values=c('grey46', 'black') ) +
  theme(panel.background = element_blank())  + 
  xlab("") +
  ylab("meanFst") +
  ylim(0, 0.4) +
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 
dev.off()


###################################################################
#Size outlier windows
###################################################################
#hopefully Tyler can find his plotting function?
#Deciding to try chromoMap because I like the asthetics and vingette 

library(chromoMap)
#the first thing I need is the general start and stop of the chromosomes for this genome assembly 

#I need the start and stop of the genome. This should be easier to google than it is. Will ask for more help later. 

# library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
# gr <- genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
# test<-seqinfo(gr)[c("chr2L","chr2R","chr3L","chr3R","chr4","chrX")]
#going to make this into a ref file with tabs. 

#I don't know what the problem is here. It reads in fine and it looks the way I expect it to without extra characters. 
chr.file="dmel_v6_chrLimits.txt"
chr.data <- read.delim("dmel_v6_chrLimits.txt", header = FALSE)

##the problem is that these elements dont have unique identifiers
#this is female data. 
ann.file="Spen2022SizeOutlier_v6.txt"
pool.ann.data <- read.delim("Spen2022SizeOutlier_v6.txt", header = FALSE)

names(pool.ann.data) <- c("feat", "chr", "start", "stop")
pool.ann.data$type <- "Pool"

#Other than the weird numbers on the side, this doesn't suck 
chromoMap(list(chr.data),list(pool.ann.data),
          chr.2D.plot = T)


#Now to create this for my crosses. first, I need to pull the outlier regions of my chromosomes. Going to use mean +3sd 

head(zi192ef96.size)
zi192ef96.size.mean <- mean(zi192ef96.size$male)
zi192ef96.size.sd <- sd(zi192ef96.size$male)
zi192ef96.size.cutoff <- zi192ef96.size.mean + 3*sd(zi192ef96.size$male)

ggplot(zi192ef96.size, aes(x = number, y = male, col = chr)) + 
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  scale_x_discrete(limits=c(chrlabel.19296.size),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  theme(panel.background = element_blank())  + 
  xlab("") +
  ylab("meanFst") +
  geom_hline(yintercept = zi192ef96.size.cutoff, alpha = 0.5) +
  ylim(0, 0.6) +
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 

zi192ef96.size.outliers <- filter(zi192ef96.size, male >= zi192ef96.size.cutoff)

dim(zi192ef96.size.outliers)

#now I want to condense this into smaller chunks. I've done this before with bumphunter. 
library(bumphunter)

zi192ef96.size.outliers$Lwindow <- zi192ef96.size.outliers$window - 2500
zi192ef96.size.outliers$Rwindow <- zi192ef96.size.outliers$window + 2500

zi192ef96.size.c1 <- clusterMaker(zi192ef96.size.outliers$chr, zi192ef96.size.outliers$window, maxGap = 100000)
table(c1)
100000/5000 #max 12 windows between to be called different. 

zi192ef96.size.outliers$group <- zi192ef96.size.c1

#Now I need to make a table with the start and end of each group.

max(zi192ef96.size.outliers$group)
#make somethign to fill in. 
zi192ef96.size.widows <- data.frame(matrix(NA, nrow = 45, ncol = 4))
names(zi192ef96.size.widows) <- c("feat", "chr", "start", "stop")

for (i in 1:45) {
  zi192ef96.size.widows[i,1] <- paste0("zi192ef96_", i)
  temp <- filter(zi192ef96.size.outliers, group == i)
  zi192ef96.size.widows[i,2] <- temp$chr[1]
  zi192ef96.size.widows[i,3] <- min(temp$Lwindow)
  zi192ef96.size.widows[i,4] <- min(temp$Rwindow)
}

zi192ef96.size.widows$type <- "Zi192Ef96"

############################################
zi192ef81.size.mean <- mean(zi192ef81.size$male)
zi192ef81.size.sd <- sd(zi192ef81.size$male)
zi192ef81.size.cutoff <- zi192ef81.size.mean + 3*sd(zi192ef81.size$male)

ggplot(zi192ef81.size, aes(x = number, y = male, col = chr)) + 
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  scale_x_discrete(limits=c(chrlabel.19296.size),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  theme(panel.background = element_blank())  + 
  xlab("") +
  ylab("meanFst") +
  geom_hline(yintercept = zi192ef81.size.cutoff, alpha = 0.5) +
  ylim(0, 0.6) +
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 

zi192ef81.size.outliers <- filter(zi192ef81.size, male >= zi192ef81.size.cutoff)

dim(zi192ef81.size.outliers)

#now I want to condense this into smaller chunks.

zi192ef81.size.outliers$Lwindow <- zi192ef81.size.outliers$window - 2500
zi192ef81.size.outliers$Rwindow <- zi192ef81.size.outliers$window + 2500

zi192ef81.size.c1 <- clusterMaker(zi192ef81.size.outliers$chr, zi192ef81.size.outliers$window, maxGap = 100000)
table(c1)
100000/5000 #max 20 windows between to be called different. 

zi192ef81.size.outliers$group <- zi192ef81.size.c1

#Now I need to make a table with the start and end of each group.

max(zi192ef81.size.outliers$group)
#make somethign to fill in. 
zi192ef81.size.widows <- data.frame(matrix(NA, nrow = 37, ncol = 4))
names(zi192ef81.size.widows) <- c("feat", "chr", "start", "stop")

for (i in 1:37) {
  zi192ef81.size.widows[i,1] <- paste0("zi192ef81_", i)
  temp <- filter(zi192ef81.size.outliers, group == i)
  zi192ef81.size.widows[i,2] <- temp$chr[1]
  zi192ef81.size.widows[i,3] <- min(temp$Lwindow)
  zi192ef81.size.widows[i,4] <- min(temp$Rwindow)
}

zi192ef81.size.widows$type <- "Zi192Ef81"

##############################################
zi418ef43.size.mean <- mean(zi418ef43.size$male)
zi418ef43.size.sd <- sd(zi418ef43.size$male)
zi418ef43.size.cutoff <- zi418ef43.size.mean + 3*sd(zi418ef43.size$male)

ggplot(zi418ef43.size, aes(x = number, y = male, col = chr)) + 
  geom_point(size=1, show.legend = F, alpha = 0.5) + 
  scale_x_discrete(limits=c(chrlabel.19296.size),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  theme(panel.background = element_blank())  + 
  xlab("") +
  ylab("meanFst") +
  geom_hline(yintercept = zi418ef43.size.cutoff, alpha = 0.5) +
  ylim(0, 0.6) +
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 

zi418ef43.size.outliers <- filter(zi418ef43.size, male >= zi418ef43.size.cutoff)

dim(zi418ef43.size.outliers)

#now I want to condense this into smaller chunks.

zi418ef43.size.outliers$Lwindow <- zi418ef43.size.outliers$window - 2500
zi418ef43.size.outliers$Rwindow <- zi418ef43.size.outliers$window + 2500

zi418ef43.size.c1 <- clusterMaker(zi418ef43.size.outliers$chr, zi418ef43.size.outliers$window, maxGap = 100000)
table(c1)
100000/5000 #max 20 windows between to be called different. 

zi418ef43.size.outliers$group <- zi418ef43.size.c1

#Now I need to make a table with the start and end of each group.

max(zi418ef43.size.outliers$group)
#make somethign to fill in. 
zi418ef43.size.widows <- data.frame(matrix(NA, nrow = 52, ncol = 4))
names(zi418ef43.size.widows) <- c("feat", "chr", "start", "stop")

for (i in 1:52) {
  zi418ef43.size.widows[i,1] <- paste0("zi418ef43_", i)
  temp <- filter(zi418ef43.size.outliers, group == i)
  zi418ef43.size.widows[i,2] <- temp$chr[1]
  zi418ef43.size.widows[i,3] <- min(temp$Lwindow)
  zi418ef43.size.widows[i,4] <- min(temp$Rwindow)
}

zi418ef43.size.widows$type <- "zi418ef43"


##############################################

anno.dat <- rbind(pool.ann.data, zi192ef96.size.widows, zi192ef81.size.widows, zi418ef43.size.widows)

#Now I want to add this to my previous plot.
#I don't know if I love this but it exists.

#I need to figure this out still. 
#I am almost 100% sure that the Pool lab made their plots without code. 
chromoMap(list(chr.data),list(anno.dat),
          chr.2D.plot = T,
          plot_filter = list(c("col","byCategory")),
          ch2D.colors = c("gray", "blue", "green", "red"), 
          chr_width = 10,
          canvas_height = 650, 
          canvas_width = 650, left_margin = 110, 
          text_font_size = 20,
          lg_x = -400,
          export.options = T)

#I could just make this a shiny that people can click on 
#Maybe for a later Katie problem. 
# library(shiny)
# 
# # Define UI for application that draws chromoMap
# ui <- fluidPage(
#   
#   # Application title
#   titlePanel("Size Mapping QTL overlap"),
#   
#   # you can use GUI controls for your chromoMap
#   #sidebarLayout(
#    # sidebarPanel(
#       #some code
#    # ),
#     
#     # Show a plot of the generated distribution
#     mainPanel(
#       chromoMapOutput("myChromoMap")
#     )
#   )
# )
# 
# # Define server logic required to draw chromoMap
# server <- function(input, output) {
#   
#   output$myChromoMap <- renderChromoMap({
#     chromoMap(list(chr.data),list(anno.dat),
#               chr.2D.plot = T,
#               plot_filter = list(c("col","byCategory")),
#               ch2D.colors = c("red", "blue", "green", "gray"))
#   })
# }
# 
# 
# # Run the application 
# shinyApp(ui = ui, server = server)

#######################################################################
#I do want to get the %overlap. 

zi418ef43.size.outliers$chr.pos <- paste(zi418ef43.size.outliers$chr, zi418ef43.size.outliers$window, sep = "_")
zi192ef96.size.outliers$chr.pos <- paste(zi192ef96.size.outliers$chr, zi192ef96.size.outliers$window, sep = "_")
zi192ef81.size.outliers$chr.pos <- paste(zi192ef81.size.outliers$chr, zi192ef81.size.outliers$window, sep = "_")

#8
sum(zi418ef43.size.outliers$chr.pos %in% zi192ef96.size.outliers$chr.pos)

#4
sum(zi418ef43.size.outliers$chr.pos %in% zi192ef81.size.outliers$chr.pos)

#112
sum(zi192ef96.size.outliers$chr.pos %in% zi192ef81.size.outliers$chr.pos)

#0
sum(zi192ef96.size.outliers$chr.pos %in% zi192ef81.size.outliers$chr.pos %in% zi418ef43.size.outliers$chr.pos)

#306
nrow(zi418ef43.size.outliers)
#421
nrow(zi192ef96.size.outliers)
#552
nrow(zi192ef81.size.outliers)



###################################################################
#Shape outlier windows
###################################################################

zi192ef96.shape.mean <- mean(zi192ef96.shape$male)
zi192ef96.shape.sd <- sd(zi192ef96.shape$male)
zi192ef96.shape.cutoff <- zi192ef96.shape.mean + 3*sd(zi192ef96.shape$male)

ggplot(zi192ef96.shape, aes(x = number, y = male, col = chr)) + 
  geom_point(shape=1, show.legend = F, alpha = 0.5) + 
  scale_x_discrete(limits=c(chrlabel.19296.shape),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  theme(panel.background = element_blank())  + 
  xlab("") +
  ylab("meanFst") +
  geom_hline(yintercept = zi192ef96.shape.cutoff, alpha = 0.5) +
  ylim(0, 0.6) +
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 

zi192ef96.shape.outliers <- filter(zi192ef96.shape, male >= zi192ef96.shape.cutoff)

dim(zi192ef96.shape.outliers)

#now I want to condense this into smaller chunks. I've done this before with bumphunter. 
library(bumphunter)

zi192ef96.shape.outliers$Lwindow <- zi192ef96.shape.outliers$window - 2500
zi192ef96.shape.outliers$Rwindow <- zi192ef96.shape.outliers$window + 2500

zi192ef96.shape.c1 <- clusterMaker(zi192ef96.shape.outliers$chr, zi192ef96.shape.outliers$window, maxGap = 100000)
table(zi192ef96.shape.c1)
100000/5000 #max 12 windows between to be called different. 

zi192ef96.shape.outliers$group <- zi192ef96.shape.c1

#Now I need to make a table with the start and end of each group.

max(zi192ef96.shape.outliers$group)
#make somethign to fill in. 
zi192ef96.shape.widows <- data.frame(matrix(NA, nrow = 47, ncol = 4))
names(zi192ef96.shape.widows) <- c("feat", "chr", "start", "stop")

for (i in 1:47) {
  zi192ef96.shape.widows[i,1] <- paste0("zi192ef96_", i)
  temp <- filter(zi192ef96.shape.outliers, group == i)
  zi192ef96.shape.widows[i,2] <- temp$chr[1]
  zi192ef96.shape.widows[i,3] <- min(temp$Lwindow)
  zi192ef96.shape.widows[i,4] <- min(temp$Rwindow)
}

zi192ef96.shape.widows$type <- "Zi192Ef96"

############################################
zi192ef81.shape.mean <- mean(zi192ef81.shape$male)
zi192ef81.shape.sd <- sd(zi192ef81.shape$male)
zi192ef81.shape.cutoff <- zi192ef81.shape.mean + 3*sd(zi192ef81.shape$male)

ggplot(zi192ef81.shape, aes(x = number, y = male, col = chr)) + 
  geom_point(shape=1, show.legend = F, alpha = 0.5) + 
  scale_x_discrete(limits=c(chrlabel.19296.shape),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  theme(panel.background = element_blank())  + 
  xlab("") +
  ylab("meanFst") +
  geom_hline(yintercept = zi192ef81.shape.cutoff, alpha = 0.5) +
  ylim(0, 0.6) +
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 

zi192ef81.shape.outliers <- filter(zi192ef81.shape, male >= zi192ef81.shape.cutoff)

dim(zi192ef81.shape.outliers)

#now I want to condense this into smaller chunks.

zi192ef81.shape.outliers$Lwindow <- zi192ef81.shape.outliers$window - 2500
zi192ef81.shape.outliers$Rwindow <- zi192ef81.shape.outliers$window + 2500

zi192ef81.shape.c1 <- clusterMaker(zi192ef81.shape.outliers$chr, zi192ef81.shape.outliers$window, maxGap = 100000)
table(zi192ef81.shape.c1)
100000/5000 #max 20 windows between to be called different. 

zi192ef81.shape.outliers$group <- zi192ef81.shape.c1

#Now I need to make a table with the start and end of each group.

max(zi192ef81.shape.outliers$group)
#make somethign to fill in. 
zi192ef81.shape.widows <- data.frame(matrix(NA, nrow = 14, ncol = 4))
names(zi192ef81.shape.widows) <- c("feat", "chr", "start", "stop")

for (i in 1:14) {
  zi192ef81.shape.widows[i,1] <- paste0("zi192ef81_", i)
  temp <- filter(zi192ef81.shape.outliers, group == i)
  zi192ef81.shape.widows[i,2] <- temp$chr[1]
  zi192ef81.shape.widows[i,3] <- min(temp$Lwindow)
  zi192ef81.shape.widows[i,4] <- min(temp$Rwindow)
}

zi192ef81.shape.widows$type <- "Zi192Ef81"

##############################################
zi418ef43.shape.mean <- mean(zi418ef43.shape$male)
zi418ef43.shape.sd <- sd(zi418ef43.shape$male)
zi418ef43.shape.cutoff <- zi418ef43.shape.mean + 3*sd(zi418ef43.shape$male)

ggplot(zi418ef43.shape, aes(x = number, y = male, col = chr)) + 
  geom_point(shape=1, show.legend = F, alpha = 0.5) + 
  scale_x_discrete(limits=c(chrlabel.19296.shape),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  theme(panel.background = element_blank())  + 
  xlab("") +
  ylab("meanFst") +
  geom_hline(yintercept = zi418ef43.shape.cutoff, alpha = 0.5) +
  ylim(0, 0.6) +
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20)) 

zi418ef43.shape.outliers <- filter(zi418ef43.shape, male >= zi418ef43.shape.cutoff)

dim(zi418ef43.shape.outliers)

#now I want to condense this into smaller chunks.

zi418ef43.shape.outliers$Lwindow <- zi418ef43.shape.outliers$window - 2500
zi418ef43.shape.outliers$Rwindow <- zi418ef43.shape.outliers$window + 2500

zi418ef43.shape.c1 <- clusterMaker(zi418ef43.shape.outliers$chr, zi418ef43.shape.outliers$window, maxGap = 100000)
table(c1)
100000/5000 #max 20 windows between to be called different. 

zi418ef43.shape.outliers$group <- zi418ef43.shape.c1

#Now I need to make a table with the start and end of each group.

max(zi418ef43.shape.outliers$group)
#make somethign to fill in. 
zi418ef43.shape.widows <- data.frame(matrix(NA, nrow = 9, ncol = 4))
names(zi418ef43.shape.widows) <- c("feat", "chr", "start", "stop")

for (i in 1:9) {
  zi418ef43.shape.widows[i,1] <- paste0("zi418ef43_", i)
  temp <- filter(zi418ef43.shape.outliers, group == i)
  zi418ef43.shape.widows[i,2] <- temp$chr[1]
  zi418ef43.shape.widows[i,3] <- min(temp$Lwindow)
  zi418ef43.shape.widows[i,4] <- min(temp$Rwindow)
}

zi418ef43.shape.widows$type <- "zi418ef43"


##############################################

anno.dat.shape <- rbind(zi192ef96.shape.widows, zi192ef81.shape.widows, zi418ef43.shape.widows)

#Now I want to add this to my previous plot.
#I don't know if I love this but it exists.

#I need to figure this out still. 
#I am almost 100% sure that the Pool lab made their plots without code. 
chromoMap(list(chr.data),list(anno.dat.shape),
          chr.2D.plot = T,
          plot_filter = list(c("col","byCategory")),
          ch2D.colors = c("blue", "green", "red"), 
          chr_width = 10,
          canvas_height = 650, 
          canvas_width = 650, left_margin = 110, 
          text_font_size = 20,
          lg_x = -400,
          export.options = T)

#######################################################################
#I do want to get the %overlap. 

zi418ef43.shape.outliers$chr.pos <- paste(zi418ef43.shape.outliers$chr, zi418ef43.shape.outliers$window, sep = "_")
zi192ef96.shape.outliers$chr.pos <- paste(zi192ef96.shape.outliers$chr, zi192ef96.shape.outliers$window, sep = "_")
zi192ef81.shape.outliers$chr.pos <- paste(zi192ef81.shape.outliers$chr, zi192ef81.shape.outliers$window, sep = "_")

#0
sum(zi418ef43.shape.outliers$chr.pos %in% zi192ef96.shape.outliers$chr.pos)

#0
sum(zi418ef43.shape.outliers$chr.pos %in% zi192ef81.shape.outliers$chr.pos)

#257
sum(zi192ef96.shape.outliers$chr.pos %in% zi192ef81.shape.outliers$chr.pos)


zi192ef96.shape.outliers[which(zi192ef96.shape.outliers$chr.pos %in% zi192ef81.shape.outliers$chr.pos),]

#0
sum(zi192ef96.shape.outliers$chr.pos %in% zi192ef81.shape.outliers$chr.pos %in% zi418ef43.shape.outliers$chr.pos)

#700
nrow(zi418ef43.shape.outliers)
#366
nrow(zi192ef96.shape.outliers)
#722
nrow(zi192ef81.shape.outliers)
