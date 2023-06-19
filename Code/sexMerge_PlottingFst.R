
library(data.table)
library(tidyverse)
source("KP_genomescan_source.R")

zi192ef96 <- fread("../Data/sexmerge_z192ef96_5000windows.fst")
zi192ef96 <- zi192ef96[,c(1:5,6,11)]

#where the fuck is that one population and why do I suck so bad. 
zi192ef96 <-populaionFst_cleanup(zi192ef96, x = c("shape", "size"))

fixed_zi192ef96 <- chrNumbering(zi192ef96)

chrlabel <- middleChr(fixed_zi192ef96)

zi192ef96sizePlot <- ggplot(data = fixed_zi192ef96, aes(x=number, y=size, color=chr)) + 
  geom_point(size=1, show.legend = F, alpha = 0.6) + 
  theme(panel.background = element_blank()) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.9) + 
  scale_x_discrete(limits=c(chrlabel),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20), 
        axis.text.y= element_text(size=20))

png("../Figures/zi192ef96_sexMergeSize_5000bp.png", width = 1060, height = 412, units = "px")
zi192ef96sizePlot
dev.off()

zi192ef96shapePlot <- ggplot(data = fixed_zi192ef96, aes(x=number, y=shape, color=chr)) +
  geom_point(size=1, show.legend = F, alpha = 0.6) +
  theme(panel.background = element_blank()) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.9) +
  scale_x_discrete(limits=c(chrlabel),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20))

png("../Figures/zi192ef96_sexMergeShape_5000bp.png", width = 1060, height = 412, units = "px")
zi192ef96shapePlot
dev.off()

#this one is missing for size. 
zi192ef81 <- fread("../Data/sexmerge_z192ef81shape_5000windows.fst")
#zi192ef81 <- zi192ef81[,c(1:5,8)]


zi192ef81 <-populaionFst_cleanup(zi192ef81, x = c('shape'))

fixed_zi192ef81 <- chrNumbering(zi192ef81)

chrlabel <- middleChr(zi192ef81)

zi192ef81shapePlot <- ggplot(data = fixed_zi192ef81, aes(x=number, y = shape, color=chr)) +
  geom_point(size=1, show.legend = F, alpha = 0.6) +
  theme(panel.background = element_blank()) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.9) +
  scale_x_discrete(limits=c(chrlabel),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20),
        axis.text.y= element_text(size=20))

png("../Figures/zi192ef96_shape_5000bp.png", width = 1060, height = 412, units = "px")
zi192ef81shapePlot
dev.off()

#add this in? talk to ID?

# zi192ef81MsizePlot <- ggplot(data = fixed_zi192ef81, aes(x=number, y = sizeM, color=chr)) +
#   geom_point(size=1, show.legend = F, alpha = 0.6) +
#   theme(panel.background = element_blank()) +
#   xlab("Chromosome") +
#   ylab("meanFst") +
#   ylim(0, 0.9) +
#   scale_x_discrete(limits=c(chrlabel),
#                    labels = c("X","2L", "2R", '3L', '3R', '4')) +
#   scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
#   theme(text = element_text(size=25),
#         axis.text.x= element_text(size=20),
#         axis.text.y= element_text(size=20))
# 
# png("zi192ef81_sizeM_1000bp.png", width = 1060, height = 412, units = "px")
# zi192ef81MsizePlot
# dev.off()


zi418ef43 <- fread("../Data/sexmerge_zi418ef43_withShape2.fst")
zi418ef43 <- zi418ef43[,c(1:5, 7, 10)]

zi418ef43 <-populaionFst_cleanup(zi418ef43, x = c('shape', 'size'))

fixed_zi418ef43 <- chrNumbering(zi418ef43)

chrlabel <- middleChr(fixed_zi418ef43)

zi418ef43shapePlot <- ggplot(data = fixed_zi418ef43, aes(x=number, y=shape, color=chr)) + 
  geom_point(size=1, show.legend = F, alpha = 0.6) + 
  theme(panel.background = element_blank()) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.9) + 
  scale_x_discrete(limits=c(chrlabel),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20), 
        axis.text.y= element_text(size=20))

png("../Figures/zi418ef43_sexMergeShape_5000bp.png", width = 1060, height = 412, units = "px")
zi418ef43shapePlot
dev.off()

zi418ef43sizePlot <- ggplot(data = fixed_zi418ef43, aes(x=number, y=size, color=chr)) + 
  geom_point(size=1, show.legend = F, alpha = 0.6) + 
  theme(panel.background = element_blank()) +
  xlab("Chromosome") +
  ylab("meanFst") +
  ylim(0, 0.9) +
  scale_x_discrete(limits=c(chrlabel),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  theme(text = element_text(size=25),
        axis.text.x= element_text(size=20), 
        axis.text.y= element_text(size=20))

png("../Figures/zi418ef43_sexMergeSize_5000bp.png", width = 1060, height = 412, units = "px")
zi418ef43sizePlot
dev.off()

# zi418ef43sizeFPlot <- ggplot(data = fixed_zi418ef43, aes(x=number, y=sizeF, color=chr)) + 
#   geom_point(size=1, show.legend = F, alpha = 0.6) + 
#   theme(panel.background = element_blank()) +
#   xlab("Chromosome") +
#   ylab("meanFst") +
#   ylim(0, 0.9) +
#   scale_x_discrete(limits=c(chrlabel),
#                    labels = c("X","2L", "2R", '3L', '3R', '4')) +
#   scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
#   theme(text = element_text(size=25),
#         axis.text.x= element_text(size=20), 
#         axis.text.y= element_text(size=20))
# 
# png("zi418ef43_sizeF_5000bp.png", width = 1060, height = 412, units = "px")
# zi418ef43sizeFPlot
# dev.off()

####################
library(cowplot)
#adding plots together for poster. 


png("shapeM_5000bp.png", width = 1060, height = (3*412), units = "px")
plot_grid(zi192ef96shapeMPlot, zi192ef81MshapePlot, zi418ef43shapeMPlot,
          ncol = 1, label_size = 25,
          labels = c("ZI192 x EF96", "ZI192 x EF81", "ZI418 x EF43"))
dev.off()

png("sizeM_5000bp.png", width = 1060, height = (3*412), units = "px")
plot_grid(zi192ef96sizeMPlot, zi192ef81MsizePlot, zi418ef43sizeMPlot,
          ncol = 1, label_size = 25,
          labels = c("ZI192 x EF96", "ZI192 x EF81", "ZI418 x EF43"))
dev.off()



#########3R#####################
zi418ef433R <- filter(fixed_zi418ef43, chr == "3R")

quantile(zi418ef433R$shapeM)

#this tiny peak is not the same as the 3R peak for the other pop. 
crap <- filter(zi418ef433R, shapeM >= 0.06)


#####I want to make BED files to look at the intersection between regions of high differentation (and include pool data for size)####


####Just M size for now because I need to triage my life##### 

zi192ef81_sizeRegions <- filter(fixed_zi192ef81, sizeM >= 0.15)

zi192ef96_sizeRegions <- filter(fixed_zi192ef96, sizeM >= 0.15)

zi418ef43_sizeRegions <- filter(fixed_zi418ef43, sizeM >= 0.15)


#A bed file has the chromosome then feature start and feature end so, here window +/- 

length(intersect(zi192ef81_sizeRegions$window, 
                 zi192ef96_sizeRegions$window))/ mean(length(zi192ef96_sizeRegions$window),length(zi192ef81_sizeRegions$window))


length(intersect(zi192ef81_sizeRegions$window, 
                 zi418ef43_sizeRegions$window))/  mean(length(zi418ef43_sizeRegions$window),length(zi192ef81_sizeRegions$window))

length(intersect(zi192ef96_sizeRegions$window, 
                 zi418ef43_sizeRegions$window))/ mean(length(zi418ef43_sizeRegions$window),length(zi192ef96_sizeRegions$window))

#Now to subset 

int9681 <- zi192ef81_sizeRegions[zi192ef81_sizeRegions$window %in%
                                    intersect(zi192ef81_sizeRegions$window, 
                                              zi192ef96_sizeRegions$window),]
int9643 <- zi192ef96_sizeRegions[zi192ef96_sizeRegions$window %in%
                                   intersect(zi192ef96_sizeRegions$window, 
                                             zi418ef43_sizeRegions$window),]

int8143 <- zi192ef81_sizeRegions[zi192ef81_sizeRegions$window %in%
                                   intersect(zi192ef81_sizeRegions$window, 
                                             zi418ef43_sizeRegions$window),]

length(intersect(zi192ef96_sizeRegions$window, 
                 intersect(zi192ef81_sizeRegions$window,
                 zi418ef43_sizeRegions$window)))


intall <- zi192ef81_sizeRegions[zi192ef81_sizeRegions$window %in%
                                   intersect(zi192ef96_sizeRegions$window, 
                                             intersect(zi192ef81_sizeRegions$window,
                                                       zi418ef43_sizeRegions$window)),]
