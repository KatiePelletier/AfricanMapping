library(data.table)
library(tidyverse)

source('KP_genomescan_source.R')

shape.cmh <- fread("allShapeMales_cmh_acer_catGenome.csv")

#reordering the genome so the X is first and adding numbers for plotting. 

shape.cmh.fix <- chrNumbering(shape.cmh)

shape.cmh.middle <- middleChr(shape.cmh)

#need -logp for plotting.

#don't know why it is like this but it is
#I think it has to do with how I wrote then read this in. 
shape.cmh.fix$logp <- -log(as.numeric(shape.cmh.fix$pval2L))


#filtering the sig sites for plotting. 
hits <- filter(shape.cmh.fix, padj2L <= 0.05)
hits
#38768
#lololol. So many sites still. 
#going to just make this a line. 
nrow(hits)

#I also realized this is calculated within chr. which is nice and all but I want the whole genome at once. 
shape.cmh.fix$All.fdr <- p.adjust(shape.cmh.fix$pval2L, "fdr")
hits2 <- filter(shape.cmh.fix, All.fdr <= 0.01)
#20224. going with this. Its unsuprisng that this is polygenic.  
nrow(hits2)





cmhPlot <- ggplot(data = shape.cmh.fix, aes(x=number, y=logp, color=chr)) + 
  geom_point(size=1, show.legend = F, alpha = 0.6) + 
  theme(panel.background = element_blank()) +
  xlab("Chromosome") +
  ylab("-log(p-val)") +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  #geom_point(data = hits2, aes(x=number, y=logp), col = "red")+
  scale_x_discrete(limits=c(shape.cmh.middle),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  theme(text = element_text(size=15),
        axis.text.x= element_text(size=12), 
        axis.text.y= element_text(size=12))


cmhPlot3 <- cmhPlot + 
  geom_hline(yintercept = min(hits2$logp), color = "red", 
             linetype = "dashed", linewidth= 1.5, alpha = 0.5) 


#I like this plot the best. 
png("shape_males_cmh_ACER_zeroGen_withFDR_line.png", width=1060,height=412,units="px")
cmhPlot3
dev.off()

#I want to write out the hit table so I can annotate it with genes. 

write.csv(hits2, "shape_males_hits_FDR001_forAnnotation.csv", row.names = FALSE, quote = FALSE)

#####Quick GO analysis using TopGO 

library(topGO)
library(org.Dm.eg.db)
library(GenomicFeatures)

#Reading in the annotated sig sites table 
#this adds an extra col with nothing but the data looks just fine. 
shape.sites <- fread("shape_males_hits_FDR001_genes.bed")

#first question. how many unique genes 
#this is specifically called because unique is masked by one of the other packages
#6395 lolol
length(unlist(base::unique(shape.sites[,9])))

#45
sum(unlist(shape.sites[,9]) == "loco")

#shape.sites[unlist(shape.sites[,9] == "loco",)]

#20
sum(unlist(shape.sites[,9]) == "wge")

#3
sum(unlist(shape.sites[,9]) == "Takl2")

#86
sum(unlist(shape.sites[,9]) == "wake")

#shape.sites[shape.sites[,9] == "wake",]


#51
sum(unlist(shape.sites[,9]) == "Efa6")

#just one. thats cool 
sum(unlist(shape.sites[,9]) == "btn")

#all genes with sig sites in my ROI
shape.sites <- data.frame(shape.sites)
names(shape.sites) <- c("chr", "pos", "end", "FDR", "chr2", "geneStart", "geneStop", "FBID", "gene")

shape.sites.ROI <- (shape.sites %>%
  filter(chr == "3R") %>%
    filter(pos > 22580970) %>%
    filter(pos < 22733819)
  )

base::unique(shape.sites.ROI$gene)


shape.sites.ROI[shape.sites.ROI$FDR == max(shape.sites.ROI$FDR), ]

shape.siggenes <- shape.sites$FBID


#vector of genes. Hopefully the duplicates don't mess this up. I dont think they will. 
shape.siggenes

#loading in the GO index we made 
gene_GO <- readMappings("fly_to_GO.delim")


#I need to make a list of all possible genes (all genes in fly) and then classify if these are in my data set or not. 
allgenes <- data.frame(names(gene_GO))
names(allgenes) <- "gene_id"

#this is all the genes that are near diff sites 
allgenes$diff <- as.numeric(allgenes$gene_id %in% shape.siggenes)


all.coded <- as.factor(allgenes$diff)
names(all.coded) <- allgenes$gene_id
head(all.coded)

diff.sites <- as.factor(rep(1, length(shape.siggenes)))
names(diff.sites) <- shape.siggenes
diff.sites

#test for sites 
gene_filter <- function(allScore){
  return(allScore == 1)
}

#testing for enrichment
allgenes <- new("topGOdata",
                ontology = "BP", 
                allGenes = all.coded,
                annotationFun = annFUN.gene2GO, 
                gene2GO = gene_GO
)


resultFisher <- runTest(allgenes, algorithm = "classic", statistic = "fisher")

#So general and truly useless. 
shape.termsTop20 <-GenTable(allgenes, classic = resultFisher, ranksOf = "classic", topNodes = 20)


shape.termsTop50 <-GenTable(allgenes, classic = resultFisher, ranksOf = "classic", topNodes = 50)

write.csv(shape.termsTop50, file = "shapeGOanalysis_CMH_FDR0001_top50Terms.csv")


