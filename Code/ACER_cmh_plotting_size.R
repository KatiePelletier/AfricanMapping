library(data.table)
library(tidyverse)

source('KP_genomescan_source.R')

size.cmh <- fread("allSizeMales_cmh_acer_catGenome.csv")

#reordering the genome so the X is first and adding numbers for plotting. 

size.cmh.fix <- chrNumbering(size.cmh)

size.cmh.middle <- middleChr(size.cmh)

#need -logp for plotting.

#don't know why it is like this but it is
#I think it has to do with how I wrote then read this in. 
size.cmh.fix$logp <- -log(as.numeric(size.cmh.fix$pval2L))


#filtering the sig sites for plotting. 
hits <- filter(size.cmh.fix, padj2L <= 0.05)
hits
#107025
#lololol. So many sites still. 
#going to just make this a line. 
nrow(hits)

#I also realized this is calculated within chr. which is nice and all but I want the whole genome at once. 
size.cmh.fix$All.fdr <- p.adjust(size.cmh.fix$pval2L, "fdr")
hits2 <- filter(size.cmh.fix, All.fdr <= 0.01)
#46802. going with this. Its unsuprisng that this is polygenic.  
nrow(hits2)





cmhPlot <- ggplot(data = size.cmh.fix, aes(x=number, y=logp, color=chr)) + 
  geom_point(size=1, show.legend = F, alpha = 0.6) + 
  theme(panel.background = element_blank()) +
  xlab("Chromosome") +
  ylab("-log(p-val)") +
  scale_colour_manual(values=c('black', 'grey46', 'black', 'grey46', 'black','grey46')) +
  #geom_point(data = hits2, aes(x=number, y=logp), col = "red")+
  scale_x_discrete(limits=c(size.cmh.middle),
                   labels = c("X","2L", "2R", '3L', '3R', '4')) +
  theme(text = element_text(size=15),
        axis.text.x= element_text(size=12), 
        axis.text.y= element_text(size=12))


cmhPlot3 <- cmhPlot + 
  geom_hline(yintercept = min(hits2$logp), color = "red", 
             linetype = "dashed", linewidth= 1.5, alpha = 0.5) 


#I like this plot the best. 
png("size_males_cmh_ACER_zeroGen_withFDR_line.png", width=1060,height=412,units="px")
cmhPlot3
dev.off()

#I want to write out the hit table so I can annotate it with genes. 

write.csv(hits2, "size_males_hits_FDR001_forAnnotation.csv", row.names = FALSE, quote = FALSE)

#####Quick GO analysis using TopGO 

library(topGO)
library(org.Dm.eg.db)
library(GenomicFeatures)

#Reading in the annotated sig sites table 
#this adds an extra col with nothing but the data looks just fine. 
size.sites <- fread("size_males_hits_FDR001_genes.bed")

#first question. how many unique genes 
#this is specifically called because unique is masked by one of the other packages
#2769 honestly, fewer than I would have expected. 
length(unlist(base::unique(size.sites[,9])))

size.sites <- data.frame(size.sites)
names(size.sites) <- c("chr", "pos", "end", "FDR", "chr2", "geneStart", "geneStop", "FBID", "gene")


size.siggenes <- size.sites$FBID

#vector of genes. Hopefully the duplicates don't mess this up. I dont think they will. 
size.siggenes

#loading in the GO index we made 
gene_GO <- readMappings("fly_to_GO.delim")

#I need to make a list of all possible genes (all genes in fly) and then classify if these are in my data set or not. 
allgenes <- data.frame(names(gene_GO))
names(allgenes) <- "gene_id"

#this is all the genes that are near diff sites 
allgenes$diff <- as.numeric(allgenes$gene_id %in% size.siggenes)


all.coded <- as.factor(allgenes$diff)
names(all.coded) <- allgenes$gene_id
head(all.coded)

diff.sites <- as.factor(rep(1, length(size.siggenes)))
names(diff.sites) <- size.siggenes
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
size.termsTop20 <-GenTable(allgenes, classic = resultFisher, ranksOf = "classic", topNodes = 20)


size.termsTop50 <-GenTable(allgenes, classic = resultFisher, ranksOf = "classic", topNodes = 50)

write.csv(shape.termsTop50, file = "sizeGOanalysis_CMH_FDR0001_top50Terms.csv")





