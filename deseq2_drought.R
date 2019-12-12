#### Looking at differentially expressed genes 

## Want to use DESeq2 on RSEM output. 
library(colorspace)
library(tximport)
library(DESeq2)
library(readr)
library(tidyverse)
library(plyr)
library(dplyr)
library(zoo)
library(topGO)
# Use tximport to import rsem RAW FILES
# these DO NOT HAVE LOWLY EXPRESSED GENES REMOVED

# You will need some objects frim transcript_info.R
# source("transcript_info.R")

counts_dir <- as.character("./rsem_results/")

## sample names 
samples <- c("YWW-1", "YWW-2", "YWW1-3", "YWW1-4", "YD1-1", "YD1-2", "YD1-3",
              "YRW-1", "YRW-2", "YWW2-3", "YWW2-4", "YD2-1", "YD2-2","YD2-3",
              "YD2-4",
              "SWW-1", "SWW-2", "SWW1-3", "SWW1-4", "SD1-1", "SD1-2", "SD1-3", 
              "SRW-1", "SRW-2", "SWW2-3", "SWW2-4", "SD2-1","SD2-2" ,"SD2-3",
              "SD2-4","SD2-5")

files <- paste0("./rsem_results/",
                paste0(samples, ".genes.results") ) 
names(files) <- samples

txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE) 

txi.rsem_lenadj <- txi.rsem

## using the median transcript length for fpkm calculation
median_len <- apply(txi.rsem$length, 1, median) 
med_len <- cbind(median_len,median_len,median_len,median_len,median_len,median_len,
                 median_len, median_len, median_len, median_len, median_len, median_len, 
                 median_len, median_len, median_len, median_len, median_len, median_len, median_len,
                 median_len,median_len,median_len,median_len,median_len,median_len,
                 median_len, median_len, median_len, median_len, median_len, median_len)
colnames(med_len) <- colnames(txi.rsem$length)
txi.rsem_lenadj$length <- med_len

## import counts  into DESeq2.
## we have genes with length 0, because not identified in experiment...change to 1 
txi.rsem$length[txi.rsem$length == 0] <- 1
txi.rsem_lenadj$length[txi.rsem_lenadj$length == 0] <- 1
samples[c(1,2,8,9,16,17,23,24)] <- c("YWW1-1", "YWW1-2", "YWW2-1", "YWW2-2",
                                            "SWW1-1", "SWW1-2", "SWW2-1", "SWW2-2") 

reps <- as.factor(c("a", "a", "b", "b","a", "a", "b","a", "a", "b", "b","a", "a", "b", "b","a", "a", "b", "b",
          "a", "a", "b","a", "a", "b", "b", "a", "a", "b", "b","b"))
sampleTable <- data.frame(library = samples, 
                          ecocond = factor(substr(samples,1,nchar(samples)-2)),
                          condition = factor(c(substr(samples[1:31], 2, nchar(samples[1:31])-2))),
                          ecotype = factor(c(rep("Yukon", 15),rep("Shandong",16))),
                          reps = reps)
                       
sampleTable$ecocond

# we are not using a separate batch correction step because of this paper:
# Methods that remove batch effects while retaining group differences may
#lead to exaggerated confidence in downstream analyses
dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~reps + ecocond)

# check the dataframe
head(assay(dds))

## filter out lowly expressed genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

## correct for batch effect from library prep protocols
dds$batch <- factor(dds$reps, levels=c("a", "b"))
dds <- DESeq(dds)
resultsNames(dds)
rnms <- resultsNames(dds)
res <- results(dds)

# genes differentially expressed between batches
batches <- results(dds, contrast=c("reps", "a", "b"))

## biologically relevent comparisons
SWW.SD1 <- results(dds, contrast=c("ecocond", "SD1", "SWW1"))
SWW.SD1sh <- lfcShrink(dds,contrast=c("ecocond", "SD1", "SWW1"), res=SWW.SD1)
SD1.SRW <- results(dds, contrast=c("ecocond", "SWW2", "SD1"))
SD1.SRWsh <- lfcShrink(dds, contrast=c("ecocond", "SWW2", "SD1"), res=SD1.SRW)
SRW.SD2 <- results(dds, contrast=c("ecocond", "SD2", "SWW2"))
SRW.SD2sh <- lfcShrink(dds, contrast=c("ecocond", "SD2", "SWW2"), res=SRW.SD2)
SD2.SWW <- results(dds, contrast=c("ecocond", "SWW1", "SD2"))
SD2.SWWsh <- lfcShrink(dds, contrast=c("ecocond", "SWW1", "SD2"), res=SD2.SWW)
SRW.SWW <- results(dds, contrast=c("ecocond", "SWW1", "SWW2"))
SRW.SWWsh <- lfcShrink(dds, contrast=c("ecocond", "SWW1", "SWW2"), res=SRW.SWW)
YWW.YD1 <- results(dds, contrast=c("ecocond", "YD1", "YWW1"))
YWW.YD1sh <- lfcShrink(dds,contrast=c("ecocond", "YD1", "YWW1"), res=YWW.YD1)
YD1.YRW <- results(dds, contrast=c("ecocond", "YWW2", "YD1"))
YD1.YRWsh <- lfcShrink(dds, contrast=c("ecocond", "YWW2", "YD1"), res=YD1.YRW)
YRW.YD2 <- results(dds, contrast=c("ecocond", "YD2", "YWW2"))
YRW.YD2sh <- lfcShrink(dds, contrast=c("ecocond", "YD2", "YWW2"), res=YRW.YD2)
YD2.YWW <- results(dds, contrast=c("ecocond", "YWW1", "YD2"))
YD2.YWWsh <- lfcShrink(dds, contrast=c("ecocond", "YWW1", "YD2"), res=YD2.YWW)
YRW.YWW <- results(dds, contrast=c("ecocond", "YWW1", "YWW2"))
YRW.YWWsh <- lfcShrink(dds, contrast=c("ecocond", "YWW1", "YWW2"), res=YRW.YWW)


# other interesting comparisons
# we have SWW.SD1
SWW.SRW <- results(dds, contrast=c("ecocond", "SWW2", "SWW1"))
SWW.SRWsh <- lfcShrink(dds,contrast=c("ecocond", "SWW2", "SWW1"), res = SWW.SRW)
SWW.SD2 <- results(dds, contrast=c("ecocond", "SD2", "SWW1"))
SWW.SD2sh <- lfcShrink(dds,contrast=c("ecocond", "SD2", "SWW1"), res=SWW.SD2)
# we have YWW.YD1
YWW.YRW <- results(dds, contrast=c("ecocond", "YWW2", "YWW1"))
YWW.YRWsh <- lfcShrink(dds,contrast=c("ecocond", "YWW2", "YWW1"), res = YWW.YRW)
YWW.YD2 <- results(dds, contrast=c("ecocond", "YD2", "YWW1"))
YWW.YD2sh <- lfcShrink(dds,contrast=c("ecocond", "YD2", "YWW1"), res=YWW.YD2)

## genes for rt-qpcr comparison
genes <- c(grep("Thhalv10004906m.g", rownames(SWW.SD1)), grep("Thhalv10015427m.g", rownames(SWW.SD1)),
           grep("Thhalv10004736m.g",rownames(SWW.SD1)), grep("Thhalv10012591m.g", rownames(SWW.SD1)))


SD1 <- SWW.SD1[genes,c(2,3,6)]
SWW2 <- SWW.SRW[genes,c(2,3,6)]
SD2 <- SWW.SD2[genes,c(2,3,6)]
SWW1 <- data.frame(log2FoldChange = rep(0, 4), lfcSE = rep("0", 4), padj = rep("NA", 4))
rownames(SWW1) <- rownames(SD1)

YD1 <- YWW.YD1[genes,c(2,3,6)]
YWW2 <- YWW.YRW[genes,c(2,3,6)]
YD2 <- YWW.YD2[genes,c(2,3,6)]
YWW1 <- data.frame(log2FoldChange = rep(0, 4), lfcSE = rep("0", 4),padj = rep("NA", 4))
rownames(YWW1) <- rownames(YD1)

conditions <- c(rep("WW1", 4),
                rep("D1", 4), 
                rep("WW2", 4), 
                rep("D2", 4))

fcData <- data.frame(rbind(SWW1, SD1, SWW2, SD2,YWW1, YD1, YWW2, YD2), 
                     genes = as.factor(rep(rownames(SD1), 8)),
                     g_names = as.factor(rep(c("EsRAB18", "EsRD29A", "EsRD22", "EsERD1"),8)),
                     Ecotype = as.factor(c(rep("Shandong", 16), 
                                           rep("Yukon", 16))),
                     Conditions = rep(conditions, 2))



qpcr <- read.csv("./data_files/rtqpr_data.csv")
# melt is depreciated
#qpcr <- melt(qpcr, id=c("X", "Gene"))
qpcr <- qpcr %>%
  pivot_longer(cols=c(-X, -Gene), names_to = "variable", values_to = "value") # to stay consistent

se <- function(x) sd(x)/sqrt(length(x))

qpcrData <- qpcr %>% group_by(X, Gene, variable) %>% summarize(mean = mean(value), se=se(value)) %>% 
  data.frame() %>% mutate(g_names = c(rep(c("EsRD22","EsRD22", "EsRAB18", "EsRAB18", "EsERD1", "EsERD1", 
                          "EsRD29A", "EsRD29A"),4)), padj = rep("NA", 32))

qpcrData <- qpcrData[,c("mean", "se", "padj", "Gene", "g_names", "variable", "X")]
colnames(qpcrData) <- colnames(fcData)
qpcrData$Source <- rep("RT-qPCR", nrow(qpcrData))
fcData$Source <- rep("RNASeq", nrow(fcData))


fc_comb <- rbind(qpcrData, fcData)
fc_comb$Conditions <- factor(fc_comb$Conditions, c("WW1", "D1", "WW2", "D2"))
fc_comb$g_names <- factor(fc_comb$g_names, c("EsRAB18", "EsRD29A", "EsRD22", "EsERD1"))
fc_comb$Ecotype<- factor(fc_comb$Ecotype, levels = c("Shandong", "Yukon"))



fc_comb$sig <- ifelse(fc_comb$padj <= 0.05, "sig", NA)
fc_comb$lfcSE <- as.numeric(as.character(fc_comb$lfcSE))
(fc_combplot <- ggplot(fc_comb, aes(x = Conditions, y = log2FoldChange,colour = Ecotype,
                                   shape = Ecotype, group = Ecotype)) + 
    geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), width=.1) +
    geom_point(size = 4) +
    scale_shape_manual(values = c(16,15)) +
    geom_line(aes(linetype=Ecotype))  +
    geom_hline(yintercept=c(0), linetype="dotted") +
    scale_colour_grey(end=0.5) +
    coord_cartesian(ylim=c(-5,5)) +
    geom_text(data = subset(fc_comb, sig == "sig"), aes(x = Conditions, y = log2FoldChange+lfcSE), label = "*", size = 5, colour="red", 
               nudge_y = 0.1) +
    facet_grid(g_names ~ Source ) +
    theme_bw() +
    theme(strip.text = element_text(face = "italic"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    ylab("Log2 Fold Change") + xlab("Condition")) +
  theme(strip.text.y = element_text(face = "italic"), text=element_text(size=15))
ggsave("./figs/allfc_relevel.pdf", fc_combplot,
       height = 8, width = 7, units = "in")

#qpcrData$X <- factor(qpcrData$X, c("WW1", "D1", "WW2", "D2"))
#qpcrData$g_names <- factor(qpcrData$g_names, c("EsRAB18", "EsRD29A", "EsRD22", "EsERD1"))
#qpcrData$variable<- factor(qpcrData$variable, c("Shandong", "Yukon"))
#qpcrData$min <- qpcrData$mean-qpcrData$se
#qpcrData$sig <- ifelse(qpcrData$padj <= 0.05, "sig", NA)
#colnames(qpcrData)[3] <- "Ecotype"
#
#(qpcr_plot <- ggplot(qpcrData, aes(x = X, y = mean ,colour = Ecotype,
#                              shape = Ecotype, group = Ecotype)) + 
#  geom_point(size = 3) +
#  geom_line(aes(linetype=Ecotype))  +
#  geom_hline(yintercept=c(0), linetype="dotted") +
#  scale_colour_grey(end=0.5) +
#  coord_cartesian(ylim=c(-5,5)) +
#  facet_wrap( ~ g_names, ncol=2) +
#  theme_bw() +
#  theme(strip.text = element_text(face = "italic"),
#        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#  ylab("Log2 Fold Change") + xlab("Condition") + 
#  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1))

fcData$lfcSE <- as.character(fcData$lfcSE)
fcData$Conditions <- factor(fcData$Conditions, c("WW1", "D1", "WW2", "D2"))
fcData$g_names <- factor(fcData$g_names, c("EsRAB18", "EsRD29A", "EsRD22", "EsERD1"))
fcData$Ecotype<- factor(fcData$Ecotype, c("Yukon", "Shandong"))
fcData$min <- fcData$log2FoldChange-fcData$lfcSE
fcData$padj <- as.character(fcData$padj)
fcData$sig <- ifelse(fcData$padj <= 0.05, "sig", NA)
#(fc_plot <- ggplot(fcData, aes(x = Conditions, y = log2FoldChange,colour = Ecotype,
#                              shape = Ecotype, group = Ecotype)) + 
#  geom_point(size = 3) +
#  geom_line(aes(linetype=Ecotype))  +
#  geom_hline(yintercept=c(0), linetype="dotted") +
#  scale_colour_grey(end=0.5) +
#  coord_cartesian(ylim=c(-5,5)) +
#  facet_wrap( ~ g_names, ncol=2) +
#  theme_bw() +
#  theme(strip.text = element_text(face = "italic"),
#        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#  scale_y_continuous(name = "Log2 Fold Change") +
#  geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), width=.1))

fpkm <- read.csv("./data_files/2018-10-15-drought_geneFPKM_droughtpaper.txt", sep = '\t')
unique_fpkm <- fpkm[row.names(unique(fpkm[,c(1,10:52)])),]

rownames(unique_fpkm) <- unique_fpkm$gene_id
unique_fpkm <- unique_fpkm[,-1]

lncrna <- read.csv("./data_files/2018-09-17-drought_geneFPKM.txt", sep="\t")[,1:5]
## want to only count which DEGs with > 1fpkm in EACH of the samples from the paired conditions....
## also want to return the # of lncrnas for a particular ecotype...

## how many new genes? how many lncrna?
novel_dr <- unique_fpkm[which(substr(rownames(unique_fpkm),1,7) == "DROUGHT"),]
novel_exp <- novel_dr[apply(novel_dr[,c(9:15, 17:21, 23:25, 29:35, 37:42, 44:46)] > 1, 1, all),c(1:8,9:15, 17:21, 23:25, 29:35, 37:42, 44:46,52:54)]

which(!is.na(novel_exp$Ath.loci)) %>% length()
which(novel_exp$y_lnc == 1 | novel_exp$s_lnc ==1) %>% length()

total_expr <- unique_fpkm[apply(unique_fpkm[,c(9:15, 17:21, 23:25, 29:35, 37:42, 44:46)] > 1, 1, all),c(1:8,9:15, 17:21, 23:25, 29:35, 37:42, 44:46,52:54)]
which(unique_fpkm$y_lnc == 1 | unique_fpkm$s_lnc ==1) %>% length()
which(total_expr$y_lnc == 1 | total_expr$s_lnc ==1) %>% length()


## function to save DEGs in an organized list
getDEGS <- function(deg, condcolnames.1, condcolnames.2, eco){ 
  listname <- list()
  pos<- as.data.frame(deg[which(deg$padj < 0.05  
                                & deg$log2FoldChange > 0),])
  neg <- as.data.frame(deg[which(deg$padj < 0.05  
                                 & deg$log2FoldChange < 0),])
  #lncrna_eco <- lncrna[,eco]
  pos <- merge(pos, unique_fpkm, by=0)
  neg <- merge(neg, unique_fpkm,  by=0)
  listname[["pos"]] <- pos[apply(pos[,condcolnames.2] > 1, 1, all),] # 2nd condition
  listname[["neg"]] <- neg[apply(neg[,condcolnames.1] > 1, 1, all),] # 1st condition 
  lnc_pos <- listname[["pos"]] %>% filter(get(eco) == 1)
  lnc_neg <- listname[["neg"]] %>% filter(get(eco) == 1)
  listname[["pos_lnc"]] <- lnc_pos
  listname[["neg_lnc"]] <- lnc_neg
  perc_pos <- round((nrow(lnc_pos)/nrow(listname[["pos"]]))*100,2)
  perc_neg <- round((nrow(lnc_neg)/nrow(listname[["neg"]]))*100,2)
  total_deg <- (nrow(listname[["pos"]])+nrow(listname[["neg"]]))
  total_perc <- round(((nrow(lnc_pos)+nrow(lnc_neg))/(nrow(listname[["pos"]])+nrow(listname[["neg"]])))*100,2)                          
  print(paste("Number of upregulated genes:", nrow(listname[["pos"]])))
  print(paste("Upregulated AND lncRNA:", nrow(lnc_pos)))
  print(paste0(perc_pos, "%"))
  print(paste("Number of downregulated genes:", nrow(listname[["neg"]])))
  print(paste("Downregulated AND lncRNA:", nrow(lnc_neg)))
  print(paste0(perc_neg, "%"))
  print(paste("Total number of DEGs", total_deg))
  print(paste("Total % of lncRNA in DEGs:", paste0(total_perc, "%")))
  return(listname)
}

SWW1.SD1.degs <- getDEGS(SWW.SD1,  condcolnames.1 =  c("SWW.1", "SWW.2", "SWW1.3", 
                               "SWW1.4"), condcolnames.2=c("SD1.1", "SD1.2", "SD1.3"), eco=c("s_lnc"))
SD1.SWW2.degs <- getDEGS(SD1.SRW,condcolnames.1 = c("SD1.1", "SD1.2", "SD1.3"), 
                         condcolnames.2 = c("SRW.1", "SRW.2", "SWW2.3", "SWW2.4"), eco=c("s_lnc"))
SWW2.SD2.degs <- getDEGS(SRW.SD2, condcolnames.1 =  c("SRW.1", "SRW.2", "SWW2.3", "SWW2.4"),
                         condcolnames.2 =  c("SD2.1", "SD2.2", "SD2.3", "SD2.4", "SD2.5"), eco=c("s_lnc"))

YWW1.YD1.degs <- getDEGS(YWW.YD1,  condcolnames.1 =  c("YWW.1", "YWW.2", "YWW1.3", "YWW1.4"),  
                         condcolnames.2= c("YD1.1", "YD1.2", "YD1.3"), eco="y_lnc")
YD1.YWW2.degs <- getDEGS(YD1.YRW, condcolnames.1 = c("YD1.1", "YD1.2", "YD1.3"), 
                                  condcolnames.2 = c("YRW.1", "YRW.2", "YWW2.3", "YWW2.4"), eco="y_lnc")
YWW2.YD2.degs <- getDEGS(YRW.YD2,  condcolnames.1 =  c("YRW.1", "YRW.2", "YWW2.3", "YWW2.4"),
                         condcolnames.2 = c("YD2.1", "YD2.2", "YD2.3", "YD2.4"), eco="y_lnc")


SWW1SD1lnc <- rbind(SWW1.SD1.degs[['pos_lnc']], SWW1.SD1.degs[['neg_lnc']])
SD1SWW2lnc <- rbind(SD1.SWW2.degs[['pos_lnc']], SD1.SWW2.degs[['neg_lnc']])
SWW2SD2lnc <- rbind(SWW2.SD2.degs[['pos_lnc']], SWW2.SD2.degs[['neg_lnc']])
YWW1YD1lnc <- rbind(YWW1.YD1.degs[['pos_lnc']], YWW1.YD1.degs[['neg_lnc']])
YD1YWW2lnc <- rbind(YD1.YWW2.degs[['pos_lnc']], YD1.YWW2.degs[['neg_lnc']])
YWW2YD2lnc <- rbind(YWW2.YD2.degs[['pos_lnc']], YWW2.YD2.degs[['neg_lnc']])

## how many of SD2 lncrnas upreg are shandong only loci
## found in drought_wgcna.R (sourced in this file)
which(rownames(slnc_only) %in% SWW2.SD2.degs[['pos_lnc']]$Row.names) %>% length()

ylncdeg <- rbind(YWW1YD1lnc,YD1YWW2lnc,YWW2YD2lnc )
slncdeg <- rbind(SWW1SD1lnc,SD1SWW2lnc,SWW2SD2lnc)

ylncdeg <- data.frame(unique(ylncdeg$Row.names))
slncdeg <- data.frame(unique(slncdeg$Row.names))


slncdeg$shandong <- 1
ylncdeg$yukon <- 1

lncdegs <- merge(ylncdeg, slncdeg, by.x=1, by.y=1 , all = TRUE)
lncdegs$yukon[is.na(lncdegs$yukon)] <- 0
lncdegs$shandong[is.na(lncdegs$shandong)] <- 0

lncdegs <- lncdegs %>% dplyr::mutate(both = ifelse(yukon == 1 & shandong == 1, 1, 0))


#write.csv(lncdegs, "drought_deg_lncs.csv")

getDEGoverlap <- function(deg1, deg2){
  pos1 <- deg1[["pos"]]
  pos2 <- deg2[["pos"]]
  neg1 <- deg1[["neg"]]
  neg2 <- deg2[["neg"]]
  pos_lnc1 <- deg1[["pos_lnc"]]
  pos_lnc2 <- deg2[["pos_lnc"]]
  neg_lnc1 <- deg1[["neg_lnc"]]
  neg_lnc2 <- deg2[["neg_lnc"]]
  pos_int <- intersect(toupper(pos1$Row.names), toupper(pos2$Row.names)) 
  neg_int <- intersect(toupper(neg1$Row.names), toupper(neg2$Row.names)) 
  pos_lnc_int <- intersect(toupper(pos_lnc1$Row.names), toupper(pos_lnc2$Row.names)) 
  neg_lnc_int <- intersect(toupper(neg_lnc1$Row.names), toupper(neg_lnc2$Row.names))
  print(paste("upregulated intersect", length(pos_int)))
  print(paste("downregulated intersect", length(neg_int)))
  print(paste("up lnc intersect", length(pos_lnc_int)))
  print(paste("down lnc intersect", length(neg_lnc_int)))
  }

getOverlap_genes <- function(deg1, deg2){
  return_list <- list()
  pos1 <- deg1[["pos"]]
  pos2 <- deg2[["pos"]]
  neg1 <- deg1[["neg"]]
  neg2 <- deg2[["neg"]]
  pos_lnc1 <- deg1[["pos_lnc"]]
  pos_lnc2 <- deg2[["pos_lnc"]]
  neg_lnc1 <- deg1[["neg_lnc"]]
  neg_lnc2 <- deg2[["neg_lnc"]]
  #print(toupper(rownames(pos1)))
  pos_int <- intersect(toupper(pos1$Row.names), toupper(pos2$Row.names)) 
  neg_int <- intersect(toupper(neg1$Row.names), toupper(neg2$Row.names)) 
  pos_lnc_int <- intersect(toupper(pos_lnc1$Row.names), toupper(pos_lnc2$Row.names)) 
  neg_lnc_int <- intersect(toupper(neg_lnc1$Row.names), toupper(neg_lnc2$Row.names))
  print(head(pos_lnc_int))
  x <- unique_fpkm
  rownames(x) <- rownames(x) %>% toupper(.)
  #pos <- merge(pos_int, x, by=0)
  #neg <- merge(neg_int, x,  by=0)
  #print(length(pos_int))
  pos <- x[rownames(x) %in% pos_int,]
  neg <- x[rownames(x) %in% neg_int,]
  pos_lnc <- x[rownames(x) %in% pos_lnc_int,]
  neg_lnc <- x[rownames(x) %in% neg_lnc_int,] 
  #print(dim(pos))
  return_list[["pos"]] <- pos
  return_list[["neg"]] <- neg
  return_list[["pos_lnc"]] <- pos_lnc
  return_list[["neg_lnc"]] <- neg_lnc
  print(paste("upregulated intersect", length(pos_int)))
  print(paste("downregulated intersect", length(neg_int)))
  print(paste("up lnc intersect", length(pos_lnc_int)))
  print(paste("down lnc intersect", length(neg_lnc_int)))
  return(return_list)
  }


YWW1.YD1.degs[["pos"]]$Row.names %>% intersect(., YD1.YWW2.degs[['neg']]$Row.names) %>% length()
YWW1.YD1.degs[["neg"]]$Row.names %>% intersect(., YD1.YWW2.degs[['pos']]$Row.names) %>% length()


getDEGoverlap(SWW1.SD1.degs, YWW1.YD1.degs)
getDEGoverlap(SD1.SWW2.degs, YD1.YWW2.degs)
getDEGoverlap(SWW2.SD2.degs, YWW2.YD2.degs)
getDEGoverlap(SWW1.SD1.degs, SWW2.SD2.degs)
getDEGoverlap(YWW1.YD1.degs, YWW2.YD2.degs)


WW1_D1 <- getOverlap_genes(SWW1.SD1.degs, YWW1.YD1.degs)
D1_WW2 <- getOverlap_genes(SD1.SWW2.degs, YD1.YWW2.degs)
WW2_D2 <- getOverlap_genes(SWW2.SD2.degs, YWW2.YD2.degs)
YD1.YD2 <- getOverlap_genes(YWW1.YD1.degs, YWW2.YD2.degs)
YD1.SD2.overlap <- getOverlap_genes(YWW1.YD1.degs, SWW2.SD2.degs)

## 
deg_lncrnas_elizabeth <- c(rownames(D1_WW2[['pos_lnc']]), rownames(D1_WW2[['neg_lnc']]), 
  rownames(WW2_D2[['pos_lnc']]), rownames(WW2_D2[['neg_lnc']])) %>% unique()
deg_lncrnas_dup <- c(rownames(D1_WW2[['pos_lnc']]), rownames(D1_WW2[['neg_lnc']]), 
                           rownames(WW2_D2[['pos_lnc']]), rownames(WW2_D2[['neg_lnc']])) %>% duplicated()

#write.csv(fpkm_liz, "../tables/droughtlncrnas.csv", row.names = F, quote = F, col.names = T)

## YD1.YD2 lnc....
yd1.yd2.pos.lnc <- rownames(YD1.YD2[['pos_lnc']])
yd1.yd2.neg.lnc <- rownames(YD1.YD2[['neg_lnc']])



## downregulated in YWW2
seven_lncs <- YD1.YWW2.degs[['neg_lnc']][YD1.YWW2.degs[['neg_lnc']]$Row.names %in% yd1.yd2.pos.lnc, "Row.names"] 
## 
YD1.YWW2.degs[['pos_lnc']][YD1.YWW2.degs[['pos_lnc']]$Row.names %in% yd1.yd2.neg.lnc,]

## up in both YD1 and YD2 AND in SD2
SWW2.SD2.degs[['pos_lnc']][SWW2.SD2.degs[['pos_lnc']]$Row.names %in% yd1.yd2.pos.lnc,]

## up in BOTH YD1 and YD2, down in YWW2, up in SD2
SWW2.SD2.degs[['pos_lnc']][SWW2.SD2.degs[['pos_lnc']]$Row.names %in% seven_lncs,]


## DEG patterns of interest.... 
getPattern_genes <- function(overlap1, overlap2, overlap3, dir1, dir2, dir3){
  return_list <- list()
  genes1 <- rownames(overlap1[[dir1]])
  genes2 <- rownames(overlap2[[dir2]])
  genes3 <- rownames(overlap3[[dir3]])
  overlap12 <- intersect(genes1, genes2)
  overlap123 <- intersect(overlap12, genes3)
  x <- unique_fpkm
  rownames(x) <- rownames(x) %>% toupper(.)
  overlap_genes <- x[rownames(x) %in% overlap123,]
  return_list[['wide']] <- overlap_genes
  overlap_genes <- overlap_genes[,c(9:15, 17:21, 23:25, 29:35, 37:42, 44:46)]
  colnames(overlap_genes)[c(1,2,8,9,16,17,23,24)] <- c("YWW1.1", "YWW1.2", "YWW2.1", "YWW2.2",
                                                    "SWW1.1", "SWW1.2", "SWW2.1", "SWW2.2")
  overlap_genes$genes <- rownames(overlap_genes)
  melted <- melt(overlap_genes, id.vars="genes")
  melted$ecocond <- c(rep("YWW", 4*nrow(overlap_genes)), 
                        rep("YD1", 3*nrow(overlap_genes)),
                        rep("YRW", 4*nrow(overlap_genes)),
                        rep("YD2", 4*nrow(overlap_genes)),
                        rep("SWW", 4*nrow(overlap_genes)), 
                        rep("SD1", 3*nrow(overlap_genes)),
                         rep("SRW", 4*nrow(overlap_genes)),
                         rep("SD2", 5*nrow(overlap_genes)))
  melted$condition <- as.factor(substr(melted$ecocond, 2, nchar(melted$ecocond)))
  melted$condition <- factor(melted$condition, levels=c("WW", "D1", "RW", "D2")) 
  melted$ecotype <- ifelse(substr(melted$ecocond, 1,1) == "Y", "Yukon", "Shandong")
  return_list[['long']] <- melted
  return(return_list)
}

getPattern_genes.degs <- function(degs1, overlap2, overlap3, dir1, dir2, dir3){
  return_list <- list()
  genes1 <- toupper(degs1[[dir1]]$Row.names)
  genes2 <- rownames(overlap2[[dir2]])
  genes3 <- rownames(overlap3[[dir3]])
  overlap12 <- intersect(genes1, genes2)
  overlap123 <- intersect(overlap12, genes3)
  x <- unique_fpkm
  rownames(x) <- rownames(x) %>% toupper(.)
  overlap_genes <- x[rownames(x) %in% overlap123,]
  return_list[['wide']] <- overlap_genes
  overlap_genes <- overlap_genes[,c(9:15, 17:21, 23:25, 29:35, 37:42, 44:46)]
  colnames(overlap_genes)[c(1,2,8,9,16,17,23,24)] <- c("YWW1.1", "YWW1.2", "YWW2.1", "YWW2.2",
                                                       "SWW1.1", "SWW1.2", "SWW2.1", "SWW2.2")
  overlap_genes$genes <- rownames(overlap_genes)
  melted <- melt(overlap_genes, id.vars="genes")
  melted$ecocond <- c(rep("YWW", 4*nrow(overlap_genes)), 
                      rep("YD1", 3*nrow(overlap_genes)),
                      rep("YRW", 4*nrow(overlap_genes)),
                      rep("YD2", 4*nrow(overlap_genes)),
                      rep("SWW", 4*nrow(overlap_genes)), 
                      rep("SD1", 3*nrow(overlap_genes)),
                      rep("SRW", 4*nrow(overlap_genes)),
                      rep("SD2", 5*nrow(overlap_genes)))
  melted$condition <- as.factor(substr(melted$ecocond, 2, nchar(melted$ecocond)))
  melted$condition <- factor(melted$condition, levels=c("WW", "D1", "RW", "D2")) 
  melted$ecotype <- ifelse(substr(melted$ecocond, 1,1) == "Y", "Yukon", "Shandong")
  return_list[['long']] <- melted
  return(return_list)
}

# Up D1 D WW2 D D2
pattern1 <- getPattern_genes(overlap1=WW1_D1, overlap2=D1_WW2, overlap3=WW2_D2, 
                             dir1="pos", dir2="neg", dir3="neg") # 0 genes
# U D1 U WW2 U D2 
pattern2 <- getPattern_genes(overlap1=WW1_D1, overlap2=D1_WW2, overlap3=WW2_D2, 
                             dir1="pos", dir2="pos", dir3="pos") # 0

# Up D1, down WW2, up D2
pattern3 <- getPattern_genes(overlap1=WW1_D1, overlap2=D1_WW2, overlap3=WW2_D2, 
             dir1="pos", dir2="neg", dir3="pos") # 8 genes
# Down D1, up WW2, down D2
pattern4 <- getPattern_genes(overlap1=WW1_D1, overlap2=D1_WW2, overlap3=WW2_D2, 
             dir1="neg", dir2="pos", dir3="neg") #4 genes

# U D1 U WW2 D D2
pattern5 <- getPattern_genes(overlap1=WW1_D1, overlap2=D1_WW2, overlap3=WW2_D2, 
                             dir1="pos", dir2="pos", dir3="neg") #0 genes

# D D1 DWW2 U D2
pattern6 <- getPattern_genes(overlap1=WW1_D1, overlap2=D1_WW2, overlap3=WW2_D2, 
                             dir1="neg", dir2="neg", dir3="pos") #0 genes
# d d1 dww2 dd2
pattern7 <- getPattern_genes(overlap1=WW1_D1, overlap2=D1_WW2, overlap3=WW2_D2, 
                             dir1="neg", dir2="neg", dir3="neg") #0 genes

# d d1 u ww2 u d2
pattern8 <- getPattern_genes(overlap1=WW1_D1, overlap2=D1_WW2, overlap3=WW2_D2, 
                             dir1="neg", dir2="pos", dir3="pos") #0 genes


## pattern 9

pattern9 <- getPattern_genes.degs(degs1=YWW1.YD1.degs, overlap2=D1_WW2, overlap3=WW2_D2,
                             dir1="pos", dir2="neg", dir3="pos")

## pattern 10
pattern10 <- getPattern_genes.degs(degs1=YWW1.YD1.degs, overlap2=D1_WW2, overlap3=WW2_D2,
                                   dir1="neg", dir2="pos", dir3="neg")

data3 <- data.frame(pattern3[['long']])
data4 <- data.frame(pattern4[['long']])
#library(plyr)
data3$genes <- revalue(data3$genes, c("THHALV10023585M.G"="EsGolS4", "THHALV10024122M.G"="EsGolS2",
                                    "THHALV10024601M.G"="EsARF16"))


pattern3genes <- row.names(data.frame(pattern3[['wide']]))
pattern4genes <- row.names(data.frame(pattern4[["wide"]]))

data3$direction <- rep("Upregulated by drought", nrow(data3))
data4$direction <- rep("Downregulated by drought", nrow(data4))
bothdir <- rbind(data3,data4)
    bothdir$condition <- revalue(bothdir$condition, c("WW"="WW1", "RW"="WW2"))
bothdir$direction <- factor(bothdir$direction, 
                            levels=c("Upregulated by drought", "Downregulated by drought" ))
## plotting both patterns together
matchColours <- colorRampPalette(c("#E7B800","#67A7C1", "#FF6F59", "#292F36"))(13)
## wanted to turn matching colours into more pastel-like...but...it looks better as original colours
# turning to RGB
pastelly <- readhex(file = textConnection(paste(matchColours, collapse = "\n")),
                 class = "RGB")
#transform to hue/lightness/saturation colorspace
pastelly <- as(pastelly, "HLS")
#additive decrease of lightness
pastelly@coords[, "L"] <- pmax(0, pastelly@coords[, "L"] + 0.2)
#going via rgb seems to work better  
pastelly <- as(pastelly, "RGB") %>% hex()

(pat34plot <- ggplot(bothdir, aes(x=condition, y=value,color=genes, group=genes)) + 
    stat_summary(fun.data = 'mean_se', geom = "errorbar",colour="darkgrey", width=.2) +
    stat_summary(aes(x = condition, y = value),  color="black",
                 fun.y = 'mean', geom = 'line', 
                 position = 'dodge') +
    stat_summary(aes(x = condition, y = value, fill=genes),
                 fun.y = 'mean', geom = 'point', 
                 position = 'dodge',colour="black",pch=21, size=5) + 
    facet_grid(direction~ecotype) + xlab("Condition") +
    #scale_fill_brewer(palette = "Set3") + 
    scale_fill_manual(values=matchColours, name="Gene name") +
    ylab("Mean expression estimate (FPKM)") + theme_bw() +
    theme(text=element_text(size=18), legend.title=element_blank()) +
    theme(legend.position="bottom", legend.text=element_text(size=10))+
    guides(fill = guide_legend(nrow = 3)))

ggsave("./figs/pattern3and4update.pdf", plot=pat34plot, width=11, height=7, units="in")


possigDEG <- function(resultofdeg){
  x<- as.data.frame(resultofdeg[which(resultofdeg$padj < 0.05  
                                      & resultofdeg$log2FoldChange > 0) 
                                #& substr(rownames(resultofdeg), 1,6) == "Thhalv")
                                ,])
  return(x)
}

negsigDEG <- function(resultofdeg){
  x<- as.data.frame(resultofdeg[which(resultofdeg$padj < 0.05  
                                      & resultofdeg$log2FoldChange < 0)
                                #& substr(rownames(resultofdeg), 1,6) == "Thhalv")
                                ,])
  return(x)
}



allDEG <- function(resultofdeg){
  x<- as.data.frame(resultofdeg[which(resultofdeg$padj < 0.05 ),])
  return(rownames(x))
}

sww1.sd1d <- allDEG(SWW.SD1)
sd1.sww2d <- allDEG(SD1.SRW)
sww2.sd2d <- allDEG(SRW.SD2)
sd2.sww1d <- allDEG(SD2.SWW)
sww2.sww1d <- allDEG(SRW.SWW)
yww1.yd1d <- allDEG(YWW.YD1)
yd1.yww2d <- allDEG(YD1.YRW)
yww2.yd2d <- allDEG(YRW.YD2)
yd2.yww1d <- allDEG(YD2.YWW)
yww2.yww1d <- allDEG(YRW.YWW)

alldegs <- c(sww1.sd1d, sd1.sww2d,sww2.sd2d, yww1.yd1d, yd1.yww2d, yww2.yd2d) %>% unique()
shandeg <-  c(sww1.sd1d, sd1.sww2d,sww2.sd2d) %>% unique()
yukdeg <- c(yww1.yd1d, yd1.yww2d, yww2.yd2d)%>% unique()


intersect(sd1.sww2d, yd1.yww2d) %>% length()
c(sd1.sww2d, yd1.yww2d) %>% unique() %>% length()
c(sww1.sd1d, yww1.yd1d) %>% unique() %>% length()
c(yww2.yd2d, sww2.sd2d) %>% unique() %>% length()





