library(ggplot2)
library(ggrepel)
library(viridis)
library(DESeq2)
library(readr)
library(dplyr)



## same way to import FPKM values
fpkm <- read.csv("./data_files/2018-10-15-drought_geneFPKM_droughtpaper.txt", sep="\t") 

unique_fpkm <- fpkm[row.names(unique(fpkm[,c(1,10:52)])),]

rownames(unique_fpkm) <- unique_fpkm$gene_id

unique_fpkm <- unique_fpkm[,-1]
y.lnc <- rownames(unique_fpkm[which(unique_fpkm[,"y_lnc", drop=F] == 1),])
s.lnc <- rownames(unique_fpkm[which(unique_fpkm[,"s_lnc", drop=F] == 1),])
rep_lib <- unique_fpkm[,c(9:15, 17:22, 23:25, 29:35, 37:43, 44:46)]
fpkmunique <- unique_fpkm[,c(9:15, 17:21, 23:25, 29:35, 37:42, 44:46)]
colnames(fpkmunique)[c(1,2,8,9,16,17,23,24)] <- c("YWW1.1", "YWW1.2", "YWW2.1", "YWW2.2",
                                                  "SWW1.1", "SWW1.2", "SWW2.1", "SWW2.2")
fpkmlog <- log2(fpkmunique+1)


####################################################
## looking for replication in re-sequenced samples #
####################################################

colnames(rep_lib)[c(1,2,8,9,17,18,24,25)] <- c("YWW1.1", "YWW1.2", "YWW2.1", "YWW2.2",
                                                  "SWW1.1", "SWW1.2", "SWW2.1", "SWW2.2")
pcarep <- prcomp(log2(rep_lib+1), scale=F)
exprVals<-data.frame(pcarep$x)
sampleVals<-data.frame(pcarep$rotation)


dim(exprVals) 
dim(sampleVals)

samples <- substring(colnames(rep_lib),1,nchar(colnames(rep_lib))-2)
ecotype <- c(rep("Yukon", 16), rep("Shandong", 17))


## extrating all PCA data for ggplot
coords<-data.frame(sampleVals, 
                   Samples = factor(samples), Name = colnames(rep_lib),
                   Condition = factor(substr(samples, 2,nchar(samples))), 
                   Ecotype = ecotype)
coords$Condition[c(13,30)] <- c("D2","D2")
coords$Condition <- factor(coords$Condition, c("WW1","D1", "WW2", "D2"))


PoV <- (pcarep$sdev^2/sum(pcarep$sdev^2)) * 100
PC1per <- paste0("(", round(PoV[1], 2), "%", ")")
PC2per <- paste0("(", round(PoV[2], 2), "%", ")")
PC3per <- paste0("(", round(PoV[3],2) , "%", ")")
PC6per <- paste0("(", round(PoV[6],2) , "%", ")")
PC4per <- paste0("(", round(PoV[4],2) , "%", ")")
PC5per <- paste0("(", round(PoV[5],2) , "%", ")")

                    
(pc12plot_rep<-ggplot(exprVals, aes(PC2, PC1)) + 
    geom_point(data = coords, aes(x=PC2, y=PC1, fill = Condition, shape = Ecotype), 
               size=5) +  #, show.legend=F) +
    coord_cartesian(xlim=c(-0.55,0.55), ylim=c(-0.55,0.55)) +
    scale_x_continuous(name = paste("PC2", PC2per)) +
    scale_y_continuous(name = paste("PC1", PC1per)) +
    geom_vline(xintercept=c(0), linetype="dotted") +
    geom_hline(yintercept=c(0), linetype="dotted") +
    #geom_text_repel(data=coords,aes(PC2,PC1,label=Name), force=1, point.padding = 0.5,size=4) + 
    #scale_fill_brewer(palette = "Set3") +
    scale_fill_manual(values = c("#E7B800","#67A7C1", "#FF6F59", "#292F36")) +
    theme_classic() +
    scale_shape_manual(values=c(21,22)) +
    theme(text=element_text(size=18)) +
    guides(fill = guide_legend(override.aes=list(shape=21))))


(pc24plot_rep<-ggplot(exprVals, aes(PC2, PC4)) + 
    geom_point(data = coords, aes(x=PC2, y=PC4, fill = Condition, shape = Ecotype), 
               size=5) +  #, show.legend=F) +
    coord_cartesian(xlim=c(-0.55,0.55), ylim=c(-0.55,0.55)) +
    scale_x_continuous(name = paste("PC2", PC2per)) +
    scale_y_continuous(name = paste("PC4", PC4per)) +
    geom_vline(xintercept=c(0), linetype="dotted") +
    geom_hline(yintercept=c(0), linetype="dotted") +
    geom_text_repel(data=coords,aes(PC2,PC4,label=Name), force=1, point.padding = 0.7,size=4) + 
    #scale_fill_brewer(palette = "Set3") +
    scale_fill_manual(values = c("#E7B800","#67A7C1", "#FF6F59", "#292F36")) +
    theme_classic() +
    scale_shape_manual(values=c(21,22)) +
    theme(text=element_text(size=18)) +
    guides(fill = guide_legend(override.aes=list(shape=21))))
ggsave("./figs/pc24plot_rep.pdf", plot=pc24plot_rep, 
       width=7, height=5, units="in")

##################################################
# Looking for influence of library prep protocol #
##################################################
pca <- prcomp(fpkmlog, center=T, scale=FALSE)
plot(pca, type = "l")
x.var <- pca$sdev ^ 2
x.pvar <- x.var/sum(x.var)

plot(x.pvar,xlab="Principal component", ylab="Proportion of variance explained", ylim=c(0,1), type='b')

summary(pca)

exprVals<-data.frame(pca$x)
sampleVals<-data.frame(pca$rotation)


dim(exprVals) 
dim(sampleVals)

samples <- substring(colnames(fpkmlog),1,nchar(colnames(fpkmlog))-2)
ecotype <- c(rep("Yukon", 15), rep("Shandong", 16))


## extrating all PCA data for ggplot
coords<-data.frame(sampleVals, 
                   Samples = factor(samples), Name = colnames(fpkmlog),
                   Condition = factor(substr(samples, 2,nchar(samples))), 
                   Ecotype = ecotype)

coords$Condition <- factor(coords$Condition, c("WW1","D1", "WW2", "D2"))


PoV <- (pca$sdev^2/sum(pca$sdev^2)) * 100
PC1per <- paste0("(", round(PoV[1], 2), "%", ")")
PC2per <- paste0("(", round(PoV[2], 2), "%", ")")
PC3per <- paste0("(", round(PoV[3],2) , "%", ")")
PC6per <- paste0("(", round(PoV[6],2) , "%", ")")
PC4per <- paste0("(", round(PoV[4],2) , "%", ")")
PC5per <- paste0("(", round(PoV[5],2) , "%", ")")

coords$EcotypeLib <- c("Yukon-A", "Yukon-A", "Yukon-B", "Yukon-B", "Yukon-A", "Yukon-A", "Yukon-B",
                           "Yukon-A", "Yukon-A", "Yukon-B", "Yukon-B", "Yukon-A", "Yukon-A", "Yukon-B", "Yukon-B",
                           "Shandong-A", "Shandong-A", "Shandong-B", "Shandong-B", "Shandong-A", "Shandong-A", "Shandong-B",
                           "Shandong-A", "Shandong-A", "Shandong-B", "Shandong-B", "Shandong-A", "Shandong-A", "Shandong-B", "Shandong-B", "Shandong-B")
(pc12plot<-ggplot(exprVals, aes(PC1, PC2)) + 
    geom_point(data = coords, aes(x=PC1, y=PC2, fill = Condition, shape = EcotypeLib), 
               size=5) +  #, show.legend=F) +
    #coord_cartesian(xlim=c(-0.55,0.55), ylim=c(-0.55,0.55)) +
    scale_x_continuous(name = paste("PC1", PC1per)) +
    scale_y_continuous(name = paste("PC2", PC2per)) +
    geom_vline(xintercept=c(0), linetype="dotted") +
    geom_hline(yintercept=c(0), linetype="dotted") +
    geom_text_repel(data=coords,aes(PC1,PC2,label=Name), force=1, point.padding = 0.5,size=4) + 
    #scale_fill_brewer(palette = "Set3") +
    scale_fill_manual(values = c("#E7B800","#67A7C1", "#FF6F59", "#292F36")) +
    theme_classic() +
    scale_shape_manual(values=c(21,24, 22,25)) +
    theme(text=element_text(size=18)) +
    guides(fill = guide_legend(override.aes=list(shape=21))))

(pc13plot_new<-ggplot(exprVals, aes(PC1, PC3)) + 
    geom_point(data = coords, aes(x=PC1, y=PC3, fill = Condition, shape = Ecotype), 
               size=5) +  #, show.legend=F) +
    # coord_cartesian(xlim=c(-0.55,0.55), ylim=c(-0.55,0.55)) +
    scale_x_continuous(name = paste("PC1", PC1per)) +
    scale_y_continuous(name = paste("PC3", PC3per)) +
    geom_vline(xintercept=c(0), linetype="dotted") +
    geom_hline(yintercept=c(0), linetype="dotted") +
    geom_text_repel(data=coords,aes(PC1,PC3,label=Name), force=1, point.padding = 0.5,size=4) + 
    #scale_fill_brewer(palette = "Set3") +
    scale_fill_manual(values = c("#E7B800","#67A7C1", "#FF6F59", "#292F36")) +
    theme_classic() +
    scale_shape_manual(values=c(21,22)) +
    theme(text=element_text(size=18),
          panel.border = element_rect(colour = "black", fill=NA, size=2)) +
    guides(fill = guide_legend(override.aes=list(shape=21))))

(pc23plot<-ggplot(exprVals, aes(PC2, PC3)) + 
  geom_point(data = coords, aes(x=PC2, y=PC3, fill = Condition, shape = EcotypeLib), 
             size=5) +  #, show.legend=F) +
  #coord_cartesian(xlim=c(-0.55,0.55), ylim=c(-0.55,0.55)) +
  scale_x_continuous(name = paste("PC2", PC2per)) +
  scale_y_continuous(name = paste("PC3", PC3per)) +
  geom_vline(xintercept=c(0), linetype="dotted") +
  geom_hline(yintercept=c(0), linetype="dotted") +
  geom_text_repel(data=coords,aes(PC2,PC3,label=Name), force=1, point.padding = 0.5,size=4) + 
  #scale_fill_brewer(palette = "Set3") +
  scale_fill_manual(values = c("#E7B800","#67A7C1", "#FF6F59", "#292F36")) +
  theme_classic() +
  scale_shape_manual(values=c(21,24, 22,25)) +
  theme(text=element_text(size=18)) +
    guides(fill = guide_legend(override.aes=list(shape=21)))) #+ annotate("text", x = .4, y = .5, label = "All other transcripts", size=8)
ggsave("./figs/pc23plot.pdf", plot=pc23plot, width=7, height=5, units="in")

(pc24plot<-ggplot(exprVals, aes(PC2, PC4)) + 
    geom_point(data = coords, aes(x=PC2, y=PC4, fill = Condition, shape = Ecotype), 
               size=5) +  #, show.legend=F) +
   # coord_cartesian(xlim=c(-0.55,0.55), ylim=c(-0.55,0.55)) +
    scale_x_continuous(name = paste("PC2", PC2per)) +
    scale_y_continuous(name = paste("PC4", PC4per)) +
    geom_vline(xintercept=c(0), linetype="dotted") +
    geom_hline(yintercept=c(0), linetype="dotted") +
    geom_text_repel(data=coords,aes(PC2,PC4,label=Name), force=1, point.padding = 0.5,size=4) + 
    #scale_fill_brewer(palette = "Set3") +
    scale_fill_manual(values = c("#E7B800","#67A7C1", "#FF6F59", "#292F36")) +
    theme_classic() +
    scale_shape_manual(values=c(21,22)) +
    theme(text=element_text(size=18),
          panel.border = element_rect(colour = "black", fill=NA, size=2)) +
    guides(fill = guide_legend(override.aes=list(shape=21))))
ggsave("./figs/pc24plot.pdf", plot=pc24plot, width=9, height=7, units="in")


 PC34plot<-ggplot(exprVals, aes(PC3, PC4)) + 
  #geom_point(data = coords, aes(x=PC3, y=PC4,color = Condition, shape = Ecotype), size=5) +
  geom_point(data=coords, aes(x=PC3, y=PC4, fill = "black", shape = Ecotype), size = 6, show.legend=F) +
  geom_point(data = coords, aes(x=PC3, y=PC4, color = Condition, shape = Ecotype), 
             size=5) + #, show.legend=F) +
  coord_cartesian(xlim=c(-0.55,0.55), ylim=c(-0.55,0.55)) +
  scale_x_continuous(name = paste("PC3", PC3per)) +
  scale_y_continuous(name = paste("PC4", PC4per)) +
  geom_vline(xintercept=c(0), linetype="dotted") +
  geom_hline(yintercept=c(0), linetype="dotted") +
  geom_text_repel(data=coords,aes(PC3,PC4,label=Name), force=1, point.padding = 0.5, size=4) + 
  scale_color_brewer(palette = "Set3") + 
  #viridis::scale_color_viridis(discrete = TRUE, option = "C") +
  #scale_y_continuous(minor_breaks = seq(-0.5 , 0.5, 0.25), breaks = seq(-0.5, 0.5, 0.5)) +
  theme_classic() + theme(text=element_text(size=15)) 
 ggsave("./figs/pc34plot.pdf", plot=PC34plot, width=9, height=7, units="in")
 





