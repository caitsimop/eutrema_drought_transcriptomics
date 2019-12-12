library(tidyverse)

## how many drought lncrnas were identified by other method?
fpkm <- read.csv("./data_files/2018-10-15-drought_geneFPKM_droughtpaper.txt", sep="\t") # 
unique_fpkm <- fpkm[row.names(unique(fpkm[,c(1,10:52)])),]

rownames(unique_fpkm) <- unique_fpkm$gene_id

unique_fpkm <- unique_fpkm[,-1]

lncrna_classcodes <- read.csv("./data_files/cmp_to_signal.droughtlncrna.gtf.tmap", sep='\t')

fpkm_codes <- merge(fpkm[,1:51], lncrna_classcodes, by.x='X', by.y = 'qry_id', all.y=T)
fpkm_codes <- fpkm_codes[!duplicated(fpkm_codes$gene_id), ]

fpkm_codes <- fpkm_codes[,-c(17,23,37,44, 27:29, 48:51)]
rownames(fpkm_codes) <- fpkm_codes$gene_id
code_expr <- fpkm_codes[,10:40]
code_expr[code_expr < 1] <- 0
code_expr <- code_expr[ rowSums(code_expr)!=0, ] 
code_expr <- merge(code_expr, fpkm_codes[,c('gene_id', 'class_code')], by.x='row.names', by.y='gene_id', all.x=T)




# Gene names repeated for lncRNA prediction only
unique_fpkm <- fpkm[row.names(unique(fpkm[,c(1,10:52)])),]
rownames(unique_fpkm) <- unique_fpkm$gene_id
unique_fpkm <- unique_fpkm[,-1]

## predicted lncRNAs in each ecotype
y.lnc <- rownames(unique_fpkm[which(unique_fpkm[,"y_lnc", drop=F] == 1),])
s.lnc <- rownames(unique_fpkm[which(unique_fpkm[,"s_lnc", drop=F] == 1),])

# same lncRNA in each ecotype
samepredictlncrna <- intersect(y.lnc, s.lnc)
fpkmunique <- unique_fpkm[,c(9:15, 17:21, 23:25, 29:35, 37:42, 44:46)]
colnames(fpkmunique)[c(1,2,8,9,16,17,23,24)] <- c("YWW1.1", "YWW1.2", "YWW2.1", "YWW2.2",
                                                  "SWW1.1", "SWW1.2", "SWW2.1", "SWW2.2")

y.lnc.exp <- fpkmunique[rownames(fpkmunique) %in% y.lnc,] 
s.lnc.exp <- fpkmunique[rownames(fpkmunique) %in% s.lnc,] 


## lncrnas ecotype specific
y.lnc.exp[apply(y.lnc.exp[,1:15] > 1, 1, any) &  
            apply(y.lnc.exp[,16:31] < 1, 1, all),] %>% dim()
ylnc_only <- y.lnc.exp[apply(y.lnc.exp[,1:15] > 1, 1, any) &  
                         apply(y.lnc.exp[,16:31] < 1, 1, all),]

s.lnc.exp[apply(s.lnc.exp[,16:31] > 1, 1, any) & 
            apply(s.lnc.exp[,1:15] < 1, 1, all),] %>% dim()
slnc_only <- s.lnc.exp[apply(s.lnc.exp[,16:31] > 1, 1, any) & 
                         apply(s.lnc.exp[,1:15] < 1, 1, all),]
## total unique lncrnas 
## these are lncRNAs that are expressed at any time in Y and S
ylnc_1 <- y.lnc.exp[apply(y.lnc.exp[,1:15] > 1, 1, any),] 
slnc_1 <- s.lnc.exp[apply(s.lnc.exp[,16:31] > 1, 1, any),]
unique_lncrna <- c(rownames(ylnc_1), rownames(slnc_1)) %>% base::unique()
both_lncrna <- intersect(rownames(ylnc_1),rownames(slnc_1))

## lncRNAs that are expressed only in stress treated libraries 

unique_lncrna_expr <- fpkmunique[rownames(fpkmunique) %in% unique_lncrna,]
x <- unique_lncrna_expr[(apply(unique_lncrna_expr[,c(5:7)] > 1, 1, any) | 
                           apply(unique_lncrna_expr[,c(13:15)] > 1, 1, any) |
                           apply(unique_lncrna_expr[,c(20:22)] > 1, 1, any) |
                           apply(unique_lncrna_expr[,c(27:31)] > 1, 1, any)) &
                          apply(unique_lncrna_expr[,-c(5:7,13:15,20:22,27:31)] < 1, 1, all),]


## do lncRNAs have more isoforms on average than protein coding genes?
lncrna_iso <- fpkm[fpkm$gene_id %in% unique_lncrna,] %>% dplyr::select(gene_id, X) %>%
  group_by(gene_id) %>% tally() ## how many isoforms does each predicted lncrna have?
lncrna_iso$n %>% mean() ## the mean number of isoforms of lncRNAs, 1.75 

nonlncrna_iso <- fpkm[!(fpkm$gene_id %in% unique_lncrna),] %>% dplyr::select(gene_id, X) %>%
  group_by(gene_id) %>% tally() ## how many isoforms does each predicted lncrna have?
nonlncrna_iso$n %>% mean() #1.14
wilcox.test(lncrna_iso$n, nonlncrna_iso$n, alternative = "g") #significant with p < 2.2e-16


## resting lncRNAs 
# lncrnas expressed at control levels 
ylnc_control <- y.lnc.exp[apply(y.lnc.exp[,1:4] > 1, 1, any),] #345
slnc_control <- s.lnc.exp[apply(s.lnc.exp[,16:19] > 1, 1, any),] #375


## what about those that have expression above 0 but below 1...in ant,,,
ylnc_control0 <- y.lnc.exp[apply((y.lnc.exp[,1:4] > 0 & y.lnc.exp[,1:4] < 1 ), 1, any),] 
ylnc_control0 <- ylnc_control0[apply(ylnc_control0[,1:4] < 1, 1, all),] #329
slnc_control0 <- s.lnc.exp[apply((s.lnc.exp[,16:19] > 0 & s.lnc.exp[,16:19] < 1 ), 1, any),] 
slnc_control0 <- slnc_control0[apply(slnc_control0[,16:19] < 1, 1, all),]


## what is the total of the "control" lncRNAs?
c(rownames(ylnc_control), rownames(slnc_control)) %>% unique() %>% length()

# lncrnas that were not expressed in control at all...but are expressed in D1!!!!
ylnc_D1 <- y.lnc.exp[apply(y.lnc.exp[,1:4] < 1, 1, all) &
                       apply(y.lnc.exp[,5:7] > 1, 1, any) &
                       apply(y.lnc.exp[,8:11] < 1, 1, all) &
                       apply(y.lnc.exp[,12:15] < 1, 1, all),]

ylnc_D1.0 <- y.lnc.exp[apply(y.lnc.exp[,1:4] < 1, 1, all) &
                         apply(y.lnc.exp[,5:7] > 0, 1, any) &
                         apply(y.lnc.exp[,8:11] < 1, 1, all) &
                         apply(y.lnc.exp[,12:15] < 1, 1, all),]
ylnc_D1.0 <- ylnc_D1.0[apply(ylnc_D1.0[,5:7] < 1, 1, all),] #237

slnc_D1 <- s.lnc.exp[apply(s.lnc.exp[,16:19] < 1, 1, all) &
                       apply(s.lnc.exp[,20:22] > 1, 1, any) &
                       apply(s.lnc.exp[,23:26] < 1, 1, all) &
                       apply(s.lnc.exp[,27:31] < 1, 1, all),]
slnc_D1.0 <- s.lnc.exp[apply(s.lnc.exp[,16:19] < 1, 1, all) &
                         apply(s.lnc.exp[,20:22] > 0, 1, any) &
                         apply(s.lnc.exp[,23:26] < 1, 1, all) &
                         apply(s.lnc.exp[,27:31] < 1, 1, all),]
slnc_D1.0 <- slnc_D1.0[apply(slnc_D1.0[,20:22] < 1, 1, all),]


# lncrnas that were not expressed in control AND D1!!! but are expressed in WW2
ylnc_WW2 <- y.lnc.exp[apply(y.lnc.exp[,1:4] < 1, 1, all) &
                        apply(y.lnc.exp[,5:7] < 1, 1, all) &
                        apply(y.lnc.exp[,8:11] > 1, 1, any) &
                        apply(y.lnc.exp[,12:15] < 1, 1, all),]

ylnc_WW2.0 <- y.lnc.exp[apply(y.lnc.exp[,1:4] < 1, 1, all) &
                          apply(y.lnc.exp[,5:7] < 1, 1, all) &
                          apply(y.lnc.exp[,8:11] > 0, 1, any) &
                          apply(y.lnc.exp[,12:15] < 1, 1, all),]
ylnc_WW2.0 <- ylnc_WW2.0[apply(ylnc_WW2.0[,8:11] < 1, 1, all),]

slnc_WW2 <- s.lnc.exp[apply(s.lnc.exp[,16:19] < 1, 1, all) &
                        apply(s.lnc.exp[,20:22] < 1, 1, all) &
                        apply(s.lnc.exp[,23:26] > 1, 1, any) &
                        apply(s.lnc.exp[,27:31] < 1, 1, all),]
slnc_WW2.0 <- s.lnc.exp[apply(s.lnc.exp[,16:19] < 1, 1, all) &
                          apply(s.lnc.exp[,20:22] < 1, 1, all) &
                          apply(s.lnc.exp[,23:26] > 0, 1, any) &
                          apply(s.lnc.exp[,27:31] < 1, 1, all),]
slnc_WW2.0 <- slnc_WW2.0[apply(slnc_WW2.0[,23:26] < 1, 1, all),]


# lncrnas that were not expressed in control, D1 AND WW2, but are expressed in D2
ylnc_D2 <- y.lnc.exp[apply(y.lnc.exp[,1:4] < 1, 1, all) &
                       apply(y.lnc.exp[,5:7] < 1, 1, all) &
                       apply(y.lnc.exp[,8:11] < 1, 1, all) &
                       apply(y.lnc.exp[,12:15] > 1, 1, any),]

ylnc_D2.0 <- y.lnc.exp[apply(y.lnc.exp[,1:4] < 1, 1, all) &
                         apply(y.lnc.exp[,5:7] < 1, 1, all) &
                         apply(y.lnc.exp[,8:11] < 1, 1, all) &
                         apply(y.lnc.exp[,12:15] > 0, 1, any),]
ylnc_D2.0 <- ylnc_D2.0[apply(ylnc_D2.0[,12:15] < 1, 1, all),]

slnc_D2 <- s.lnc.exp[apply(s.lnc.exp[,16:19] < 1, 1, all) &
                       apply(s.lnc.exp[,20:22] < 1, 1, all) &
                       apply(s.lnc.exp[,23:26] < 1, 1, all) &
                       apply(s.lnc.exp[,27:31] > 1, 1, any),]
slnc_D2.0 <- s.lnc.exp[apply(s.lnc.exp[,16:19] < 1, 1, all) &
                         apply(s.lnc.exp[,20:22] < 1, 1, all) &
                         apply(s.lnc.exp[,23:26] < 1, 1, all) &
                         apply(s.lnc.exp[,27:31] > 0, 1, any),]
slnc_D2.0 <- slnc_D2.0[apply(slnc_D2.0[,27:31] < 1, 1, all),]



all_lncrna_codes <- code_expr[code_expr$Row.names %in% unique_lncrna,] 
imperf_lncrna <- all_lncrna_codes %>% dplyr::filter(., (class_code != "u" & class_code != "=")) 
imperf_lncrna$Row.names %in% rownames(ylnc_only)
imperf_lncrna$Row.names %in% rownames(slnc_only)

## all transcripts ecotype specific
fpkmunique[apply(fpkmunique[,1:15] > 1, 1, any) & 
             apply(fpkmunique[,16:31] < 1, 1, all),] %>% dim()

fpkmunique[apply(fpkmunique[,16:31] > 1, 1, any) & 
             apply(fpkmunique[,1:15] < 1, 1, all),] %>% dim()  


yin <- read.csv("./data_files/yin_simopoulos_overlap.txt", header=F)
yin <- fpkmunique[rownames(fpkmunique) %in% yin$V1,]
drought <- fpkmunique[substr(rownames(fpkmunique), 1,3) == "DRO",]
drought <- drought[apply(drought > 1, 1, any),]
nondrought <- fpkmunique[substr(rownames(fpkmunique), 1,3) != "DRO",]

allgenes <- fpkmunique

all.libs <- colnames(allgenes)
yuk.lib <- all.libs[which(substr(all.libs,1,1) == "Y" )]
shan.lib <- all.libs[which(substr(all.libs,1,1) == "S" )]

dim(fpkmunique[apply(fpkmunique[,1:15] > 1, 1, any) & apply(fpkmunique[,16:31] < 1, 1, all),])



for (x in all.libs){
  print(x)
}
for (x in all.libs){
  len <- length(which(fpkmunique[,x] > 1))
  print(paste(x, len))}


for (x in all.libs){
  xlocs <- fpkmunique[substr(rownames(fpkmunique), 1,4) == "XLOC",]
  len <- length(which(xlocs[,x] > 1))
  print(paste(x, len))}

for (x in all.libs){
  droughts <- fpkmunique[substr(rownames(fpkmunique), 1,7) == "DROUGHT",]
  len <- length(which(droughts[,x] > 1))
  print(paste(x, len))}

## num of Yin genes 
for (x in all.libs){
  len <- length(which(yin[,x] > 1))
  print(paste(x, len))
}


## num of lncRNAs needs ecotype
for (x in yuk.lib){
  ylncs <- fpkmunique[rownames(fpkmunique) %in% y.lnc,]
  len <- length(which(ylncs[,x] > 1))
  print(paste(x, len))}

for (x in shan.lib){
  ylncs <- fpkmunique[rownames(fpkmunique) %in% s.lnc,]
  len <- length(which(ylncs[,x] > 1))
  print(paste(x, len))}




for (x in yuk.lib){
  xlocs <- fpkmunique[substr(rownames(fpkmunique), 1,4) == "XLOC",]
  ylncs <- xlocs[rownames(xlocs) %in% y.lnc,]
  len <- length(which(ylncs[,x] > 1))
  print(paste(x, len))}

for (x in yuk.lib){
  xlocs <- fpkmunique[substr(rownames(fpkmunique), 1,7) == "DROUGHT",]
  ylncs <- xlocs[rownames(xlocs) %in% y.lnc,]
  len <- length(which(ylncs[,x] > 1))
  print(paste(x, len))}




for (x in shan.lib){
  xlocs <- fpkmunique[substr(rownames(fpkmunique), 1,4) == "XLOC",]
  slncs <- xlocs[rownames(xlocs) %in% s.lnc,]
  len <- length(which(slncs[,x] > 1))
  print(paste(x, len))}

for (x in shan.lib){
  xlocs <- fpkmunique[substr(rownames(fpkmunique), 1,7) == "DROUGHT",]
  slncs <- xlocs[rownames(xlocs) %in% s.lnc,]
  len <- length(which(slncs[,x] > 1))
  print(paste(x, len))}


dim(fpkmunique[apply(fpkmunique[,1:15] > 1, 1, any) & apply(fpkmunique[,16:31] < 1, 1, all),])
dim(fpkmunique[apply(fpkmunique[,1:15] < 1, 1, all) & apply(fpkmunique[,16:31] > 1, 1, any),])

fpkmunique[apply(fpkmunique > 1,1, any),] %>% dim() 
xlocs[apply(xlocs > 1,1, any),] %>% dim() 
droughts[apply(droughts > 1,1, any),] %>% dim() 

fpkmunique[apply(fpkmunique[,shan.lib] > 1,1, all),] %>% dim() 

fpkmunique[apply(fpkmunique > 1, 1, all) %>% names(.) %in% s.lnc & 
             rownames(fpkmunique) %>% substr(., 1,1) == "T",] %>% dim()
fpkmunique[apply(fpkmunique > 1, 1, all) %>% names(.) %in% y.lnc,] %>% dim()



