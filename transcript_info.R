library(tidyverse)
setwd("~/iclouddrive/Documents/Manuscripts/drought/eutrema_drought_transcriptomics/")
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

## Exporting predicted lncRNAs for FileS2

lncinfo <- fpkm[,c(1,2,3,5)]

## for FileS2
lncinfo <- lncinfo %>% filter((y_lnc == 1 | s_lnc == 1)) %>% select(-X) %>% distinct()




## predicted lncRNAs in each ecotype
y.lnc <- rownames(unique_fpkm[which(unique_fpkm[,"y_lnc", drop=F] == 1),])
s.lnc <- rownames(unique_fpkm[which(unique_fpkm[,"s_lnc", drop=F] == 1),])

# same lncRNA in each ecotype
samepredictlncrna <- intersect(y.lnc, s.lnc)
fpkmunique <- unique_fpkm[,c(9:15, 17:21, 23:25, 29:35, 37:42, 44:46)]
colnames(fpkmunique)[c(1,2,8,9,16,17,23,24)] <- c("YWW1.1", "YWW1.2", "YWW2.1", "YWW2.2",
                                                  "SWW1.1", "SWW1.2", "SWW2.1", "SWW2.2")


## lncRNAs actually expressed 
## work on this 
ylncexp <- lncinfo %>% select(gene_id, y_lnc) %>% merge(., fpkmunique, by.x = "gene_id", by.y = 0) %>%
  select(-starts_with("S"))
ylncexp <- ylncexp[apply(ylncexp[,3:17] > 1, 1, any),]
slncexp <- lncinfo %>% select(gene_id, s_lnc) %>% merge(., fpkmunique, by.x = "gene_id", by.y = 0) %>%
  select(-starts_with("Y"))
slncexp <- slncexp[apply(slncexp[,3:18] > 1, 1, any),]

## lncrnas prediected and expressed in their appropriate ecotype
lncexp_ys <- c(as.character(ylncexp$gene_id), as.character(slncexp$gene_id)) %>% unique()
lncexp_ys <- merge(lncinfo, data.frame(gene_id=lncexp_ys), all.y=T) %>% unique()
## the issue is that there are lncRNAs predicted differently for eahch ecotype...
## but since we are using gene level info, we should change this to just be predicted as lncrnas 
## in Y and S
## example:
lncexp_ys %>% group_by(gene_id) %>% filter(n() > 1)

## remmove the duplicated ones 
dup <- lncexp_ys %>% group_by(gene_id) %>% filter(n() > 1) %>% select(gene_id) %>% unique() 
lncexp_ys <- lncexp_ys %>% mutate(y_lnc = case_when(gene_id %in% dup$gene_id ~ "1",
                                            TRUE ~ as.character(y_lnc)),
                          s_lnc = case_when(gene_id %in% dup$gene_id ~ "1",
                                            TRUE ~ as.character(s_lnc))) %>% unique()
lncexp_ys %>% filter(y_lnc == 1 & s_lnc == 1) %>% dim()
grep("Thhalv", lncexp_ys$gene_id) %>% length()
76/1007 *100  # percentage of Thhalv lncrnas
#write.table(lncexp_ys, "./data_files/lnc_expressed_out.txt", sep="\t", quote =F, col.names = F)

## redo the y and s lnc dataframes using updated dataframe...
ylncexp <- lncexp_ys %>% filter(y_lnc ==1) %>% merge(., fpkmunique, by.x = "gene_id", by.y = 0) %>%
  select(-starts_with("S"))

slncexp <- lncexp_ys %>% filter(s_lnc ==1) %>% merge(., fpkmunique, by.x = "gene_id", by.y = 0) %>%
  select(-starts_with("Y"))


## all lncRNAs predicted, no matter their expression patterns
y.lnc.exp <- fpkmunique[rownames(fpkmunique) %in% y.lnc,] 
s.lnc.exp <- fpkmunique[rownames(fpkmunique) %in% s.lnc,] 


## lncrnas ecotype specific, but have no consideration for expression
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
#ylnc_1 <- y.lnc.exp[apply(y.lnc.exp[,1:15] > 1, 1, any),] 
#slnc_1 <- s.lnc.exp[apply(s.lnc.exp[,16:31] > 1, 1, any),]
#unique_lncrna <- c(rownames(ylnc_1), rownames(slnc_1)) %>% base::unique()
#both_lncrna <- intersect(rownames(ylnc_1),rownames(slnc_1))

## lncRNAs that are expressed only in stress treated libraries 
#unique_lncrna_expr <- fpkmunique[rownames(fpkmunique) %in% unique_lncrna,]

## updated code 2020-04-03 
unique_lncrna_expr <- fpkmunique[rownames(fpkmunique) %in% lncexp_ys$gene_id,]

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
##updated code 2020-04-03
## ylncexp
## slncexp 
ylnc_control <- ylncexp[apply(ylncexp[,3:6] > 1, 1, any),] #627
slnc_control <- slncexp[apply(slncexp[,3:6] > 1, 1, any),] #704


## stopped here 
## what about those that have expression above 0 but below 1...in ant,,,
ylnc_control0 <- ylncexp[apply(( ylncexp[,3:6] > 0 &  ylncexp[,3:6] < 1 ), 1, any),] 
ylnc_control0 <- ylnc_control0[apply(ylnc_control0[,3:6] < 1, 1, all),] #329
slnc_control0 <- slncexp[apply((slncexp[,3:6] > 0 & slncexp[,3:6] < 1 ), 1, any),] 
slnc_control0 <- slnc_control0[apply(slnc_control0[,3:6] < 1, 1, all),]


## what is the total of the "control" lncRNAs?
c(rownames(ylnc_control), rownames(slnc_control)) %>% unique() %>% length()


### function for finding lncrna expression at specific treatments and ecptypes
#lnccond_eco <- function(ecotype, condition){
#  if (ecotype == "Yukon"){
#    if (condition == "WW1"){
#      
#    } elif (condition == "D1") {
#      
#    } elif (condition == "WW2") {
#      
#    } elif (condition == "D2") {
#      
#    }
#   cond_exp <- ylncexp[apply(ylncexp[,3:6] < 1, 1, all) &
#              apply(ylncexp[,7:9] > 1, 1, any) &
#              apply(ylncexp[,10:13] < 1, 1, all) &
#              apply(ylncexp[,14:17] < 1, 1, all),]
#    
#    #d1 failure
#    ylnc_D1.0 <- ylncexp[apply(ylncexp[,3:6] < 1, 1, all) &
#                           apply(ylncexp[,7:9] > 0, 1, any) &
#                           apply(ylncexp[,10:13] < 1, 1, all) &
#                           apply(ylncexp[,14:17] < 1, 1, all),]
#    ylnc_D1.0 <- ylnc_D1.0[apply(ylnc_D1.0[,7:9] < 1, 1, all),] #237
#  } else { ##ecotype == "Shandong"
#    
#    
#  }
#  
#}


# lncrnas that were not expressed in control at all...but are expressed in D1!!!!
## ww1 3:6
## d1 7:9
## ww2 10:13
## d2 14:17
ylnc_D1 <- ylncexp[apply(ylncexp[,3:6] < 1, 1, all) &
                       apply(ylncexp[,7:9] > 1, 1, any) &
                       apply(ylncexp[,10:13] < 1, 1, all) &
                       apply(ylncexp[,14:17] < 1, 1, all),]

#d1 failure
ylnc_D1.0 <- ylncexp[apply(ylncexp[,3:6] < 1, 1, all) &
                         apply(ylncexp[,7:9] > 0, 1, any) &
                         apply(ylncexp[,10:13] < 1, 1, all) &
                         apply(ylncexp[,14:17] < 1, 1, all),]
ylnc_D1.0 <- ylnc_D1.0[apply(ylnc_D1.0[,7:9] < 1, 1, all),] #237


## ww1 3:6
## d1 7:9
## ww2 10:13
## d2 14:18
slnc_D1 <- slncexp[apply(slncexp[,3:6] < 1, 1, all) &
                       apply(slncexp[,7:9] > 1, 1, any) &
                       apply(slncexp[,10:13] < 1, 1, all) &
                       apply(slncexp[,14:18] < 1, 1, all),]
slnc_D1.0 <- slncexp[apply(slncexp[,3:6] < 1, 1, all) &
                         apply(slncexp[,7:9] > 0, 1, any) &
                         apply(slncexp[,10:13] < 1, 1, all) &
                         apply(slncexp[,14:18] < 1, 1, all),]
slnc_D1.0 <- slnc_D1.0[apply(slnc_D1.0[,7:9] < 1, 1, all),]


# lncrnas that were not expressed in control AND D1!!! but are expressed in WW2
## y
## ww1 3:6
## d1 7:9
## ww2 10:13
## d2 14:17

# s
## ww1 3:6
## d1 7:9
## ww2 10:13
## d2 14:18
ylnc_WW2 <- ylncexp[apply(ylncexp[,3:6] < 1, 1, all) &
                        apply(ylncexp[,7:9] < 1, 1, all) &
                        apply(ylncexp[,10:13] > 1, 1, any) &
                        apply(ylncexp[,14:17] < 1, 1, all),]

ylnc_WW2.0 <- ylncexp[apply(ylncexp[,3:6] < 1, 1, all) &
                          apply(ylncexp[,7:9] < 1, 1, all) &
                          apply(ylncexp[,10:13] > 0, 1, any) &
                          apply(ylncexp[,14:17] < 1, 1, all),]
ylnc_WW2.0 <- ylnc_WW2.0[apply(ylnc_WW2.0[,10:13] < 1, 1, all),]

slnc_WW2 <- slncexp[apply(slncexp[,3:6] < 1, 1, all) &
                        apply(slncexp[,7:9] < 1, 1, all) &
                        apply(slncexp[,10:13] > 1, 1, any) &
                        apply(slncexp[,14:18] < 1, 1, all),]
slnc_WW2.0 <- slncexp[apply(slncexp[,3:6] < 1, 1, all) &
                          apply(slncexp[,7:9] < 1, 1, all) &
                          apply(slncexp[,10:13] > 0, 1, any) &
                          apply(slncexp[,14:18] < 1, 1, all),]
slnc_WW2.0 <- slnc_WW2.0[apply(slnc_WW2.0[,10:13] < 1, 1, all),]


# lncrnas that were not expressed in control, D1 AND WW2, but are expressed in D2
## y
## ww1 3:6
## d1 7:9
## ww2 10:13
## d2 14:17

# s
## ww1 3:6
## d1 7:9
## ww2 10:13
## d2 14:18

ylnc_D2 <- ylncexp[apply(ylncexp[,3:6] < 1, 1, all) &
                       apply(ylncexp[,7:9] < 1, 1, all) &
                       apply(ylncexp[,10:13] < 1, 1, all) &
                       apply(ylncexp[,14:17] > 1, 1, any),]

ylnc_D2.0 <- ylncexp[apply(ylncexp[,3:6] < 1, 1, all) &
                         apply(ylncexp[,7:9] < 1, 1, all) &
                         apply(ylncexp[,10:13] < 1, 1, all) &
                         apply(ylncexp[,14:17] > 0, 1, any),]
ylnc_D2.0 <- ylnc_D2.0[apply(ylnc_D2.0[,14:17] < 1, 1, all),]

slnc_D2 <- slncexp[apply(slncexp[,3:6] < 1, 1, all) &
                       apply(slncexp[,7:9] < 1, 1, all) &
                       apply(slncexp[,10:13] < 1, 1, all) &
                       apply(slncexp[,14:18] > 1, 1, any),]
slnc_D2.0 <- slncexp[apply(slncexp[,3:6] < 1, 1, all) &
                         apply(slncexp[,7:9] < 1, 1, all) &
                         apply(slncexp[,10:13] < 1, 1, all) &
                         apply(slncexp[,14:18] > 0, 1, any),]
slnc_D2.0 <- slnc_D2.0[apply(slnc_D2.0[,13:18] < 1, 1, all),]



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



