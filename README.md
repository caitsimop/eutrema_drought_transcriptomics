# R scripts for the analysis of our RNAseq data
Code to accompany "Coding and long non-coding RNAs provide evidence of distinct transcriptional reprogramming for two ecotypes of the extremophile plant Eutrema salsugineum undergoing water deficit stress"

See our manuscript in BMC Genomics [here](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-06793-7)



## Information on files
-  `data_files`: files needed for analysis are in this directory
- `deseq2_drought`: contains code for identifying differentially expressed genes and fold change figure
- `drought_wgcna.R`: contains code for WGCNA analysis, GO enrichment and heatmap figure
- `pca.R`: contains code for all PCAs and accompanying biplots
- `transcript_info.R`: contains code for identifying various stats about the transcripts (ie how many are only expressed in a certain condition).
