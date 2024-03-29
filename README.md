# OSA_RNAseq
This repository contains the scripts used for processing canine osteosarcoma RNAseq data in the manuscript titled "Transcriptomic Analysis of Canine Osteosarcoma from a Precision Medicine Perspective Reveals Limitations of Differential Gene Expression Studies" (Nance et al, 2022). 

The scripts were used sequentially as follows:  
1_Quality.sh :: assess quality of raw sequence data using FastQC  
2_Trim.sh :: trim reads using Trimmomatic and re-assess the quality using FastQC  
3_Map.sh :: map the reads onto the reference genome using HiSat2 and count the mapped reads using Stringtie   
PrepDE.py :: used in 3_Map.sh to convert files in ballgown folder to a count matrix     
4_DESeq2 :: R script to analyze DEGs using DESeq2, prepare data for GSEA and Cytoscape, and generate figures/plots/files for further analysis and publication

The data files include:  
gene_count_matrix :: this is the raw counts matrix generated from Stringtie    
res_all :: list of all differentially expressed genes outputted from DEseq2    
res_sig :: list of significant differentially expressed genes, as defined by padj<0.05, base mean<10, log2FC>1,<-1    
res_indiv_prefilter :: list of log2FC values for each dog for every gene listed in res_sig     
