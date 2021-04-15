# OSA_RNAseq
This repository contains the scripts used for processing canine osteosarcoma RNAseq data in the proposed manuscript titled "Transcriptomic Analysis of Canine Osteosarcoma".

The scripts were used sequentially as follows:  
1_Quality.sh :: assess quality of raw sequence data using FastQC  
2_Trim.sh :: trim reads using Trimmomatic and re-assess the quality using FastQC  
3_Map.sh :: map the reads onto the reference genome using HiSat2 and count the mapped reads using Stringtie   
