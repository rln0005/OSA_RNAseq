
###############  Differential Gene Expression Analysis for Canine OSA RNAseq data  ###############
###########################  Written By: Rebecca L. Nance, 08/24/2021  ###########################  

#set working directory 
#load required programs
library(DESeq2)
library(genefilter)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(ellipsis)
library(biomaRt)
library(dbplyr)
library(tidyverse)
library(org.Cf.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(ggplotify)
library(ggforce)
library(gg.gap)

###################################  DESeq2 Analysis  ###################################  
#import count data from pythonDE ballgown -- 30951 genes
countdata <- as.matrix(read.csv("gene_count_matrix_NEW.csv", row.names="gene_id"))
#import meta data file
coldata <-(read.table("pheno_NEW.txt", header=TRUE, row.names=1))
#check all sample IDs in colData are also in CountData and match their orders
all(rownames(coldata) %in% colnames(countdata))
all(rownames(coldata) == colnames(countdata))
#create DEseq2 dataset
dds <- DESeqDataSetFromMatrix(countData = countdata, colData=coldata,  design = ~Dog + Source)
#pre-filter to exclude genes with less than 1 count
dds <- dds[ rowSums(counts(dds)) > 1, ]
#differential expression analysis
dds <- DESeq(dds)

###################################  Extract DESeq2 Results  ###################################  
#extract results
res <- results(dds, alpha=0.05)
#how many significant genes with p<0.05? 
sum(res$padj < 0.05, na.rm=TRUE)
#summary of results
summary(res)
#create a tibble of ALL of the results with a column labeled "gene" containing ENSEMBL IDs
res_all <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
#create a tibble of the SIGNIFICANT results (padj<0.05, baseMean>10, FC>2,<-2) with a column labeled "gene" containing ENSEMBL IDs
res_sig <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  filter(padj < 0.05) %>%
  filter(baseMean > 10) %>%
  filter(log2FoldChange > 1 | log2FoldChange < -1) %>%
  as_tibble()

###################################  Annotate ENSEMBL IDs with gene symbols  ###################################  
#extract a list of ALL ensembl gene IDs from the analysis
genelist <- rownames(countdata) %>%
  data.frame()
#create annotation library
annotations <- AnnotationDbi::select(org.Cf.eg.db,
                                         keys = genelist[,1], 
                                         columns = c("SYMBOL", "GENENAME"),
                                         keytype = "ENSEMBL") 
#determine indices for NON-DUPLICATED genes and return only the non-duplicated genes 
non_dupl_idx <- which(duplicated(annotations$ENSEMBL) == FALSE)
annotations <- annotations[non_dupl_idx, ]
#determine indices for NON-NA genes and return only the genes with gene symbol annotations--not using this
  #non_na_idx <- which(is.na(annotations$SYMBOL) == FALSE)
  #annotations_no_na <- annotations[non_na_idx, ]
#merge together the significant results with the annotated symbols for both results
res_all <- merge(res_all, annotations, by.x="gene", by.y="ENSEMBL")
res_sig <- merge(res_sig, annotations, by.x="gene", by.y="ENSEMBL")
#reorder the significant results based on smallest padj
res_all <- res_all[order(res_all$padj),]
res_sig <- res_sig[order(res_sig$padj),]
#write results to a file
write.csv(as.data.frame(res_all), file="res_all_NEW.csv", row.names=FALSE)
write.csv(as.data.frame(res_sig), file="res_sig_NEW.csv", row.names=FALSE)

###################################  Making Figures for Paper  ###################################  
## Volcano plot--FIGURE 2A
ylab <- expression(paste(-Log[10], "(adjusted p-value)"))
volcano <- EnhancedVolcano(res_all,
                           lab=res_all$SYMBOL,
                           x='log2FoldChange',
                           y='padj',
                           title="Volcano Plot of Differentially Expressed Genes in Tumor vs Bone",
                           ylim=c(0,35),
                           xlim=c(-11,13),
                           titleLabSize=12,
                           axisLabSize=10,
                           legendLabSize=10,
                           legendIconSize=4,
                           subtitleLabSize=0.01,
                           captionLabSize=6,
                           labSize=2.5, 
                           subtitle= "",
                           ylab=ylab,
                           pCutoff=0.05,
                           gridlines.major=FALSE,
                           gridlines.minor=FALSE,
                           pointSize=0.6,
                           legendLabels=c('Not sig.', 'log2(FC)', 'p-value', 'p-value & log2(FC)'))
## Heatmap--FIGURE 2B
vst <- varianceStabilizingTransformation(dds, blind=FALSE)
vstDF <- assay(vst)
vstDF <- vstDF %>%
  data.frame() %>%
  rownames_to_column(var="gene")
vstDF <- vstDF %>%
  filter(gene %in% res_sig$gene)
rownames(vstDF) <- vstDF$gene  
vstDF <- vstDF[,2:15]
label_col <- c("1-Bone","2-Bone","3-Bone","4-Bone","5-Bone","6-Bone","7-Bone","1-Tumor","2-Tumor","3-Tumor","4-Tumor","5-Tumor","6-Tumor","7-Tumor")
colnames(vstDF) <- label_col
anno <-(read.table("heatmap_anno_NEW.txt", header=TRUE, row.names=1))
rownames(anno) <- colnames(vstDF)
heat_col = list(Source = c(Bone = "firebrick", Tumor="darkblue"))
heatmap <- as.ggplot(pheatmap(vstDF,
                              labels_col=label_col,
                              annotation_col=anno,
                              annotation_names_col=FALSE,
                              annotation_colors=heat_col,
                              color = colorRampPalette(c("blue2", "white", "red2"))(50),
                              show_rownames=FALSE,
                              cluster_rows=TRUE,
                              cluster_cols=TRUE,
                              scale="row",
                              border_color=NA,
                              treeheight_row = 0,
                              treeheight_col= 15,
                              fontsize=7,
                              fontsize_col=8,
                              legend=T,
                              main="Heatmap of Significant DEGs in Tumor vs Normal",
                              annotation_legend=FALSE,
                              silent=TRUE))
## PCA plot--FIGURE 2C
PCAdata <- plotPCA(vst, intgroup=c("Source"), returnData=TRUE)
percentVar <- round(100 * attr(PCAdata, "percentVar"))
PCA <- ggplot(PCAdata, aes(PC1, PC2, color=Source))+
  geom_point(size=2) +
  ggtitle("Principal Component Analysis")+
  xlab(paste0("PC1: ",percentVar[1], "% variance"))+
  ylab(paste0("PC2: ",percentVar[2], "% variance"))+
  theme(legend.position="top", legend.title=element_blank(), axis.title=element_text(size=9), plot.title=element_text(size=8.5, face="bold"))
### Combine the Volcano plot, Heatmap, and PCA plot in one figure--FIGURE 2
plot<- ggarrange(volcano,
                 ggarrange(heatmap, PCA, ncol=2, widths=c(2,1), labels=c("B", "C")),
                 nrow=2, heights=c(2,1),
                 labels="A")
ggsave("Figure2.png", plot, width=180, height=180, units="mm", dpi=600)
ggsave("Figure2.pdf", plot, width=180, height=180, units="mm", dpi=600)

### Plot individual genes--FIGURE 3
## Plot the #1 upregulated gene (with biggest LFC) from group data--FIGURE 3A
upplot1 <- plotCounts(dds, normalized=TRUE, gene="ENSCAFG00000044295", intgroup =c("Source", "Dog"), returnData=TRUE)
upplot1 <- upplot1 %>%
  ggplot(aes(x=Source,y=count,color=Dog)) +
  geom_line(aes(group=Dog), position=position_dodge(0.3), size=0.2)+
  geom_point(size=0.3, position=position_dodge(0.3))+
  labs(title="#1 Upregulated Gene (biggest LFC): \nENSCAFG00000044295", y="Normalized Counts")+
  annotate("text", x=1.5, y=Inf, vjust=1.3, parse=TRUE, size=1.5, label= "'FDR=1.3'%*%10^-4")+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(color="black"),
        panel.border=element_rect(fill=NA),
        plot.title=element_text(hjust=0.5, size=5),
        axis.title=element_text(size=4),
        axis.text.y=element_text(size=4),
        axis.text.x=element_text(size=4),
        axis.title.x=element_blank(),
        legend.text=element_text(size=5),
        legend.title=element_text(size=6),
        legend.key.size=unit(1.5, 'mm'))
## Plot the #1 downregulated gene (with biggest LFC) from group data--FIGURE 3C
downplot1 <-plotCounts(dds, normalized=TRUE, gene="ENSCAFG00000028799", intgroup=c("Source", "Dog"), returnData=TRUE)
downplot1 <- downplot1 %>%
  ggplot(aes(x=Source,y=count,color=Dog)) +
  geom_line(aes(group=Dog), position=position_dodge(0.3), size=0.2)+
  geom_point(size=0.3, position=position_dodge(0.3))+
  labs(title="#1 Downregulated Gene (biggest LFC): \nARHGEF1", y="Normalized Counts")+
  annotate("text", x=1.5, y=Inf, vjust=1.3, parse=TRUE, size=1.5, label= "'FDR=2.9'%*%10^-3")+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(color="black"),
        panel.border=element_rect(fill=NA),
        plot.title=element_text(hjust=0.5, size=5),
        axis.title=element_text(size=4),
        axis.text.y=element_text(size=4),
        axis.text.x=element_text(size=4),
        axis.title.x=element_blank(),
        legend.text=element_text(size=5),
        legend.title=element_text(size=6),
        legend.key.size=unit(1.5, 'mm'))
## Plot the #1 upregulated gene (with smallest padj) from group data--FIGURE 3B
upplot2 <- plotCounts(dds, normalized=TRUE, gene="ENSCAFG00000000782", intgroup =c("Source", "Dog"), returnData=TRUE)
upplot2 <- upplot2 %>%
  ggplot(aes(x=Source,y=count,color=Dog)) +
  geom_line(aes(group=Dog), position=position_dodge(0.3), size=0.2)+
  geom_point(size=0.3, position=position_dodge(0.3))+
  labs(title="#1 Upregulated Gene (smallest padj): \nGTSE1", y="Normalized Counts")+
  annotate("text", x=1.5, y=Inf, vjust=1.3, parse=TRUE, size=1.5, label= "'FDR=8.5'%*%10^-31")+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(color="black"),
        panel.border=element_rect(fill=NA),
        plot.title=element_text(hjust=0.5, size=5),
        axis.title=element_text(size=4),
        axis.text.y=element_text(size=4),
        axis.text.x=element_text(size=4),
        axis.title.x=element_blank(),
        legend.text=element_text(size=5),
        legend.title=element_text(size=6),
        legend.key.size=unit(1.5, 'mm'))
## Plot the #1 downregulated gene (with smallest padj) from group data--FIGURE 3D
downplot2 <-plotCounts(dds, normalized=TRUE, gene="ENSCAFG00000011986", intgroup=c("Source", "Dog"), returnData=TRUE)
downplot2 <- downplot2 %>%
  ggplot(aes(x=Source,y=count,color=Dog)) +
  geom_line(aes(group=Dog), position=position_dodge(0.3), size=0.2)+
  geom_point(size=0.3, position=position_dodge(0.3))+
  labs(title="#1 Downregulated Gene (smallest padj): \nPLIN1", y="Normalized Counts")+
  annotate("text", x=1.5, y=Inf, vjust=1.3, parse=TRUE, size=1.5, label= "'FDR=1.9'%*%10^-35")+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(color="black"),
        panel.border=element_rect(fill=NA),
        plot.title=element_text(hjust=0.5, size=5),
        axis.title=element_text(size=4),
        axis.text.y=element_text(size=4),
        axis.text.x=element_text(size=4),
        axis.title.x=element_blank(),
        legend.text=element_text(size=5),
        legend.title=element_text(size=6),
        legend.key.size=unit(1.5, 'mm'))
## Combine all the plots into one figure--FIGURE 3
comboplot <- ggarrange(upplot1, upplot2, downplot1, downplot2, labels=c("A","B","C","D"), font.label=list(color="black", size=8),ncol=2, nrow=2, align=c("hv"))
comboplot <- annotate_figure(comboplot, top=text_grob("Top Genes from Group Analysis", face="bold", size=8))
ggsave("Figure4.png", comboplot, width=90, height=80, units="mm", dpi=600)
ggsave("Figure4.pdf", comboplot, width=90, height=80, units="mm", dpi=600)

###################################  Individual Analysis  ###################################  
### Individual FC Pre-filtered Results
#extract variance stabilized transformation 
vst <- varianceStabilizingTransformation(dds, blind=FALSE)
vst <- assay(vst)
vst <- vst %>%
  as.data.frame() %>%
  rownames_to_column("gene")
#filter genes based on significant DEGs (padj<0.05, baseMean>10, FC>2,<-2)
vst_filter <- vst %>%
  as.data.frame() %>%
  filter(gene %in% res_sig$gene)
#generate L2FC values for each patient (since VST is already on log2 scale, subtract them instead of divide)
pt1 <- as.data.frame(vst_filter[,9]-vst_filter[,2], rownames<- vst_filter$gene)
pt2 <- as.data.frame(vst_filter[,10]-vst_filter[,3], rownames<- vst_filter$gene)
pt3 <- as.data.frame(vst_filter[,11]-vst_filter[,4], rownames<- vst_filter$gene)
pt4 <- as.data.frame(vst_filter[,12]-vst_filter[,5], rownames<- vst_filter$gene)
pt5 <- as.data.frame(vst_filter[,13]-vst_filter[,6], rownames<- vst_filter$gene)
pt6 <- as.data.frame(vst_filter[,14]-vst_filter[,7], rownames<- vst_filter$gene)
pt7 <- as.data.frame(vst_filter[,15]-vst_filter[,8], rownames<- vst_filter$gene)
res_indiv_filter <- cbind(pt1,pt2,pt3,pt4,pt5,pt6,pt7)
colnames(res_indiv_filter) <- c("pt1","pt2","pt3","pt4","pt5","pt6","pt7")
res_indiv_filter <- res_indiv_filter %>%
  data.frame() %>%
  rownames_to_column("gene")
res_indiv_filter <- merge(res_indiv_filter, annotations, by.x="gene", by.y="ENSEMBL")
res_indiv_filter[,2:8][res_indiv_filter[,2:8] == 0] <- "na"
#write results to a file
write.csv(as.data.frame(res_indiv_filter), file="res_indiv_prefilter_NEW.csv", row.names=FALSE)

###################################  Figures For Individual Analysis  ###################################  
##### FIGURE 5--Top Upregulated Gene in Each Patient 
#pt 1 top upreg gene==ENSCAFG00000041995
pt1plot <- plotCounts(dds, normalized=TRUE, gene="ENSCAFG00000041995", intgroup =c("Source", "Dog"), returnData=TRUE)
pt1plot <- pt1plot %>%
  ggplot(aes(x=Source,y=count,color=Dog)) +
  geom_line(aes(group=Dog), position=position_dodge(0.3), size=0.4)+
  geom_point(size=0.5, position=position_dodge(0.3))+
  labs(title="Patient A: \nENSCAFG00000041995", y="Normalized Counts")+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(color="black"),
        panel.border=element_rect(fill=NA),
        plot.title=element_text(hjust=0.5, size=7),
        axis.title=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6),
        axis.title.x=element_blank(),
        legend.text=element_text(size=6),
        legend.title=element_text(size=7),
        legend.key.size=unit(1.5, 'mm'))
#pt 2 top upreg gene plot==ENSCAFG00000009135 (LOC403585)
pt2plot <- plotCounts(dds, normalized=TRUE, gene="ENSCAFG00000009135", intgroup =c("Source", "Dog"), returnData=TRUE)
pt2plot <- pt2plot %>%
  ggplot(aes(x=Source,y=count,color=Dog)) +
  geom_line(aes(group=Dog), position=position_dodge(0.3), size=0.4)+
  geom_point(size=0.5, position=position_dodge(0.3))+
  labs(title="Patient B: \nLOC403585", y="Normalized Counts")+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(color="black"),
        panel.border=element_rect(fill=NA),
        plot.title=element_text(hjust=0.5, size=7),
        axis.title=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6),
        axis.title.x=element_blank(),
        legend.text=element_text(size=6),
        legend.title=element_text(size=7),
        legend.key.size=unit(1.5, 'mm'))
#pt 3,7 top upreg gene plot==ENSCAFG00000002040 (TFPI2)
pt3plot <- plotCounts(dds, normalized=TRUE, gene="ENSCAFG00000002040", intgroup =c("Source", "Dog"), returnData=TRUE)
pt3plot <- pt3plot %>%
  ggplot(aes(x=Source,y=count,color=Dog)) +
  geom_line(aes(group=Dog), position=position_dodge(0.3), size=0.4)+
  geom_point(size=0.5, position=position_dodge(0.3))+
  labs(title="Patients C, G: \nTFPI2", y="Normalized Counts")+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(color="black"),
        panel.border=element_rect(fill=NA),
        plot.title=element_text(hjust=0.5, size=7),
        axis.title=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6),
        axis.title.x=element_blank(),
        legend.text=element_text(size=6),
        legend.title=element_text(size=7),
        legend.key.size=unit(1.5, 'mm'))
#pt 4 top upreg gene plot==ENSCAFG00000019985 (COL11A1)
pt4plot <- plotCounts(dds, normalized=TRUE, gene="ENSCAFG00000019985", intgroup =c("Source", "Dog"), returnData=TRUE)
pt4plot <- pt4plot %>%
  ggplot(aes(x=Source,y=count,color=Dog)) +
  geom_line(aes(group=Dog), position=position_dodge(0.3), size=0.4)+
  geom_point(size=0.5, position=position_dodge(0.3))+
  labs(title="Patient D: \nCOL11A1", y="Normalized Counts")+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(color="black"),
        panel.border=element_rect(fill=NA),
        plot.title=element_text(hjust=0.5, size=7),
        axis.title=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6),
        axis.title.x=element_blank(),
        legend.text=element_text(size=6),
        legend.title=element_text(size=7),
        legend.key.size=unit(1.5, 'mm'))
#pt 5 top upreg gene plot==ENSCAFG00000008353 (SFRP2)
pt5plot <- plotCounts(dds, normalized=TRUE, gene="ENSCAFG00000008353", intgroup =c("Source", "Dog"), returnData=TRUE)
pt5plot <- pt5plot %>%
  ggplot(aes(x=Source,y=count,color=Dog)) +
  geom_line(aes(group=Dog), position=position_dodge(0.3), size=0.4)+
  geom_point(size=0.5, position=position_dodge(0.3))+
  labs(title="Patient E: \nSFRP2", y="Normalized Counts")+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(color="black"),
        panel.border=element_rect(fill=NA),
        plot.title=element_text(hjust=0.5, size=7),
        axis.title=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6),
        axis.title.x=element_blank(),
        legend.text=element_text(size=6),
        legend.title=element_text(size=7),
        legend.key.size=unit(1.5, 'mm'))
#pt 6 top upreg gene plot==ENSCAFG00000028460
pt6plot <- plotCounts(dds, normalized=TRUE, gene="ENSCAFG00000028460", intgroup =c("Source", "Dog"), returnData=TRUE)
pt6plot <- pt6plot %>%
  ggplot(aes(x=Source,y=count,color=Dog)) +
  geom_line(aes(group=Dog), position=position_dodge(0.3), size=0.4)+
  geom_point(size=0.5, position=position_dodge(0.3))+
  labs(title="Patient F: \nENSCAFG00000028460", y="Normalized Counts")+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(color="black"),
        panel.border=element_rect(fill=NA),
        plot.title=element_text(hjust=0.5, size=7),
        axis.title=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6),
        axis.title.x=element_blank(),
        legend.text=element_text(size=6),
        legend.title=element_text(size=7),
        legend.key.size=unit(1.5, 'mm'))
## Combine all the plots into one figure--FIGURE 5
comboplot2 <- ggarrange(pt1plot, pt2plot, pt3plot, pt4plot, pt5plot, pt6plot, labels=c("A","B","C","D","E","F"), font.label=list(color="black", size=8),ncol=3, nrow=2, align=c("hv"))
comboplot2 <- annotate_figure(comboplot2, top=text_grob("Top Upregulated Gene in Each Patient", face="bold", size=10))
ggsave("Figure5.png", comboplot2, width=180, height=100, units="mm", dpi=600)
ggsave("Figure5.pdf", comboplot2, width=180, height=100, units="mm", dpi=600)

##### FIGURE 6--Top Downregulated Gene in Each Patient 
#pt 1 top downreg gene==ENSCAFG00000015706
pt1plot <- plotCounts(dds, normalized=TRUE, gene="ENSCAFG00000015706", intgroup =c("Source", "Dog"), returnData=TRUE)
pt1plot <- pt1plot %>%
  ggplot(aes(x=Source,y=count,color=Dog)) +
  geom_line(aes(group=Dog), position=position_dodge(0.3), size=0.4)+
  geom_point(size=0.5, position=position_dodge(0.3))+
  labs(title="Patient A: \nCYTL1", y="Normalized Counts")+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(color="black"),
        panel.border=element_rect(fill=NA),
        plot.title=element_text(hjust=0.5, size=6),
        axis.title=element_text(size=4),
        axis.text.y=element_text(size=4),
        axis.text.x=element_text(size=4),
        axis.title.x=element_blank(),
        legend.text=element_text(size=5),
        legend.title=element_text(size=6),
        legend.key.size=unit(1.5, 'mm'))
#pt 2,5,7 top downreg gene plot==ENSCAFG00000034058
B <- plotCounts(dds, normalized=TRUE, gene="ENSCAFG00000034058", intgroup =c("Source", "Dog"), returnData=TRUE)
B2 <- B %>%
  ggplot(aes(x=Source,y=count,color=Dog)) +
  geom_line(aes(group=Dog), position=position_dodge(0.2), size=0.4)+
  geom_point(size=0.5, position=position_dodge(0.2))+
  labs(title="Patients B, E, G: \nENSCAFG00000034058", y="Normalized Counts")+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(color="black"),
        panel.border=element_rect(fill=NA),
        plot.title=element_text(hjust=0.5, size=6),
        axis.title=element_text(size=4),
        axis.text.y=element_text(size=4),
        axis.text.x=element_text(size=4),
        axis.title.x=element_blank(),
        legend.text=element_text(size=5),
        legend.title=element_text(size=6),
        legend.key.size=unit(1.5, 'mm'))
B2 <- B2 + facet_zoom(ylim=c(0,15000))
#pt 3 top downreg gene plot==ENSCAFG00000014924
pt3plot <- plotCounts(dds, normalized=TRUE, gene="ENSCAFG00000014924", intgroup =c("Source", "Dog"), returnData=TRUE)
pt3plot <- pt3plot %>%
  ggplot(aes(x=Source,y=count,color=Dog)) +
  geom_line(aes(group=Dog), position=position_dodge(0.3), size=0.4)+
  geom_point(size=0.5, position=position_dodge(0.3))+
  labs(title="Patient C: \nMYOC", y="Normalized Counts")+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(color="black"),
        panel.border=element_rect(fill=NA),
        plot.title=element_text(hjust=0.5, size=6),
        axis.title=element_text(size=4),
        axis.text.y=element_text(size=4),
        axis.text.x=element_text(size=4),
        axis.title.x=element_blank(),
        legend.text=element_text(size=5),
        legend.title=element_text(size=6),
        legend.key.size=unit(1.5, 'mm'))
#pt 4, 6 top upreg gene plot==ENSCAFG00000047433
pt4plot <- plotCounts(dds, normalized=TRUE, gene="ENSCAFG00000047433", intgroup =c("Source", "Dog"), returnData=TRUE)
pt4plot <- pt4plot %>%
  ggplot(aes(x=Source,y=count,color=Dog)) +
  geom_line(aes(group=Dog), position=position_dodge(0.3), size=0.4)+
  geom_point(size=0.5, position=position_dodge(0.3))+
  labs(title="Patients D, F: \nMEPE", y="Normalized Counts")+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(color="black"),
        panel.border=element_rect(fill=NA),
        plot.title=element_text(hjust=0.5, size=6),
        axis.title=element_text(size=4),
        axis.text.y=element_text(size=4),
        axis.text.x=element_text(size=4),
        axis.title.x=element_blank(),
        legend.text=element_text(size=5),
        legend.title=element_text(size=6),
        legend.key.size=unit(1.5, 'mm'))
## Combine all the plots into one figure--FIGURE 6
comboplot3 <- ggarrange(pt1plot, B2, pt3plot, pt4plot, 
                        labels=c("A","B","C","D"), font.label=list(color="black", size=8),
                        ncol=2, nrow=2, widths=c(2,2), heights=c(2,2))
comboplot3 <- annotate_figure(comboplot3, top=text_grob("Top Downregulated Gene in Each Patient", face="bold", size=8))
ggsave("Figure6.png", comboplot3, width=180, height=100, units="mm", dpi=600)
ggsave("Figure6.pdf", comboplot3, width=180, height=100, units="mm", dpi=600)



##### ZOOM OF TOP UPREG GENE IN PATIENT A
A <- plotCounts(dds, normalized=TRUE, gene="ENSCAFG00000041995", intgroup =c("Source", "Dog"), returnData=TRUE)
A <- A %>%
  ggplot(aes(x=Source,y=count,color=Dog)) +
  geom_line(aes(group=Dog), position=position_dodge(0.3), size=1.2)+
  geom_point(size=1.4, position=position_dodge(0.3))+
  labs(title="Patient A: \nENSCAFG00000041995", y="Normalized Counts")+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(color="black"),
        panel.border=element_rect(fill=NA),
        plot.title=element_text(hjust=0.5, size=15),
        axis.title=element_text(size=11),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=11),
        axis.title.x=element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.key.size=unit(5, 'mm'))
A <- A + facet_zoom(ylim=c(0,2000))




