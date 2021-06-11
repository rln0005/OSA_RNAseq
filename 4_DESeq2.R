
#load required programs
library(DESeq2)
library(genefilter)
library(pheatmap)
library("RColorBrewer")
library(EnhancedVolcano)
library(AnnotationDbi)
library("org.Cf.eg.db")
library(biomaRt)

#########################  Prepare Data/Annotate Gene IDs  #########################
#import count data from pythonDE ballgown -- 30951 genes
countdata <- as.matrix(read.csv("gene_count_matrix_reorganized2.csv", row.names="gene_id"))
#import meta data file
coldata <-(read.table("pheno_reorganized2.txt", header=TRUE, row.names=1))
#check all sample IDs in colData are also in CountData and match their orders
all(rownames(coldata) %in% colnames(countdata))
all(rownames(coldata) == colnames(countdata))
###replace row names(ensembl) in countdata with the gene ID
#obtain ensembl annotation data
ensembl <- useEnsembl(biomart="ensembl", dataset="clfamiliaris_gene_ensembl")
#may have to try this a few times to get it to work -- weird connection issue with biomaRt
genelist <- getBM(attributes=c("ensembl_gene_id", "external_gene_name","description"),
             filters = "ensembl_gene_id",
             values=row.names(countdata),
             mart=ensembl)
#check if all are present & match
all(rownames(countdata) %in% genelist[,1])
#rename blank spaces as 'blank'
genelist$external_gene_name[genelist$external_gene_name %in% c("")] <- "blank"
#replace row names with the gene symbol
row.names(countdata) <- genelist[match(row.names(countdata), genelist[,1]),2]
#replace column names with dog/sample ID
colnames(countdata) <- coldata[,1]
#remove blank rows (ie, genes that didn't map to known/annotated gene IDs)
countdata <- countdata[rownames(countdata) != "blank",]
#filter out genes with reads less than 1
countdata <- countdata[ rowSums(countdata) > 1, ]


#########################  DESeq2 Analysis  #########################
#create DEseq2 dataset
dds <- DESeqDataSetFromMatrix(countData = countdata, colData=coldata,  design = ~Source)
#set factors for statistical analysis 
dds$condition <- factor(dds$Source, levels=c("Bone","Tumor"))
#differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
res
#reorder based on smallest adjusted p value
resOrdered <- res[order(res$padj),]
resOrdered
#basic summary of the results
summary(res)
#how many adjusted p values less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)
#using a p value cutoff of 0.05 instead of 0.1
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)
resOrdered05 <- res05[order(res05$padj),]
resOrdered05
###### write results to a file
write.csv(as.data.frame(resOrdered05), file="DGESeq_results_biomart.csv") 

#########################  Preparing Data for GSEA/Cytoscape  #########################
#### Make a new column with the ranking score ranked list for GSEA -- keep the sign of fold change (positive or negative)
res2 <-  within(resOrdered05, rank <- sign(log2FoldChange) * -log10(pvalue))
res2 
#subset the results so only rank column
DGErank = subset( res2, select = c(rank) )
DGErank
### Write the file as a tab delimited file with .rnk for import into GSEA
write.table(as.data.frame(DGErank), file="DGErank_SignValues_biomart.rnk" , sep="\t", row.names=TRUE, quote = FALSE)  


####  Make a file for the normalized expression DATA
nt <- normTransform(dds) 
#compare to original count data
head(assay(nt))
head(countdata)
#make it a new dataframe
nt <- assay(nt)
ntDF<-as.data.frame(nt)
### Write the file as a tab delimited file with .txt for import into GSEA
write.table(as.data.frame(ntDF), file="NormTransExpressionData.txt", sep="\t", row.names=TRUE, quote = FALSE) 


#########################  Making Figures for Paper  #########################
## Extract transformed values 
rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)
ntd <- normTransform(dds)
## Compare all the different transformations
head(assay(rld), 3)
head(assay(vsd),3)
head(assay(dds),3)
head(res05)

## PCA plot--Figure 2
plotPCA(vsd, intgroup=c("Source"))

## Volcano Plot--Figure 3 (fold change>2, pvalue<0.05)
EnhancedVolcano(resOrdered05,
                lab=rownames(resOrdered05),
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 2, 
                labSize=4.5,
                pointSize=1.3,
                col=c('black', 'black', 'black', 'red3'), colAlpha = 1)

## Heatmap of Significant Genes--Figure 4
label_col <- c("1-Bone","2-Bone","3-Bone","4-Bone","5-Bone","6-Bone","7-Bone","1-Tumor","2-Tumor","3-Tumor","4-Tumor","5-Tumor","6-Tumor","7-Tumor")
#res=DEseq2 initial results
resDF <- as.data.frame(res)
vsd <- vst(dds)
vsd <- assay(vsd)
vsDF <- as.data.frame(vsd)
#vsDF$Gene <- rownames(vsDF) 
sigGenes <- rownames(resDF[resDF$padj <= 0.05,])
vsDF2 <- vsDF[rownames(vsDF) %in% sigGenes,]
anno <-(read.table("heatmap_anno.txt", header=TRUE, row.names=1))
pheatmap(vsDF2,
         labels_col=label_col,
         show_rownames=FALSE,
         cluster_rows=TRUE,
         cluster_cols=TRUE,
         scale="row",
         annotation_col=anno)

#############################################################################################

### Extract only the genes with p<0.05
p05genes <- head(resOrdered05,3039)
### List of most downregulated genes (within p<0.05)
FCdown <- p05genes[order(p05genes$log2FoldChange),]
### List of most upregulated genes (within p<0.05)
FCup <- p05genes[order(-p05genes$log2FoldChange),]
###this file contains the significant DEGs (p<0.05) for population analysis sorted by highest log2FoldChange (this is NOT the raw counts--should I be using raw counts??)
write.csv(as.data.frame(FCup), file="SigGenes_sortbyFCup.csv")


########  individual analysis --- NEED TO WORK ON THIS ##########
vsd <- varianceStabilizingTransformation(dds)
vsd <- assay(vsd)
colnames(vsd) <- coldata[,1]
#vsDF <- as.data.frame(vsd)
#vsDF <- colnames(vsDF)

#log2foldchange for each individual patient---(tumor/bone) 
pt1 <- as.data.frame(vsd[,8]/vsd[,1])
colnames(pt1) <- "pt1"
pt2 <- as.data.frame(vsd[,9]/vsd[,2])
colnames(pt2) <- "pt2"
pt3 <- as.data.frame(vsd[,10]/vsd[,3])
colnames(pt3) <- "pt3"
pt4 <- as.data.frame(vsd[,11]/vsd[,4])
colnames(pt4) <- "pt4"
pt5 <- as.data.frame(vsd[,12]/vsd[,5])
colnames(pt5) <- "pt5"
pt6 <- as.data.frame(vsd[,13]/vsd[,6])
colnames(pt6) <- "pt6"
pt7 <- as.data.frame(vsd[,14]/vsd[,7])
colnames(pt7) <- "pt7"
all <- cbind(pt1,pt2,pt3,pt4,pt5,pt6,pt7)
### this file contains the normalized log2FC (tumor/bone) for each patient and gene
write.csv(as.data.frame(all), file="individual_analysis.csv")

#reorder based on largest FC difference
pt1order<-pt1[order(-pt1$pt1), , drop = FALSE]
pt2order<-pt2[order(-pt2$pt2), , drop = FALSE]
pt3order<-pt3[order(-pt3$pt3), , drop = FALSE]
pt4order<-pt4[order(-pt4$pt4), , drop = FALSE]
pt5order<-pt5[order(-pt5$pt5), , drop = FALSE]
pt6order<-pt6[order(-pt6$pt6), , drop = FALSE]
pt7order<-pt7[order(-pt7$pt7), , drop = FALSE]



####### TRYING TO MAKE INDIVIDUAL GRAPHS --- NEED TO WORK ON THIS ######
library(ggplot2)
topGene <- rownames(res05)[which.min(res05$padj)]
plotCounts(dds, gene=topGene, intgroup=c("Source"))
#Plot individual genes
plotCounts(dds, gene=which.min(res05$padj), intgroup="Source")
plotCounts(dds, gene=topGene, intgroup=c("Source","Dog")) 
data <-plotCounts(dds, gene=topGene, intgroup=c("Source","Dog"), returnData=TRUE)
ggplot(data, aes(x=Source, y=count, color=Dog)) + scale_y_log10() + geom_point(position=position_jitter(width=.15, height=0))
ggplot(data) + geom_point(position=position_jitter(width=.15, height=0))
ggplot(data)






################  Alternative method of annotation--not used for publication  ################
#### Using dbi annotation database for gene IDs
#annotate ensemble gene IDs to symbols
symbols <- mapIds(org.Cf.eg.db, keys = rownames(countdata2), column = "SYMBOL", keytype = "ENSEMBL")
#make it a dataframe
df <- as.data.frame(symbols)
#label blank rows as 'blank'
df[,1][df[,1] %in% c(NA)] <- "blank"
#replace row names with the gene symbol
row.names(countdata2) <- df[match(row.names(countdata2), rownames(df)),1]
#remove blank rows (ie, genes that didn't map to known/annotated gene IDs)
countdata2 <- countdata2[rownames(countdata2) != "blank",]
#filter out genes with reads less than 1
countdata2 <- countdata2[ rowSums(countdata2) > 1, ]
