## ===dependencies==============================================================
library("DESeq2")
library("tidyverse")
library("pheatmap")
library("RColorBrewer")
library("clusterProfiler")
library("DOSE")
library("org.Mm.eg.db") 
library("magrittr")
library("ggpubr")

## ====set working directory====================================================
dir <- "/mnt/nfs/data/Kian_Lab/Xinlong/data/RNAseq_Michela/"
setwd(dir)


## ====load count reads files===================================================
### set variables "cts" and from featureCount countMatrix output
filePath_batch1 <- file.path(dir, "1st_batch", "countReads", "by_featureCounts",
                             "countMatrix_modified.count")
filePath_batch2 <- file.path(dir, "2nd_batch", "countReads", "by_featureCounts",
                             "countMatrix_modified.count")


cts1 <- read_tsv(filePath_batch1)
cts2 <- read_tsv(filePath_batch2)

glimpse(cts1) # check data structure
glimpse(cts2) # check data structure


sampleNames1 <- sub("../../mapping/","", colnames(cts1))[-1]
sampleNames2 <- sub("../mapping/","", colnames(cts2))[-1]

# store sample names (e.g.BT****.sorted.bam)
colnames(cts1)[-1] <- sub(".sorted.bam","",sampleNames1) 
colnames(cts2)[-1] <- sub(".sorted.bam","",sampleNames2) 

cts <- cts1 %>% 
  left_join(cts2, by = "Geneid")

rownames <- cts$Geneid
cts <- as.data.frame(cts, row.names = rownames)[-1]
rownames(cts) <- rownames
cts <- cts[-c(37:42)] # remove iPSC34_KO and iPSC_WT

### prepare colData table for sample information
cell_type1 <- rep(c("MEF_from_naive",
                    c(paste("naive", 
                            c("d10Ch+", "d12Ch+GFP+","d12Ch+","d14Ch+GFP+","d14Ch+"),sep = "_")),
                    "iPSC_from_naive"),3)
batches1 <- rep(c("65E1", "65E4", "65E5"),each = 7)
colData1 <- data.frame(sampleNames = sampleNames1, cell_type = cell_type1, 
                       batches = batches1, row.names = colnames(cts1)[-1])

# cell_type2 <- c(paste("total", rep(c("d10Ch+GFP+", "d10Ch+","d12Ch+GFP+","d12Ch+","d8Ch+"),3),
#                       sep = "_"), rep(c("iPS_KO","iPS_WT"), each = 3))
# remove wt and KO
cell_type2 <- paste("total", 
                    rep(c("d10Ch+GFP+", "d10Ch+","d12Ch+GFP+","d12Ch+","d8Ch+"),3),
                    sep = "_")
batches2 <- rep(c("73E6", "73E7", "73E8"),each = 5)
colData2 <- data.frame(sampleNames = sampleNames2[1:15], cell_type = cell_type2, 
                       batches = batches2, row.names = colnames(cts2)[2:16])


colData <- colData1 %>% rbind(colData2)
#colData <- rbind(colData1, colData2)[-c(37:42),] # remove iPSC34_KO and iPSC_WT 
glimpse(colData)

### check if rownames of coldata is equal to colnames of cts. It is critical
### that the columns of the count matrix and the rows of the column data
### (information about samples) are in the same order
all(rownames(colData) == colnames(cts))
all(ncol(cts) == nrow(colData))

## ====load variables to DESeq2=================================================
dds_unbiased <- DESeqDataSetFromMatrix(countData = cts,
                                       colData = colData,
                                       design = ~ 1)
#dds_unbiased <- estimateSizeFactors(dds_unbiased)

## ====pre-filter===============================================================
dds_unbiased <- dds_unbiased[ rowSums(counts(dds_unbiased)) > 1, ]

## ====Extracting transformed values============================================
rld <- rlog(dds_unbiased, blind = TRUE)
#vsd <- varianceStabilizingTransformation(dds_unbiased,blind = TRUE)
#vsd.fast <- vst(dds_unbiased,blind = TRUE)
head(assay(rld), 3)

## ====plot PCA=================================================================
plotPCA(rld, intgroup = c("cell_type", "batches"))

# It is also possible to customize the PCA plot using the ggplot function.
pcaData <- plotPCA(rld, intgroup = c("cell_type", "batches"),returnDat = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = cell_type, shape = batches)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw()

# extraction of gene weights in PCA
pca <- prcomp(t(assay(rld)))
loading <- abs(pca$rotation)
val_per_gene <- sweep(loading, 2, colSums(loading), "/")
colSums(val_per_gene)
#write.csv(val_per_gene, file = "./unsupervisedPCAdata.csv")

# plot PCA with 5 populations used for WGBS and iPSC wt/KO
cell_type_5pops_plus_iPSC <- c("total_d10Ch+" ,"total_d10Ch+GFP+", "naive_d12Ch+",
                               "naive_d12Ch+GFP+", "naive_d14Ch+GFP+",
                               "iPS_KO", "iPS_WT")
colData_5pops_plus_iPSC <-
  subset(colData, cell_type %in% cell_type_5pops_plus_iPSC)
#colData_5pops_plus_iPSC$batches <-
#  factor(c(rep(1:3, each = 4), rep(1:3, each = 2), rep(1:3, 2)))
colData_5pops_plus_iPSC$cell_type <-
  colData_5pops_plus_iPSC$cell_type[, drop = TRUE]
rld_5pops_plus_iPSC <-
  rld[, colData$cell_type %in% cell_type_5pops_plus_iPSC]

pcaData <- plotPCA(rld_5pops_plus_iPSC,
                   intgroup = c("cell_type", "batches"), returnDat = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = cell_type, shape = batches)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  theme_bw() +
  scale_shape_manual(values = c(16, 17, 15, 3, 4, 8, 11))

## ====Hierarchical Heatmap of the sample-to-sample distances===================
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)

rownames(sampleDistMatrix) <- paste(rld$batches, rld$cell_type, sep = "_")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

## ====Heatmap of the count matrix==============================================
library("vsn")
library("pheatmap")
ntd_unbiased <- normTransform(dds_unbiased)
meanSdPlot(assay(ntd_unbiased))
select <- order(rowMeans(counts(dds_unbiased,normalized = TRUE)),
                decreasing = TRUE)[1:100]
df <- as.data.frame(colData(dds_unbiased)[,c("cell_type", "batches")])
pheatmap(assay(ntd_unbiased)[select,], cluster_rows = TRUE, 
         show_rownames = FALSE,cluster_cols = FALSE, annotation_col = df)

## ====differnetial expression analysis=========================================
### select subsets for further analysis
cell_type_selected <- c("MEF_from_naive", "total_d10Ch+" ,"total_d10Ch+GFP+", "naive_d12Ch+", "naive_d12Ch+GFP+", "naive_d14Ch+GFP+", "iPSC_from_naive")
colData_selected <- subset(colData,
                           cell_type %in% cell_type_selected)
### Be careful here! cell_type is factor viriable and it would not change itself
### Therefore you have to manually drop unused levels.
colData_selected$cell_type <- colData_selected$cell_type[, drop = TRUE]

cts_selected <- cts[colData$cell_type %in% cell_type_selected]

### check if rownames of coldata is equal to colnames of cts
all(rownames(colData_selected) == colnames(cts_selected))
all(ncol(cts_selected) == nrow(colData_selected))

### load data into dds_selected
dds_selected <- DESeqDataSetFromMatrix(countData = cts_selected,
                                       colData = colData_selected, design = ~ cell_type)
dds_selected$cell_type <- factor(dds_selected$cell_type,
                        levels = cell_type_selected)

dds_selected <- dds_selected[rowSums(counts(dds_selected)) > 1, ]
dds_selected <- DESeq(dds_selected)

## ----contrast DEGs bewteen each pair of groups--------------------------------
### total_d10Ch+ vs MEF
res_Td10ChvsMEF <- 
  results(dds_selected, 
          contrast = c("cell_type", "total_d10Ch+", "MEF_from_naive"),
          alpha = 0.05)
summary(res_Td10ChvsMEF)
plotMA(res_Td10ChvsMEF, ylim = c(-2,2))
sum(res_Td10ChvsMEF$padj < 0.05, na.rm = TRUE)

### total_d10Ch+GFP+ vs total_d10Ch+
res_Td10ChGFPvsTd10Ch <- 
  results(dds_selected, 
          contrast = c("cell_type", "total_d10Ch+GFP+", "total_d10Ch+"),
          alpha = 0.05)
summary(res_Td10ChGFPvsTd10Ch)
plotMA(res_Td10ChGFPvsTd10Ch, ylim = c(-2,2))
sum(res_Td10ChGFPvsTd10Ch$padj < 0.05, na.rm = TRUE)


### naive_d12Ch+ vs total_d10Ch+GFP+
res_Nd12ChvsTd10ChGFP <- 
  results(dds_selected, 
          contrast = c("cell_type","naive_d12Ch+", "total_d10Ch+GFP+"),
          alpha = 0.05)
summary(res_Nd12ChvsTd10ChGFP)
plotMA(res_Nd12ChvsTd10ChGFP, ylim = c(-2,2))
sum(res_Nd12ChvsTd10ChGFP$padj < 0.05, na.rm = TRUE)

### naive_d12Ch+GFP+ vs naive_d12Ch+
res_Nd12ChGFPvsNd12Ch <- 
  results(dds_selected, 
          contrast = c("cell_type","naive_d12Ch+GFP+","naive_d12Ch+"),
          alpha = 0.05)
summary(res_Nd12ChGFPvsNd12Ch)
plotMA(res_Nd12ChGFPvsNd12Ch, ylim = c(-2,2))
sum(res_Nd12ChGFPvsNd12Ch$padj < 0.05, na.rm = TRUE)

### naive_d14Ch+GFP+ vs naive_d12Ch+GFP+
res_Nd14ChGFPvsNd12ChGFP <- 
  results(dds_selected, 
          contrast = c("cell_type","naive_d14Ch+GFP+","naive_d12Ch+GFP+"),
          alpha = 0.05)
summary(res_Nd14ChGFPvsNd12ChGFP)
plotMA(res_Nd14ChGFPvsNd12ChGFP, ylim = c(-2,2))
sum(res_Nd14ChGFPvsNd12ChGFP$padj < 0.05, na.rm = TRUE)

### iPSC vs naive_d14Ch+GFP+
res_iPSCvsNd14ChGFP <- 
  results(dds_selected, 
          contrast = c("cell_type","iPSC_from_naive","naive_d12Ch+GFP+"),
          alpha = 0.05)
summary(res_iPSCvsNd14ChGFP)
plotMA(res_iPSCvsNd14ChGFP, ylim = c(-2,2))
sum(res_iPSCvsNd14ChGFP$padj < 0.05, na.rm = TRUE)

res_Td10ChvsMEF <- na.omit(as.data.frame(res_Td10ChvsMEF))
res_Td10ChGFPvsTd10Ch <- na.omit(as.data.frame(res_Td10ChGFPvsTd10Ch))
res_Nd12ChvsTd10ChGFP <- na.omit(as.data.frame(res_Nd12ChvsTd10ChGFP))
res_Nd12ChGFPvsNd12Ch <- na.omit(as.data.frame(res_Nd12ChGFPvsNd12Ch))
res_Nd14ChGFPvsNd12ChGFP <- na.omit(as.data.frame(res_Nd14ChGFPvsNd12ChGFP))
res_iPSCvsNd14ChGFP <- na.omit(as.data.frame(res_iPSCvsNd14ChGFP))

# write.csv(res_Td10ChvsMEF, file = "res_Td10ChvsMEF.csv")
# write.csv(res_Td10ChGFPvsTd10Ch, file = "res_Td10ChGFPvsTd10Ch.csv")
# write.csv(res_Nd12ChvsTd10ChGFP, file = "res_Nd12ChvsTd10ChGFP.csv")
# write.csv(res_Nd12ChGFPvsNd12Ch, file = "res_Nd12ChGFPvsNd12Ch.csv")
# write.csv(res_Nd14ChGFPvsNd12ChGFP, file = "res_Nd14ChGFPvsNd12ChGFP.csv")
# write.csv(res_iPSCvsNd14ChGFP, file = "res_iPSCvsNd14ChGFP.csv")


## ----GO analsyis for 5 comparisons--------------------------------------------
# filter the DEGs across 5 comparisons with FDR < 0.05 and log2FC > 1
DEGs_filtered_Td10ChvsMEF <- res_Td10ChvsMEF %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1)
DEGs_filtered_Td10ChGFPvsTd10Ch <- res_Td10ChGFPvsTd10Ch %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1)
DEGs_filtered_Nd12ChvsTd10ChGFP <- res_Nd12ChvsTd10ChGFP %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1)
DEGs_filtered_Nd12ChGFPvsNd12Ch <- res_Nd12ChGFPvsNd12Ch %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1)
DEGs_filtered_Nd14ChGFPvsNd12ChGFP <- res_Nd14ChGFPvsNd12ChGFP %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1)
DEGs_filtered_iPSCvsNd14ChGFP <- res_iPSCvsNd14ChGFP %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1)

## Read DEGs and convert to gene ID: Symbol -> Entrezid
eg_DEGs_filtered_Td10ChvsMEF <- 
  bitr(rownames(DEGs_filtered_Td10ChvsMEF), 
       fromType = "SYMBOL", 
       toType = "ENTREZID",
       OrgDb = "org.Mm.eg.db")

eg_DEGs_filtered_Td10ChGFPvsTd10Ch <- 
  bitr(rownames(DEGs_filtered_Td10ChGFPvsTd10Ch), 
       fromType = "SYMBOL", 
       toType = "ENTREZID",
       OrgDb = "org.Mm.eg.db")

eg_DEGs_filtered_Nd12ChvsTd10ChGFP <- 
  bitr(rownames(DEGs_filtered_Nd12ChvsTd10ChGFP), 
       fromType = "SYMBOL", 
       toType = "ENTREZID",
       OrgDb = "org.Mm.eg.db")

eg_DEGs_filtered_Nd12ChGFPvsNd12Ch <- 
  bitr(rownames(DEGs_filtered_Nd12ChGFPvsNd12Ch), 
       fromType = "SYMBOL", 
       toType = "ENTREZID",
       OrgDb = "org.Mm.eg.db")

eg_DEGs_filtered_Nd14ChGFPvsNd12ChGFP <- 
  bitr(rownames(DEGs_filtered_Nd14ChGFPvsNd12ChGFP), 
       fromType = "SYMBOL", 
       toType = "ENTREZID",
       OrgDb = "org.Mm.eg.db")

eg_DEGs_filtered_iPSCvsNd14ChGFP <- 
  bitr(rownames(DEGs_filtered_iPSCvsNd14ChGFP), 
       fromType = "SYMBOL", 
       toType = "ENTREZID",
       OrgDb = "org.Mm.eg.db")

# Initialize org.Mm.eg.db
ids <- keys(org.Mm.eg.db, 'ENTREZID')
id_GO <- select(org.Mm.eg.db,keys = ids,columns = c("ENTREZID","GO"))
id_GO <- subset(id_GO,!is.na(GO)) #remove genes without GO annotation
id_GO <- select(org.Mm.eg.db,keys = ids,columns = c("ENTREZID","GO")) %>% subset(!is.na(GO))

length(unique(id_GO$ENTREZID)) #remove duplicates

# enrich GO terms for each eg
ego_DEGs_filtered_Td10ChvsMEF <- 
  enrichGO(gene = eg_DEGs_filtered_Td10ChvsMEF$ENTREZID,
           keyType = "ENTREZID",
           OrgDb = "org.Mm.eg.db",
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.01,
           qvalueCutoff = 0.05,
           readable = TRUE)

ego_DEGs_filtered_Td10ChGFPvsTd10Ch <- 
  enrichGO(gene = eg_DEGs_filtered_Td10ChGFPvsTd10Ch$ENTREZID,
           keyType = "ENTREZID",
           OrgDb = "org.Mm.eg.db",
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.01,
           qvalueCutoff = 0.05,
           readable = TRUE)

ego_DEGs_filtered_Nd12ChvsTd10ChGFP <- 
  enrichGO(gene = eg_DEGs_filtered_Nd12ChvsTd10ChGFP$ENTREZID,
           keyType = "ENTREZID",
           OrgDb = "org.Mm.eg.db",
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.01,
           qvalueCutoff = 0.05,
           readable = TRUE)

ego_DEGs_filtered_Nd12ChGFPvsNd12Ch <- 
  enrichGO(gene = eg_DEGs_filtered_Nd12ChGFPvsNd12Ch$ENTREZID,
           keyType = "ENTREZID",
           OrgDb = "org.Mm.eg.db",
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.01,
           qvalueCutoff = 0.05,
           readable = TRUE)

ego_DEGs_filtered_Nd14ChGFPvsNd12ChGFP <- 
  enrichGO(gene = eg_DEGs_filtered_Nd14ChGFPvsNd12ChGFP$ENTREZID,
           keyType = "ENTREZID",
           OrgDb = "org.Mm.eg.db",
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.01,
           qvalueCutoff = 0.05,
           readable = TRUE)

ego_DEGs_filtered_iPSCvsNd14ChGFP <- 
  enrichGO(gene = eg_DEGs_filtered_iPSCvsNd14ChGFP$ENTREZID,
           keyType = "ENTREZID",
           OrgDb = "org.Mm.eg.db",
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.01,
           qvalueCutoff = 0.05,
           readable = TRUE)

ego_DEGs_filtered_Td10ChvsMEF_df <- 
  as.data.frame(ego_DEGs_filtered_Td10ChvsMEF)  
ego_DEGs_filtered_Td10ChGFPvsTd10Ch_df <- 
  as.data.frame(ego_DEGs_filtered_Td10ChGFPvsTd10Ch)  # ego_df <- ego@result
ego_DEGs_filtered_Nd12ChvsTd10ChGFP_df <- 
  as.data.frame(ego_DEGs_filtered_Nd12ChvsTd10ChGFP)
ego_DEGs_filtered_Nd12ChGFPvsNd12Ch_df <- 
  as.data.frame(ego_DEGs_filtered_Nd12ChGFPvsNd12Ch)
ego_DEGs_filtered_Nd14ChGFPvsNd12ChGFP_df <- 
  as.data.frame(ego_DEGs_filtered_Nd14ChGFPvsNd12ChGFP)
ego_DEGs_filtered_iPSCvsNd14ChGFP_df <- 
  as.data.frame(ego_DEGs_filtered_iPSCvsNd14ChGFP)

dotplot(ego_DEGs_filtered_Td10ChvsMEF, showCategory = 30)
dotplot(ego_DEGs_filtered_Td10ChGFPvsTd10Ch, showCategory = 30)
dotplot(ego_DEGs_filtered_Nd12ChvsTd10ChGFP, showCategory = 30)
dotplot(ego_DEGs_filtered_Nd12ChGFPvsNd12Ch, showCategory = 30)
dotplot(ego_DEGs_filtered_Nd14ChGFPvsNd12ChGFP, showCategory = 30)
dotplot(ego_DEGs_filtered_iPSCvsNd14ChGFP, showCategory = 30)

write.csv(ego_DEGs_filtered_Td10ChvsMEF_df, 
          "./GO_DEGs_filtered_Td10ChvsMEF.csv")
write.csv(ego_DEGs_filtered_Td10ChGFPvsTd10Ch_df, 
          "./GO_DEGs_filtered_Td10ChGFPvsTd10Ch.csv")
write.csv(ego_DEGs_filtered_Td10ChGFPvsTd10Ch_df, 
          "./GO_DEGs_filtered_Td10ChGFPvsTd10Ch.csv")
write.csv(ego_DEGs_filtered_Nd12ChvsTd10ChGFP_df, 
          "./GO_DEGs_filtered_Nd12ChvsTd10ChGFP.csv")
write.csv(ego_DEGs_filtered_Nd12ChGFPvsNd12Ch_df, 
          "./GO_DEGs_filtered_Nd12ChGFPvsNd12Ch.csv")
write.csv(ego_DEGs_filtered_Nd14ChGFPvsNd12ChGFP_df, 
          "./GO_DEGs_filtered_Nd14ChGFPvsNd12ChGFP.csv")
write.csv(ego_DEGs_filtered_iPSCvsNd14ChGFP_df, 
          "./GO_DEGs_filtered_iPSCvsNd14ChGFP.csv")



## ----make customized heatmap with DEGs ---------------------------------------
# filter the DEGs across 5 comparisons with FDR < 0.05 and log2FC > 2
DEGs_filtered_Td10ChGFPvsTd10Ch <- res_Td10ChGFPvsTd10Ch %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 2)
DEGs_filtered_Nd12ChvsTd10ChGFP <- res_Nd12ChvsTd10ChGFP %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 2)
DEGs_filtered_Nd12ChGFPvsNd12Ch <- res_Nd12ChGFPvsNd12Ch %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 2)
DEGs_filtered_Nd14ChGFPvsNd12ChGFP <- res_Nd14ChGFPvsNd12ChGFP %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 2)
DEGs_filtered_iPSCvsNd14ChGFP <- res_iPSCvsNd14ChGFP %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 2)
  
# extract and combine the gene names of the filtered DEGs
DEGs_list <- list(DEGs_filtered_Td10ChGFPvsTd10Ch,
                  DEGs_filtered_Nd12ChvsTd10ChGFP,
                  DEGs_filtered_Nd12ChGFPvsNd12Ch,
                  DEGs_filtered_Nd14ChGFPvsNd12ChGFP,
                  DEGs_filtered_iPSCvsNd14ChGFP)

geneList_for_heatmap <- reduce(sapply(DEGs_list, rownames), union)
# return the unioned gene list without duplicates

# extract normalized count values and make matrix of above unioned gene list
# generate transformed count values scaled by log2 and library size.
# Note:this step usually has been done before but here we have to redo it 
# because the dataset is subsetted
rld_selected <- rlog(dds_selected, blind = FALSE)

idx <- match(geneList_for_heatmap, rownames(assay(rld_selected)))
filteredVarGeneMatrix <- assay(rld_selected)[idx, ]

# reorder dataset to group triplicate cell populations together
filteredVarGeneMatrix <- filteredVarGeneMatrix[, c(14, 16, 18, 13, 15, 17, 2, 6,
                                                   10, 1, 5, 9, 3, 7, 11, 4, 8, 12)]

annotation_data <- as.data.frame(colData(rld_selected))[, "cell_type", drop = FALSE]

heatmap  <- pheatmap(filteredVarGeneMatrix, 
                     annotation_col = annotation_data, # sample annotation
                     fontsize = 8,                     # general font size 
                     show_rownames = FALSE,            # hide gene names 
                     scale = 'row',                    # scale data by rows
                     cluster_rows = TRUE,              # culster data by rows
                     cutree_rows = 5,                  # split rows to 5 clusters
                     treeheight_row = 300,
                     cluster_cols = FALSE)

## ====cluster DEGs according to above heatmap===================================
# extract genes from each cluster (FC > 2 and cutree = 5)
genes_hclusted <- heatmap$tree_row
genes_hclusted_cutree5 <- cutree(genes_hclusted, 5)

# Cluster order here is different from the top-bottom cluster order in heatmap!!!
genes_clust1 <- which(genes_hclusted_cutree5 == 4) 
genes_clust2 <- which(genes_hclusted_cutree5 == 5) 
genes_clust3 <- which(genes_hclusted_cutree5 == 3) 
genes_clust4 <- which(genes_hclusted_cutree5 == 2) 
genes_clust5 <- which(genes_hclusted_cutree5 == 1)

cluster_order <- heatmap$tree_row$order
genes_cluster1_ordered <- genes_clust1[cluster_order,]
genes_cluster2_ordered <- genes_clust2[cluster_order,]
genes_cluster3_ordered <- genes_clust3[cluster_order,]
genes_cluster4_ordered <- genes_clust4[cluster_order,]
genes_cluster5_ordered <- genes_clust5[cluster_order,]
# filteredVarGeneMatrix.order <- 
#   cbind(filteredVarGeneMatrix[genes_hclusted[["order"]],],
#         cluster = cutree(genes_hclusted, k = 5)[genes_hclusted[["order"]]])

# export gene list of each clutser
write.table(names(genes_cluster1_ordered), 
            "./genes_cluster1_ordered.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(names(genes_cluster2_ordered ), 
            "./genes_cluster2_ordered.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(names(genes_cluster3_ordered), 
            "./genes_cluster3_ordered.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(names(genes_cluster4_ordered), 
            "./genes_cluster4_ordered.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(names(genes_cluster5_ordered), 
            "./genes_cluster1_ordered.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

## ----GO analysis for genes of each cluster------------------------------------
## Read DEGs and convert to gene ID: Symbol -> Entrezid
eg_culster1 <- bitr(names(genes_clust1), 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID",
                    OrgDb = "org.Mm.eg.db")

eg_culster2 <- bitr(names(genes_clust2), 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID",
                    OrgDb = "org.Mm.eg.db")

eg_culster3 <- bitr(names(genes_clust3), 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID",
                    OrgDb = "org.Mm.eg.db")

eg_culster4 <- bitr(names(genes_clust4), 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID",
                    OrgDb = "org.Mm.eg.db")

eg_culster5 <- bitr(names(genes_clust5), 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID",
                    OrgDb = "org.Mm.eg.db")

# Initialize org.Mm.eg.db
ids <- keys(org.Mm.eg.db, 'ENTREZID')
id_GO <- select(org.Mm.eg.db,keys = ids,columns = c("ENTREZID","GO"))
id_GO <- subset(id_GO,!is.na(GO)) #remove genes without GO annotation
id_GO <- select(org.Mm.eg.db,keys = ids,columns = c("ENTREZID","GO")) %>% subset(!is.na(GO))

length(unique(id_GO$ENTREZID)) #remove duplicates

# enrich GO terms for each eg
ego_cluster1 <- enrichGO(gene = eg_culster1$ENTREZID,
                         keyType = "ENTREZID",
                         OrgDb = "org.Mm.eg.db",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.01,
                         qvalueCutoff = 0.05,
                         readable = TRUE)

ego_cluster2 <- enrichGO(gene = eg_culster2$ENTREZID,
                         keyType = "ENTREZID",
                         OrgDb = "org.Mm.eg.db",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.01,
                         qvalueCutoff = 0.05,
                         readable = TRUE)

ego_cluster3 <- enrichGO(gene = eg_culster3$ENTREZID,
                         keyType = "ENTREZID",
                         OrgDb = "org.Mm.eg.db",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.01,
                         qvalueCutoff = 0.05,
                         readable = TRUE)

ego_cluster4 <- enrichGO(gene = eg_culster4$ENTREZID,
                         keyType = "ENTREZID",
                         OrgDb = "org.Mm.eg.db",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.01,
                         qvalueCutoff = 0.05,
                         readable = TRUE)

ego_cluster5 <- enrichGO(gene = eg_culster5$ENTREZID,
                         keyType = "ENTREZID",
                         OrgDb = "org.Mm.eg.db",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.01,
                         qvalueCutoff = 0.05,
                         readable = TRUE)

ego_cluster45 <- enrichGO(gene = c(eg_culster4$ENTREZID, eg_culster5$ENTREZID),
                          keyType = "ENTREZID",
                          OrgDb = "org.Mm.eg.db",
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.05,
                          readable = TRUE)

ego_cluster1_df <- as.data.frame(ego_cluster1)  # ego_df <- ego@result
ego_cluster2_df <- as.data.frame(ego_cluster2)
ego_cluster3_df <- as.data.frame(ego_cluster3)
ego_cluster4_df <- as.data.frame(ego_cluster4)
ego_cluster5_df <- as.data.frame(ego_cluster5)
ego_cluster45_df <- as.data.frame(ego_cluster45)

dotplot(ego_cluster1, showCategory = 30)
dotplot(ego_cluster2, showCategory = 30)
dotplot(ego_cluster3, showCategory = 30)
dotplot(ego_cluster4, showCategory = 30)
dotplot(ego_cluster5, showCategory = 30)
dotplot(ego_cluster45, showCategory = 30)

write.csv(ego_cluster1, "./GO_cluster1.csv")
write.csv(ego_cluster2, "./GO_cluster2.csv")
write.csv(ego_cluster3, "./GO_cluster3.csv")
write.csv(ego_cluster4, "./GO_cluster4.csv")
write.csv(ego_cluster5, "./GO_cluster5.csv")
write.csv(ego_cluster45, "./GO_cluster45.csv")
