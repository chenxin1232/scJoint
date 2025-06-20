# ---- 0. 加载所需包 ----
library(Seurat)
library(Signac)
library(caret)
library(pheatmap)
library(ggplot2)
library(SummarizedExperiment)
# ---- 1. 加载 RNA 和 ATAC Seurat 对象 ----
rna <- readRDS("scRNA-Healthy-Hematopoiesis-191120 (1).rds")
atac <- readRDS("scATAC-Cicero-GA-Hematopoiesis-191120.rds")
assayNames(atac)

# 查看每个 assay 的结构
lapply(assays(atac), class)

# 查看具体维度（是否是 matrix 或稀疏矩阵）
lapply(assays(atac), dim)
# 提取表达矩阵
rna_counts <- assays(rna)$counts  # 如果 assays(rna) 中不是 counts，请改成实际名称
atac_counts <- assays(atac)[["gA"]]
fatac_counts <- assays(atac)$counts
# 提取元数据
rna_meta <- as.data.frame(colData(rna))
atac_meta <- as.data.frame(colData(atac))
rna_seurat <- CreateSeuratObject(counts = rna_counts, meta.data = rna_meta)
atac_seurat <- CreateSeuratObject(counts = atac_counts, meta.data = atac_meta)
rna_seurat <- NormalizeData(rna_seurat)
rna_seurat <- FindVariableFeatures(rna_seurat)
rna_seurat <- ScaleData(rna_seurat)
rna_seurat <- RunPCA(rna_seurat)

# Step 2: ATAC 预处理
DefaultAssay(atac_seurat) <- "RNA"
atac_seurat <- NormalizeData(atac_seurat)
atac_seurat <- ScaleData(atac_seurat)

# Step 3: 标签迁移
anchors <- FindTransferAnchors(
  reference = rna_seurat,
  query = atac_seurat,
  reduction = "cca",
  dims = 1:30
)

# 你可以替换 celltype 为你的真实标签列名
predictions <- TransferData(
  anchorset = anchors,
  refdata = rna_seurat$BioClassification,
  dims = 1:30,
)

atac_seurat <- AddMetaData(atac_seurat, metadata = predictions)

# Step 4: 添加真实标签
# 假设你已经有 atac_labels 向量（顺序对应于 atac_seurat）
atac_seurat$true.label <- atac_labels

# Step 5: 准确率等评估
library(caret)

predicted <- factor(atac_seurat$predicted.id)
truth <- factor(atac_seurat$true.label)

all_levels <- union(levels(predicted), levels(truth))
predicted <- factor(predicted, levels = all_levels)
truth <- factor(truth, levels = all_levels)

conf_matrix <- confusionMatrix(predicted, truth)
print(conf_matrix)

cat("\n性能评估:\n")
cat("Accuracy:", round(conf_matrix$overall["Accuracy"] * 100, 2), "%\n")
cat("Recall:", round(mean(conf_matrix$byClass[,"Sensitivity"], na.rm = TRUE) * 100, 2), "%\n")
cat("Precision:", round(mean(conf_matrix$byClass[,"Precision"], na.rm = TRUE) * 100, 2), "%\n")
cat("F1 Score:", round(mean(conf_matrix$byClass[,"F1"], na.rm = TRUE) * 100, 2), "%\n")