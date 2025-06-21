# ---- 0. 加载所需包 ----
library(Seurat)
library(Signac)
library(caret)
library(pheatmap)
library(ggplot2)
library(SummarizedExperiment)

# ---- 1. 读取 RNA 和 ATAC Seurat 对象 ----
rna <- readRDS("scRNA-Healthy-Hematopoiesis-191120 (1).rds")
atac <- readRDS("scATAC-Cicero-GA-Hematopoiesis-191120.rds")

# ---- 2. 提取表达矩阵和元数据 ----
# 确认 RNA counts 层
if ("counts" %in% names(assays(rna))) {
  rna_counts <- assays(rna)[["counts"]]
} else {
  rna_counts <- assays(rna)[[1]]
}
# 确认 ATAC gene activity 层
if (!"gA" %in% names(assays(atac))) {
  stop("请确认 ATAC 对象中 gene activity 层名称，这里假设为 'gA'")
}
atac_counts <- assays(atac)[["gA"]]

rna_meta  <- as.data.frame(colData(rna))
atac_meta <- as.data.frame(colData(atac))

# ---- 3. 构建 Seurat 对象 ----
rna_seurat  <- CreateSeuratObject(counts = rna_counts,  meta.data = rna_meta)
atac_seurat <- CreateSeuratObject(counts = atac_counts, meta.data = atac_meta)

# ---- 4. RNA 数据标准流程 ----
rna_seurat <- NormalizeData(rna_seurat)
rna_seurat <- FindVariableFeatures(rna_seurat)
rna_seurat <- ScaleData(rna_seurat)
rna_seurat <- RunPCA(rna_seurat)

# ---- 5. ATAC 数据标准流程 ----
DefaultAssay(atac_seurat) <- "RNA"  # 将 gene activity 当作 RNA assay
atac_seurat <- NormalizeData(atac_seurat)
atac_seurat <- FindVariableFeatures(atac_seurat)
atac_seurat <- ScaleData(atac_seurat)
atac_seurat <- RunPCA(atac_seurat)

# ---- 6. 查找锚点 (Transfer Anchors) ----
anchors <- Seurat::FindTransferAnchors(
  reference = rna_seurat,
  query     = atac_seurat,
  normalization.method = "LogNormalize",
  reference.reduction  = "pca",
  dims = 1:30
)

# ---- 7. 执行 TransferData 并接收带预测标签的 Seurat 对象 ----
predictions <- Seurat::TransferData(
  anchorset = anchors,
  refdata   = rna_seurat$BioClassification,  # 请确认此列存在
  query     = atac_seurat,
  weight.reduction = "pca",
  dims      = 1:30
)

# predictions 本身就是一个 Seurat 对象，包含 predicted.id 等元数据
atac_seurat <- predictions

# ---- 8. 从 ATAC 原始 metadata 提取真实标签 ----
# 这里假设真标签列名为 "celltype"，请根据实际列名修改
if (!"BioClassification" %in% colnames(atac_meta)) {
  stop("请确认 atac_meta 中真标签的列名，并替换下方 'celltype'")
}
atac_labels <- atac_meta$BioClassification
stopifnot(length(atac_labels) == ncol(atac_seurat))

# 添加到 Seurat 对象中
atac_seurat$true.label <- atac_labels

# ---- 9. 性能评估 ----
predicted <- factor(atac_seurat$predicted.id)
truth     <- factor(atac_seurat$true.label)

# 对齐因子水平
all_levels <- union(levels(predicted), levels(truth))
predicted  <- factor(predicted, levels = all_levels)
truth      <- factor(truth,     levels = all_levels)

conf_matrix <- confusionMatrix(predicted, truth)
print(conf_matrix)

cat("\n性能评估指标:\n")
cat("Accuracy :", round(conf_matrix$overall["Accuracy"] * 100, 2), "%\n")
cat("Recall   :", round(mean(conf_matrix$byClass[,"Sensitivity"]) * 100, 2), "%\n")
cat("Precision:", round(mean(conf_matrix$byClass[,"Precision"]) * 100, 2), "%\n")
cat("F1 Score :", round(mean(conf_matrix$byClass[,"F1"]) * 100, 2), "%\n")
