# RNA 预处理
rna_seurat <- NormalizeData(rna_seurat)
rna_seurat <- FindVariableFeatures(rna_seurat)
rna_seurat <- ScaleData(rna_seurat)
rna_seurat <- RunPCA(rna_seurat)

# ATAC 预处理（默认使用 RNA 方式）
DefaultAssay(atac_seurat) <- "RNA"
atac_seurat <- NormalizeData(atac_seurat)
atac_seurat <- FindVariableFeatures(atac_seurat)
atac_seurat <- ScaleData(atac_seurat)
atac_seurat <- RunPCA(atac_seurat)

# 找到锚点
anchors <- FindTransferAnchors(reference = rna_seurat,
                               query = atac_seurat,
                               normalization.method = "LogNormalize",
                               reference.reduction = "pca",
                               dims = 1:30)

# 标签迁移
predictions <- TransferData(anchorset = anchors,
                            refdata = rna_seurat$BioClassification,
                            query = atac_seurat,
                            weight.reduction = "pca",
                            dims = 1:30)

# 添加预测标签和真实标签
atac_seurat <- AddMetaData(atac_seurat, metadata = predictions)
atac_seurat$true.label <- atac_seurat$BioClassification
library(caret)

predicted <- factor(atac_seurat$predicted.id)
truth     <- factor(atac_seurat$true.label)

# 对齐标签
all_levels <- union(levels(predicted), levels(truth))
predicted <- factor(predicted, levels = all_levels)
truth     <- factor(truth, levels = all_levels)

conf_matrix <- confusionMatrix(predicted, truth)
print(conf_matrix)

# 输出性能指标
cat("\n性能评估指标:\n")
cat("Accuracy :", round(conf_matrix$overall["Accuracy"] * 100, 2), "%\n")
cat("Recall   :", round(mean(conf_matrix$byClass[,"Sensitivity"], na.rm = TRUE) * 100, 2), "%\n")
cat("Precision:", round(mean(conf_matrix$byClass[,"Precision"],   na.rm = TRUE) * 100, 2), "%\n")
cat("F1 Score :", round(mean(conf_matrix$byClass[,"F1"],          na.rm = TRUE) * 100, 2), "%\n")
