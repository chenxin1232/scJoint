library(SummarizedExperiment)
library(GenomicRanges)
library(Matrix)
library(S4Vectors)

# 读数据
rna <- readRDS("scRNA-Healthy-Hematopoiesis-191120 (1).rds")
atac <- readRDS("scATAC-Cicero-GA-Hematopoiesis-191120.rds")

# 找共同基因
rna_genes <- rownames(rna)
atac_genes <- rownames(atac)
common_genes <- intersect(rna_genes, atac_genes)

# 只保留共同基因
rna_common <- rna[common_genes, ]
atac_common <- atac[common_genes, ]

# 细胞类型信息
rna_cell_types <- as.character(colData(rna)$BioClassification)
atac_cell_types <- as.character(colData(atac)$BioClassification)
print(rna_cell_types)

# 过滤细胞类型，去掉含 "25" 或 "26" 的细胞
rna_keep_idx <- !grepl("^unk$", rna_cell_types, ignore.case = TRUE)
atac_keep_idx <- !grepl("^unk$", atac_cell_types, ignore.case = TRUE)

rna_filtered <- rna_common[, rna_keep_idx]
atac_filtered <- atac_common[, atac_keep_idx]
dim(rna_filtered)
dim(atac_filtered)
# 过滤后的细胞类型
rna_types_filtered <- as.character(colData(rna_filtered)$BioClassification)
atac_types_filtered <- as.character(colData(atac_filtered)$BioClassification)

# 合并所有细胞类型并排序去重，生成标签映射
unique_types <- sort(unique(c(rna_types_filtered, atac_types_filtered)))
type_to_label <- setNames(seq_along(unique_types) - 1, unique_types)  # 标签从0开始

# 为每个细胞生成对应数字标签
rna_labels <- type_to_label[rna_types_filtered]
atac_labels <- type_to_label[atac_types_filtered]

# 创建细胞类型与标签映射表
type_label_map <- data.frame(
  cell_type = unique_types,
  label = seq_along(unique_types) - 1,
  stringsAsFactors = FALSE
)

# 创建只有标签的细胞标签表（无barcode）
rna_label_only <- data.frame(label = rna_labels)
atac_label_only <- data.frame(label = atac_labels)

# 写文件，注意路径可根据需要修改
# write.table(type_label_map, file = "celltype_label_mapping.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(rna_label_only, file = "rna_labels.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# write.table(atac_label_only, file = "atac_labels.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# rna_mat <- assay(rna_filtered)
# atac_mat <- assay(atac_filtered)

# 转置矩阵
# rna_mat_t <- t(rna_mat)
# atac_mat_t <- t(atac_mat)

# 保存为 Matrix Market 格式（稀疏矩阵）
# writeMM(rna_mat_t, file = "rna_filtered_transposed.mtx")
# writeMM(atac_mat_t, file = "atac_filtered_transposed.mtx")