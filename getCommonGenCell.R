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

# 过滤掉含 "unk" 的细胞
rna_keep_idx <- !grepl("unk", rna_cell_types, ignore.case = TRUE)
atac_keep_idx <- !grepl("unk", atac_cell_types, ignore.case = TRUE)

rna_filtered <- rna_common[, rna_keep_idx]
atac_filtered <- atac_common[, atac_keep_idx]
dim(rna_filtered) 
dim(atac_filtered)
# sum(is.na(assay(rna_filtered)))
# any(is.nan(assay(rna_filtered)))
# 获取过滤后的细胞类型
# rna_types_filtered <- as.character(colData(rna_filtered)$BioClassification)
# atac_types_filtered <- as.character(colData(atac_filtered)$BioClassification)

# # 构建合并映射规则
# merge_map <- list(
#   CD4_combined     = c("0_CD4.N1", "21_CD4.N2", "22_CD4.M"),
#   CD8_combined     = c("23_CD8.EM", "24_CD8.CM"),
#   CD14_combined    = c("11_CD14.Mono.1", "12_CD14.Mono.2"),
#   GMP_combined     = c("07_GMP", "08_GMP.Neut"),
#   Erythro_Baso     = c("04_Early.Baso", "02_Early.Eryth")
# )

# # 反向映射表（从原始类型映射到合并名）
# reverse_map <- unlist(lapply(names(merge_map), function(x) setNames(rep(x, length(merge_map[[x]])), merge_map[[x]])))

# # 将其他未归类的类型保持原样
# all_types <- unique(c(rna_types_filtered, atac_types_filtered))
# unmapped <- setdiff(all_types, names(reverse_map))
# for (t in unmapped) reverse_map[[t]] <- t

# # 应用映射，替换成合并后的类别名
# rna_types_merged <- reverse_map[rna_types_filtered]
# atac_types_merged <- reverse_map[atac_types_filtered]

# # 构建统一标签（从0开始编号）
# unique_types <- sort(unique(c(rna_types_merged, atac_types_merged)))
# type_to_label <- setNames(seq_along(unique_types) - 1, unique_types)

# # 生成标签向量
# rna_labels <- type_to_label[rna_types_merged]
# atac_labels <- type_to_label[atac_types_merged]

# # 构建映射表
# type_label_map <- data.frame(
#   cell_type = unique_types,
#   label = seq_along(unique_types) - 1,
#   stringsAsFactors = FALSE
# )

# # 标签文件（只含数字）
# rna_label_only <- data.frame(label = rna_labels)
# atac_label_only <- data.frame(label = atac_labels)

# # 写文件
# # write.table(type_label_map, file = "celltype_label_mapping.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# # write.table(rna_label_only, file = "rna_labels.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# # write.table(atac_label_only, file = "atac_labels.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# # 提取表达矩阵并转置
# rna_mat_t <- t(assay(rna_filtered))
# atac_mat_t <- t(assay(atac_filtered))

# # 保存稀疏矩阵
# writeMM(rna_mat_t, file = "rna_filtered_transposed.mtx")
# writeMM(atac_mat_t, file = "atac_filtered_transposed.mtx")
