library(tidyverse)  # 推荐
# Step 1: 读取数据
atac <- read.csv(gzfile("GSE149683_File_S6.Cicero_gene_activity_scores_by_cell_type.csv.gz"))
rna <- read.csv(gzfile("GSE156793_S6_gene_expression_celltype.txt.gz"))
# Step 2: 宽表转换，基因 × 细胞类型
dim(atac)
dim(rna)
# atac_matrix <- atac %>%
#   pivot_wider(names_from = col, values_from = x, values_fill = 0)
# colnames(atac)
# colnames(rna)
# # Step 3: 可选 - 将 row 设置为 rownames 并删除首列
# atac_mat <- atac_matrix %>% column_to_rownames("row")

# # 查看结果
# print(dim(atac_mat))
# head(atac_mat[, 1:5])  # 显示前5个细胞类型列
