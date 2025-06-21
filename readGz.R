atac <- read.csv(gzfile("GSE149683_File_S6.Cicero_gene_activity_scores_by_cell_type.csv.gz"))
# rna <- read.csv(gzfile("GSE156793_S6_gene_expression_celltype.txt.gz"))
head(atac)
library(tidyr)

# 把 row/col/x 转成 gene activity matrix
atac_matrix <- pivot_wider(atac, names_from = col, values_from = x, values_fill = 0)

# 查看前几行
head(atac_matrix)