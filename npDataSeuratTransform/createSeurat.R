library(Seurat)
library(Signac)
library(caret)

# 1. 读入表达矩阵
df <- read.csv("citeseq_rna.csv", header = TRUE, stringsAsFactors = FALSE)
# atac_counts <- read.csv("asapseq_atac.csv", row.names = 1)
head(df[,1], 20)
any(duplicated(df[,1]))  # 是否有重复
# 转换为矩阵格式
rna_counts <- as.matrix(rna_counts)
atac_counts <- as.matrix(atac_counts)

# 2. 读入细胞标签（这里假设每列对应一个细胞，每行一个标签）
rna_labels <- read.csv("citeseq_rna_labels.csv", header = FALSE)[[1]]
atac_labels <- read.csv("asapseq_atac_labels.csv", header = FALSE)[[1]]

# 确保细胞数一致
stopifnot(length(rna_labels) == ncol(rna_counts))
stopifnot(length(atac_labels) == ncol(atac_counts))

# 3. 构建 Seurat 对象（注意设置 meta.data）
rna_seurat <- CreateSeuratObject(counts = rna_counts,
                                  meta.data = data.frame(BioClassification = rna_labels))

atac_seurat <- CreateSeuratObject(counts = atac_counts,
                                   meta.data = data.frame(BioClassification = atac_labels))
