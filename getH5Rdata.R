# 加载所需库
library(Signac)
library(Seurat)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # UCSC 风格注释
library(SummarizedExperiment)

# 设置随机种子
set.seed(1234)

# Step 1: 读取 .rds 文件（SummarizedExperiment 格式）
rse <- readRDS("scATAC-Healthy-Hematopoiesis-191120.rds")

# Step 2: 提取 counts 和 peaks (GRanges)
peak_counts <- assay(rse, "counts")
peak_ranges <- rowRanges(rse)

# Step 3: 设置 genome 为 'hg38'（必须）
genome(peak_ranges) <- "hg38"

# Step 4: 构建 ChromatinAssay（此时不加 annotation，避免冲突）
chrom_assay <- CreateChromatinAssay(
  counts = peak_counts,
  ranges = peak_ranges,
  genome = "hg38"
)

# Step 5: 构建 Seurat 对象
seurat_atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = as.data.frame(colData(rse))
)

# Step 6: 添加基因注释（UCSC 风格，避免染色体命名冲突）
annotation <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
genome(annotation) <- "hg38"
Annotation(seurat_atac) <- annotation

# Step 7: 添加 fragment 文件路径（这是 GeneActivity 所必须的！）
# 请确保 fragments 文件真实存在，路径正确
seurat_atac <- SetFragments(
  object = seurat_atac,
  file = "fragments.tsv.gz"  # 替换成你的 fragments 文件路径
)

# Step 8: 设置默认 assay
DefaultAssay(seurat_atac) <- "peaks"

# Step 9: 生成基因活性矩阵（延伸 2kb 到启动子区）
gene_activities <- GeneActivity(
  object = seurat_atac,
  extend.upstream = 2000,
  extend.downstream = 0
)

# Step 10: 加入 gene activity 为新的 assay
seurat_atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene_activities)

# Step 11: 设置默认 assay 为 ACTIVITY 并归一化
DefaultAssay(seurat_atac) <- "ACTIVITY"
seurat_atac <- NormalizeData(seurat_atac)
seurat_atac <- FindVariableFeatures(seurat_atac)
seurat_atac <- ScaleData(seurat_atac)

# Step 12: 保存 gene activity 矩阵
write.csv(as.matrix(seurat_atac[["ACTIVITY"]]@counts), file = "atac_gene_activity_matrix.csv")
cat("✅ 基因活性矩阵已保存为 atac_gene_activity_matrix.csv\n")
