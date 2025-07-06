# Step 1: 读取数据
library(Matrix)
genes <- read.delim(gzfile("adultMouse/GSE126074_AdBrainCortex_SNAREseq_cDNA.genes.tsv.gz"), header = FALSE)
rna_counts <- readMM("adultMouse/GSE126074_AdBrainCortex_SNAREseq_cDNA.counts.mtx.gz")
barcodes <- read.delim(gzfile("adultMouse/GSE126074_AdBrainCortex_SNAREseq_cDNA.barcodes.tsv.gz"), header = FALSE)
head(barcodes)
head(genes)

dim(rna_counts)