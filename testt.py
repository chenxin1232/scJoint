from scipy.io import mmread
from scipy.sparse import save_npz


# 读取 Matrix Market 格式的稀疏矩阵
rna = mmread("rna_filtered_transposed.mtx").tocsr()
atac = mmread("atac_filtered_transposed.mtx").tocsr()

# 保存为 .npz 文件（稀疏格式）
save_npz("rna_filtered_t.npz", rna)
save_npz("atac_filtered_t.npz", atac)
