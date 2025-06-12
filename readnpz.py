import scipy.sparse

# 读取稀疏矩阵
sparse_matrix = scipy.sparse.load_npz('data_10x/exprs_10xPBMC_rna.npz')

# 打印基本信息
print("矩阵形状:", sparse_matrix.shape)
print("非零元素数量:", sparse_matrix.nnz)
print("稀疏矩阵格式:", type(sparse_matrix))

# 查看稀疏矩阵的前三行
print("\n前三行数据（转成密集数组显示）:")
print(sparse_matrix[:3, :].toarray())

# 查看部分非零元素的索引和值（示范）
rows, cols = sparse_matrix.nonzero()
values = sparse_matrix.data
print("\n前10个非零元素的位置和值:")
for i in range(min(10, len(values))):
    print(f"位置: ({rows[i]}, {cols[i]}), 值: {values[i]}")
