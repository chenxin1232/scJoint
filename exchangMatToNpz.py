import numpy as np
from scipy.io import loadmat
from scipy import sparse

# 1. 加载表达矩阵
expr_mat = loadmat('data1/example1A_S.mat')
X = expr_mat['X'] 
X_sparse = sparse.csr_matrix(X)
print(X_sparse)
# 2. 加载标签
label_mat = loadmat('data1/example1A_S_label.mat')
# print(label_mat)
labels = label_mat['Cx_truth']  # 替换为你实际的变量名
labels = labels.squeeze()     # 转为一维数组

# 3. 保存为 npz 和 npy
sparse.save_npz('coupleCocData/example1A_S_label.npz',X_sparse)   # 保存表达矩阵
np.savetxt('coupleCocData/example1A_S.txt', labels,fmt='%s')  # 保存标签数组
