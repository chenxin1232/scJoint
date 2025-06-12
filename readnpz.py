import numpy as np

# 读取 .npz 文件
data = np.load('data/citeseq_control_rna.npz')
print(data)

# 查看有哪些键
print(data.files)  # ['X'] 或类似

# 假设里面有个数组 'X'，写出成 .txt 文件
np.savetxt('citeseq_control_rna.txt', data['data'], fmt='%.4f')