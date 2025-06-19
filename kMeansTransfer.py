import os
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.cluster import KMeans

# 加载嵌入向量
rna_embeddings = np.loadtxt('./output/rna_filtered_t_embeddings.txt')
atac_embeddings = np.loadtxt('./output/atac_filtered_t_embeddings.txt')

# 合并嵌入向量
embeddings =  np.concatenate((rna_embeddings, atac_embeddings))

# 加载标签
rna_labels = np.loadtxt('./HematopoiesisData/rna_labels.txt')
atac_labels = np.loadtxt('./HematopoiesisData/atac_labels.txt')

# 聚类操作
n_clusters = len(np.unique(rna_labels))  # 聚类数量为已知标签的类别数
kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(embeddings)
cluster_labels = kmeans.labels_

# 确定每个簇的标签
cluster_to_label = {}
for cluster in range(n_clusters):
    # 找到属于该簇的 RNA 数据的索引
    rna_cluster_indices = np.where(cluster_labels[:len(rna_embeddings)] == cluster)[0]
    if len(rna_cluster_indices) > 0:
        # 统计该簇中 RNA 数据的标签分布
        label_counts = np.bincount(rna_labels[rna_cluster_indices].astype(int))
        # 选择出现次数最多的标签作为该簇的标签
        cluster_to_label[cluster] = np.argmax(label_counts)

# 标签迁移
atac_predicted_labels = []
for cluster in cluster_labels[len(rna_embeddings):]:
    atac_predicted_labels.append(cluster_to_label[cluster])
atac_predicted_labels = np.array(atac_predicted_labels)

# 计算准确率
common_label_cnt = 0
correct = 0
for i in range(len(atac_labels)):
    if atac_labels[i] >= 0:
        common_label_cnt += 1
        if atac_labels[i] == atac_predicted_labels[i]:
            correct += 1
accuracy = correct / common_label_cnt
print(f"标签迁移准确率: {accuracy}")

# t-SNE 降维
tsne_results = TSNE(perplexity=30, n_iter = 1000).fit_transform(embeddings)

# 创建 DataFrame
df = pd.DataFrame()
df['tSNE1'] = tsne_results[:,0]
df['tSNE2'] = tsne_results[:,1]

data_label = np.array(["CITE-seq", "ASAP-seq"])
df['data'] = np.repeat(data_label, [rna_embeddings.shape[0], atac_embeddings.shape[0]], axis=0)

# 合并 RNA 标签和迁移后的 ATAC 标签
labels =  np.concatenate((rna_labels, atac_predicted_labels))
label_to_idx = pd.read_csv('./HematopoiesisData/celltype_label_mapping.txt', sep = '\t', header = None)
label_dic = []
for i in range(label_to_idx.shape[0]):
    label_dic = np.append(label_dic, label_to_idx[0][i][:-2])
df['predicted'] = label_dic[labels.astype(int)]

# 可视化
plt.figure(figsize=(10,10))
sns.scatterplot(
    x = "tSNE1", y = "tSNE2",
    hue = "data",
    palette = sns.color_palette("tab10", 2),
    data = df,
    legend = "full",
    alpha = 0.3
)
plt.show()