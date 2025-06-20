import os
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
rna_embeddings = np.loadtxt('./output/rna_filtered_t_embeddings.txt')
atac_embeddings = np.loadtxt('./output/atac_filtered_t_embeddings.txt')
print(rna_embeddings.shape)
print(atac_embeddings.shape)
embeddings =  np.concatenate((rna_embeddings, atac_embeddings))
print(embeddings.shape)
tsne_results = TSNE(perplexity=30, n_iter = 1000).fit_transform(embeddings)
tsne_results.shape
df = pd.DataFrame()
df['tSNE1'] = tsne_results[:,0]
df['tSNE2'] = tsne_results[:,1]
rna_labels = np.loadtxt('./HematopoiesisData/rna_labels.txt')
atac_labels = np.loadtxt('./HematopoiesisData/atac_labels.txt')
atac_predictions = np.loadtxt('./output/atac_filtered_t_knn_predictions.txt')
labels =  np.concatenate((rna_labels, atac_predictions))
label_to_idx = pd.read_csv('./HematopoiesisData/celltype_label_mapping.txt', sep = '\t', header = None)
label_to_idx.shape
label_dic = []
for i in range(label_to_idx.shape[0]):
    label_dic = np.append(label_dic, label_to_idx[0][i][:-2])
    common_label_cnt = 0
correct = 0
for i in range(len(atac_labels)):
    if atac_labels[i] >= 0:
        common_label_cnt += 1
        if atac_labels[i] == atac_predictions[i]:
            correct += 1
correct/common_label_cnt
print(correct/common_label_cnt)
data_label = np.array(["CITE-seq", "ASAP-seq"])
df['data'] = np.repeat(data_label, [rna_embeddings.shape[0], atac_embeddings.shape[0]], axis=0)
df['predicted'] = label_dic[labels.astype(int)]
plt.figure(figsize=(10,10))
print(df['data'])
sns.scatterplot(
    x = "tSNE1", y = "tSNE2",
    hue = "data",
    palette = sns.color_palette("tab10", 10),
    data = df,
    legend = "full",
    alpha = 0.3
)
plt.figure(figsize=(10,10))
sns.scatterplot(
    x = "tSNE1", y = "tSNE2",
    hue = "predicted",
    palette = sns.color_palette("hls", 10),
    data = df,
    legend = "full",
    alpha = 0.3
)
plt.show()