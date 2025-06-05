import process_db
import h5py
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
import random
random.seed(1)
rna_h5_files = ["data_10x/exprs_10xPBMC_rna.h5"] 
rna_label_files = ["data_10x/cellType_10xPBMC_rna.csv"] # csv file

atac_h5_files = ["data_10x/exprs_10xPBMC_atac.h5"]
atac_label_files = []
process_db.data_parsing(rna_h5_files, atac_h5_files)
rna_label = pd.read_csv(rna_label_files[0], index_col = 0)
rna_label
print(rna_label.value_counts(sort = False))
process_db.label_parsing(rna_label_files, atac_label_files)
rna_embeddings = np.loadtxt('./output/exprs_10xPBMC_rna_embeddings.txt')
atac_embeddings = np.loadtxt('./output/exprs_10xPBMC_atac_embeddings.txt')
print(rna_embeddings.shape)
print(atac_embeddings.shape)
embeddings =  np.concatenate((rna_embeddings, atac_embeddings))
print(embeddings.shape)
tsne_results = TSNE(perplexity=30, n_iter = 1000).fit_transform(embeddings)
tsne_results.shape
df = pd.DataFrame()
df['tSNE1'] = tsne_results[:,0]
df['tSNE2'] = tsne_results[:,1]
rna_labels = np.loadtxt('./data_10x/cellType_10xPBMC_rna.txt')
atac_predictions = np.loadtxt('./output/exprs_10xPBMC_atac_knn_predictions.txt')
labels =  np.concatenate((rna_labels, atac_predictions))
label_to_idx = pd.read_csv('./data_10x/label_to_idx.txt', sep = '\t', header = None)
label_to_idx.shape
label_dic = []
for i in range(label_to_idx.shape[0]):
    label_dic = np.append(label_dic, label_to_idx[0][i][:-2])
    data_label = np.array(["scRNA-seq", "scATAC-seq"])
df['data'] = np.repeat(data_label, [rna_embeddings.shape[0], atac_embeddings.shape[0]], axis=0)
df['predicted'] = label_dic[labels.astype(int)]
plt.figure(figsize=(10,10))
sns.scatterplot(
    x = "tSNE1", y = "tSNE2",
    hue = "data",
    palette = sns.color_palette("tab10", 2),
    data = df,
    legend = "full",
    alpha = 0.3
)

plt.figure(figsize=(10,10))
sns.scatterplot(
    x = "tSNE1", y = "tSNE2",
    hue = "predicted",
    data = df,
    legend = "full",
    alpha = 0.3
)
plt.show()