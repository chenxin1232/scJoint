import pandas as pd
from scipy.sparse import load_npz

# RNA
rna_sparse = load_npz('data/citeseq_control_rna.npz')
rna_dense = rna_sparse.toarray()
pd.DataFrame(rna_dense).to_csv('citeseq_rna.csv', index=False)

# ATAC
atac_sparse = load_npz('data/asapseq_control_atac.npz')
atac_dense = atac_sparse.toarray()
pd.DataFrame(atac_dense).to_csv('asapseq_atac.csv', index=False)
# RNA 标签
rna_labels = pd.read_csv('data/citeseq_control_cellTypes.txt', header=None)
rna_labels.to_csv('citeseq_rna_labels.csv', index=False)

# ATAC 标签
atac_labels = pd.read_csv('data/asapseq_control_cellTypes.txt', header=None)
atac_labels.to_csv('asapseq_atac_labels.csv', index=False)