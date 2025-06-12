import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.neighbors import NearestNeighbors
import numpy as np

# 加载共享空间的嵌入
rna_embeddings = np.loadtxt('./output/citeseq_control_rna_embeddings.txt')
atac_embeddings = np.loadtxt('./output/asapseq_control_atac_embeddings.txt')

rna_embeddings = torch.tensor(rna_embeddings, dtype=torch.float32)
atac_embeddings = torch.tensor(atac_embeddings, dtype=torch.float32)

# 用 k-NN 找每个 scRNA 对应最近的 scATAC
k = 10
nn_model = NearestNeighbors(n_neighbors=k)
nn_model.fit(atac_embeddings.numpy())
distances, indices = nn_model.kneighbors(rna_embeddings.numpy())  # [n_rna, k]

# -------------------------------
# Top-K 最近邻加权平均构造伪配对
# -------------------------------
# 权重为 1 / distance
weights = 1.0 / (distances + 1e-8)  # 防止除0
weights = weights / weights.sum(axis=1, keepdims=True)  # 每行归一化为 1
weights = torch.tensor(weights, dtype=torch.float32)  # shape [n_rna, k]

# 构造 weighted paired atac
paired_atac = torch.stack([
    (atac_embeddings[neighbor_indices] * weight[:, None]).sum(dim=0)
    for neighbor_indices, weight in zip(indices, weights)
])
# shape: [n_rna, embedding_dim]

# -------------------------------
# 定义 GAN 网络
# -------------------------------
class Generator(nn.Module):
    def __init__(self, input_dim):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, 256),
            nn.ReLU(),
            nn.Linear(256, 256),
            nn.ReLU(),
            nn.Linear(256, input_dim)
        )

    def forward(self, x):
        return self.net(x)

class Discriminator(nn.Module):
    def __init__(self, input_dim):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, 256),
            nn.LeakyReLU(0.2),
            nn.Linear(256, 128),
            nn.LeakyReLU(0.2),
            nn.Linear(128, 1),
            nn.Sigmoid()
        )

    def forward(self, x):
        return self.net(x)

input_dim = rna_embeddings.shape[1]
G = Generator(input_dim)
D = Discriminator(input_dim)

optimizer_G = optim.Adam(G.parameters(), lr=0.0001)
optimizer_D = optim.Adam(D.parameters(), lr=0.0001)
criterion = nn.BCELoss()

# -------------------------------
# 训练 GAN
# -------------------------------
epochs = 100
batch_size = 64
num_batches = rna_embeddings.shape[0] // batch_size

for epoch in range(epochs):
    for i in range(num_batches):
        idx = slice(i * batch_size, (i + 1) * batch_size)
        rna_batch = rna_embeddings[idx]
        real_atac_batch = paired_atac[idx]

        # 训练判别器
        optimizer_D.zero_grad()
        real_labels = torch.ones((batch_size, 1))
        fake_labels = torch.zeros((batch_size, 1))

        real_output = D(real_atac_batch)
        d_real_loss = criterion(real_output, real_labels)

        fake_atac = G(rna_batch)
        fake_output = D(fake_atac.detach())
        d_fake_loss = criterion(fake_output, fake_labels)

        d_loss = d_real_loss + d_fake_loss
        d_loss.backward()
        optimizer_D.step()

        # 训练生成器
        optimizer_G.zero_grad()
        fake_output = D(fake_atac)
        g_loss = criterion(fake_output, real_labels)
        g_loss.backward()
        optimizer_G.step()

    print(f"[Epoch {epoch+1}/{epochs}] D_loss: {d_loss.item():.4f} | G_loss: {g_loss.item():.4f}")
