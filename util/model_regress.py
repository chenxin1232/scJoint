import torch
import torch.nn as nn


class Net_encoder(nn.Module):
    def __init__(self, input_size):
        super(Net_encoder, self).__init__()
        self.input_size = input_size
        self.k = 64
        self.f = 64

        self.encoder = nn.Sequential(
            nn.Linear(self.input_size, 64)
        )

    def forward(self, data):
        print("输入data的shape:", data.shape)
        print("输入data的元素总数:", data.numel())
        print("期望每个样本的特征维度:", self.input_size)
        print("元素总数除以输入维度的商和余数:", data.numel() // self.input_size, data.numel() % self.input_size)
        data = data.float().view(-1, self.input_size)
        
        embedding = self.encoder(data)

        return embedding


class Net_cell(nn.Module):
    def __init__(self, num_of_class):
        super(Net_cell, self).__init__()
        self.cell = nn.Sequential(
            nn.Linear(64, num_of_class)
        )

    def forward(self, embedding):
        cell_prediction = self.cell(embedding)

        return cell_prediction
