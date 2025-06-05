import torch
import torch.nn as nn
import torch.nn.functional as F

class ResidualBlock(nn.Module):
    """残差块实现"""
    def __init__(self, in_dim, out_dim, dropout=0.2, activation='relu'):
        super(ResidualBlock, self).__init__()
        
        # 主路径
        self.layer1 = nn.Linear(in_dim, out_dim)
        self.bn1 = nn.BatchNorm1d(out_dim)
        self.dropout = nn.Dropout(dropout)
        
        self.layer2 = nn.Linear(out_dim, out_dim)
        self.bn2 = nn.BatchNorm1d(out_dim)
        
        # 跳跃连接 (如果输入输出维度不同，需要进行投影)
        self.need_projection = (in_dim != out_dim)
        if self.need_projection:
            self.projection = nn.Sequential(
                nn.Linear(in_dim, out_dim),
                nn.BatchNorm1d(out_dim)
            )
        
        # 激活函数
        if activation == 'relu':
            self.activation = nn.ReLU()
        elif activation == 'leaky_relu':
            self.activation = nn.LeakyReLU(0.2)
        else:
            self.activation = nn.ReLU()
    
    def forward(self, x):
        identity = x
        
        # 主路径
        out = self.layer1(x)
        out = self.bn1(out)
        out = self.activation(out)
        out = self.dropout(out)
        
        out = self.layer2(out)
        out = self.bn2(out)
        
        # 跳跃连接
        if self.need_projection:
            identity = self.projection(identity)
            
        out += identity  # 残差连接的核心
        out = self.activation(out)
        
        return out

class Encoder(nn.Module):
    """使用残差块的编码器"""
    def __init__(self, input_dim, hidden_dim, embedding_dim, num_res_blocks=2, dropout=0.2):
        super(Encoder, self).__init__()
        
        # 初始投影层
        self.initial = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU()
        )
        
        # 残差块堆叠
        res_layers = []
        for i in range(num_res_blocks):
            res_layers.append(ResidualBlock(
                in_dim=hidden_dim, 
                out_dim=hidden_dim, 
                dropout=dropout
            ))
        self.res_blocks = nn.Sequential(*res_layers)
        
        # 输出层
        self.output = nn.Linear(hidden_dim, embedding_dim)
    
    def forward(self, x):
        x = self.initial(x)
        x = self.res_blocks(x)
        x = self.output(x)
        return x

class Decoder(nn.Module):
    """使用残差块的解码器"""
    def __init__(self, embedding_dim, hidden_dim, output_dim, num_res_blocks=2, dropout=0.2):
        super(Decoder, self).__init__()
        
        # 初始投影层
        self.initial = nn.Sequential(
            nn.Linear(embedding_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU()
        )
        
        # 残差块堆叠
        res_layers = []
        for i in range(num_res_blocks):
            res_layers.append(ResidualBlock(
                in_dim=hidden_dim, 
                out_dim=hidden_dim, 
                dropout=dropout
            ))
        self.res_blocks = nn.Sequential(*res_layers)
        
        # 输出层
        self.output = nn.Linear(hidden_dim, output_dim)
    
    def forward(self, x):
        x = self.initial(x)
        x = self.res_blocks(x)
        x = self.output(x)
        return x