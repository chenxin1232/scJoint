import torch
import torch.nn as nn
import torchvision.models as models
import matplotlib.pyplot as plt

class VGGUNet(nn.Module):
    def __init__(self, num_classes=2):
        super(VGGUNet, self).__init__()
        vgg = models.vgg16_bn(pretrained=False)
        features = list(vgg.features.children())

        self.enc1 = nn.Sequential(*features[:6])    # 64
        self.enc2 = nn.Sequential(*features[6:13])  # 128
        self.enc3 = nn.Sequential(*features[13:23]) # 256
        self.enc4 = nn.Sequential(*features[23:33]) # 512
        self.enc5 = nn.Sequential(*features[33:43]) # 512

        self.center = nn.Sequential(
            nn.Conv2d(512, 1024, 3, padding=1),
            nn.ReLU(),
            nn.Conv2d(1024, 512, 3, padding=1),
            nn.ReLU()
        )

        def up_block(in_c, out_c):
            return nn.Sequential(
                nn.ConvTranspose2d(in_c, out_c, kernel_size=2, stride=2),
                nn.ReLU()
            )

        self.up5 = up_block(512, 512)
        self.up4 = up_block(512, 256)
        self.up3 = up_block(256, 128)
        self.up2 = up_block(128, 64)

        self.dec5 = nn.Conv2d(1024, 512, 3, padding=1)
        self.dec4 = nn.Conv2d(512, 256, 3, padding=1)
        self.dec3 = nn.Conv2d(256, 128, 3, padding=1)
        self.dec2 = nn.Conv2d(128, 64, 3, padding=1)

        self.final = nn.Conv2d(64, num_classes, kernel_size=1)

    def forward(self, x):
        e1 = self.enc1(x)
        e2 = self.enc2(e1)
        e3 = self.enc3(e2)
        e4 = self.enc4(e3)
        e5 = self.enc5(e4)

        c = self.center(e5)

        d5 = self.up5(c)
        d5 = torch.cat([d5, e4], dim=1)
        d5 = self.dec5(d5)

        d4 = self.up4(d5)
        d4 = torch.cat([d4, e3], dim=1)
        d4 = self.dec4(d4)

        d3 = self.up3(d4)
        d3 = torch.cat([d3, e2], dim=1)
        d3 = self.dec3(d3)

        d2 = self.up2(d3)
        d2 = torch.cat([d2, e1], dim=1)
        d2 = self.dec2(d2)

        out = self.final(d2)
        return out

# ---------- 主函数 ----------
if __name__ == "__main__":
    torch.manual_seed(42)
    model = VGGUNet(num_classes=2)

    # 生成 1 张图：模拟图像 + 简单前景区域
    input_tensor = torch.randn(1, 3, 256, 256)
    target_mask = torch.zeros(1, 256, 256).long()
    target_mask[:, 100:150, 100:150] = 1  # 设置前景块

    # 损失函数和优化器
    criterion = nn.CrossEntropyLoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-5)

    # 多轮训练
    model.train()
    for epoch in range(10):
        output = model(input_tensor)
        loss = criterion(output, target_mask)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        print(f"Epoch {epoch+1}/10, Loss: {loss.item():.4f}")

    # 可视化预测
    model.eval()
    with torch.no_grad():
        output = model(input_tensor)
        pred_mask = torch.argmax(output, dim=1).squeeze().numpy()  # (H, W)
        input_image = input_tensor[0].permute(1, 2, 0).numpy()
        gt_mask = target_mask.squeeze().numpy()

    # 绘图展示
    plt.figure(figsize=(12, 4))
    plt.subplot(1, 3, 1)
    plt.imshow(input_image[..., 0], cmap='gray')  # 只展示R通道
    plt.title("Input Image (R channel)")
    plt.axis("off")

    plt.subplot(1, 3, 2)
    plt.imshow(gt_mask, cmap='gray')
    plt.title("Target Mask")
    plt.axis("off")

    plt.subplot(1, 3, 3)
    plt.imshow(pred_mask, cmap='gray')
    plt.title("Predicted Mask")
    plt.axis("off")

    plt.tight_layout()
    plt.show()
