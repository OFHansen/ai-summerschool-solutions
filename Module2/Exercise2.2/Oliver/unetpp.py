# unetpp.py
import torch
from torch import nn

class ConvBlock(nn.Module):
    def __init__(self, in_ch, out_ch):
        super().__init__()
        self.seq = nn.Sequential(
            nn.Conv2d(in_ch, out_ch, 3, padding=1),
            nn.ReLU(inplace=True),
            nn.Conv2d(out_ch, out_ch, 3, padding=1),
            nn.ReLU(inplace=True),
        )
    def forward(self, x): return self.seq(x)

class Up(nn.Module):
    def __init__(self, in_ch, out_ch):
        super().__init__()
        self.up = nn.ConvTranspose2d(in_ch, out_ch, 2, 2)
    def forward(self, x): return self.up(x)

class UNetPP(nn.Module):
    def __init__(self, in_ch=1, out_ch=1):
        super().__init__()
        ch = [64, 128, 256, 512, 1024]
        self.pool = nn.MaxPool2d(2)

        # encoder
        self.conv0_0 = ConvBlock(in_ch, ch[0])
        self.conv1_0 = ConvBlock(ch[0], ch[1])
        self.conv2_0 = ConvBlock(ch[1], ch[2])
        self.conv3_0 = ConvBlock(ch[2], ch[3])
        self.conv4_0 = ConvBlock(ch[3], ch[4])

        # upsamplers shared per stage
        self.up1 = Up(ch[1], ch[0])
        self.up2 = Up(ch[2], ch[1])
        self.up3 = Up(ch[3], ch[2])
        self.up4 = Up(ch[4], ch[3])

        # nested decoder (X^{i,j})
        self.conv0_1 = ConvBlock(ch[0] + ch[0], ch[0])
        self.conv1_1 = ConvBlock(ch[1] + ch[1], ch[1])
        self.conv2_1 = ConvBlock(ch[2] + ch[2], ch[2])
        self.conv3_1 = ConvBlock(ch[3] + ch[3], ch[3])

        self.conv0_2 = ConvBlock(ch[0]*2 + ch[0], ch[0])
        self.conv1_2 = ConvBlock(ch[1]*2 + ch[1], ch[1])
        self.conv2_2 = ConvBlock(ch[2]*2 + ch[2], ch[2])

        self.conv0_3 = ConvBlock(ch[0]*3 + ch[0], ch[0])
        self.conv1_3 = ConvBlock(ch[1]*3 + ch[1], ch[1])

        self.conv0_4 = ConvBlock(ch[0]*4 + ch[0], ch[0])

        self.final = nn.Conv2d(ch[0], out_ch, 1)

    def forward(self, x):
        x0_0 = self.conv0_0(x)
        x1_0 = self.conv1_0(self.pool(x0_0))
        x2_0 = self.conv2_0(self.pool(x1_0))
        x3_0 = self.conv3_0(self.pool(x2_0))
        x4_0 = self.conv4_0(self.pool(x3_0))

        x0_1 = self.conv0_1(torch.cat([x0_0, self.up1(x1_0)], dim=1))
        x1_1 = self.conv1_1(torch.cat([x1_0, self.up2(x2_0)], dim=1))
        x2_1 = self.conv2_1(torch.cat([x2_0, self.up3(x3_0)], dim=1))
        x3_1 = self.conv3_1(torch.cat([x3_0, self.up4(x4_0)], dim=1))

        x0_2 = self.conv0_2(torch.cat([x0_0, x0_1, self.up1(x1_1)], dim=1))
        x1_2 = self.conv1_2(torch.cat([x1_0, x1_1, self.up2(x2_1)], dim=1))
        x2_2 = self.conv2_2(torch.cat([x2_0, x2_1, self.up3(x3_1)], dim=1))

        x0_3 = self.conv0_3(torch.cat([x0_0, x0_1, x0_2, self.up1(x1_2)], dim=1))
        x1_3 = self.conv1_3(torch.cat([x1_0, x1_1, x1_2, self.up2(x2_2)], dim=1))

        x0_4 = self.conv0_4(torch.cat([x0_0, x0_1, x0_2, x0_3, self.up1(x1_3)], dim=1))

        return torch.sigmoid(self.final(x0_4))
