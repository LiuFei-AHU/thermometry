"""
使用复数拟合的方式计算相位差
V1.0 BY LiuFei 2021.09.25
"""

import torch
import torch.nn as nn
import numpy as np
from numpy.fft import fftshift, ifft2
import math
from utils.util import complex_multi


class MyLoss(nn.Module):
    def __init__(self):
        super().__init__()

    def forward(self, x, y):
        err = x - y
        # anglex = torch.acos(x[:, 0] / torch.sqrt(torch.pow(x[:, 0], 2) + torch.pow(x[:, 1], 2)))
        # angley = torch.acos(y[:, 0] / torch.sqrt(torch.pow(y[:, 0], 2) + torch.pow(y[:, 1], 2)))
        # anglex = torch.sqrt(torch.pow(x[:, 0], 2) + torch.pow(x[:, 1], 2))
        # angley = torch.sqrt(torch.pow(y[:, 0], 2) + torch.pow(y[:, 1], 2))
        # anglex = x
        # angley = y
        # print(angley.shape,anglex.shape)
        # err = torch.sqrt(torch.pow(anglex - angley, 2))
        # err = anglex - angley
        # err = complexMulti(err, err) * 0.5
        # return err
        # return torch.sum(torch.sqrt(torch.pow(err[:, 0], 2) + torch.pow(err[:, 1], 2)))
        return torch.sum((torch.pow(err[:, 0], 2) + torch.pow(err[:, 1], 2)) / 2)
        # return torch.pow(err, 2)/2
        # print(err)
        # print(torch.sum(anglex - angley))
        # return torch.sum(err)
        # return torch.mean(err)
        # # print(err.T.shape, err.shape)
        # return torch.matmul(err.T, err) * 0.5
        # return torch.abs(torch.cosine_similarity(anglex, angley, dim=0))


def train(heat, baseline, max_epoch=100, fairy=False):
    """

    :param baseline: 基线数据
    :param heat: 加热数据
    :param max_epoch: 最大迭代次数
    :param fairy:
    :return:相位差
    """

    heat = fftshift(ifft2(fftshift(heat)))
    baseline = fftshift(ifft2(fftshift(baseline)))

    # load mask data, mask kspace data before process
    # mask_1 = np.load(r'D:\personal\thermometry\codes\MRD_Parse\MRD\20210930\heat\1\image_merged\1.npy')
    # mask_0 = np.load(r'D:\personal\thermometry\codes\MRD_Parse\MRD\20210930\heat\10\image_merged\1.npy')
    # heat[mask_1 == 0] = 0j
    # baseline[mask_0 == 0] = 0j
    # print(heat)

    # 初始化theta
    theta = np.random.randn(256 * 256, 2)
    theta = torch.from_numpy(theta)
    # theta_m = torch.sqrt(torch.pow(theta[:, 0], 2) + torch.pow(theta[:, 1], 2))
    # theta[:, 0] /= theta_m
    # theta[:, 1] /= theta_m
    theta.to('cpu')
    theta.requires_grad = True

    # 数据转换成 (n^2,2) 存放实部和虚部
    bl = np.zeros((256 * 256, 2))
    bl_tmp = baseline.reshape((256 * 256, 1))
    # bl_tmp /= np.absolute(bl_tmp)
    bl[:, 0] = np.squeeze(bl_tmp.real)
    bl[:, 1] = np.squeeze(bl_tmp.imag)

    bl = torch.from_numpy(bl)
    bl.to('cpu')
    bl.requires_grad = False

    ht = np.zeros((256 * 256, 2))
    ht_tmp = heat.reshape((256 * 256, 1))
    # ht_tmp /= np.absolute(ht_tmp)
    ht[:, 0] = np.squeeze(ht_tmp.real)
    ht[:, 1] = np.squeeze(ht_tmp.imag)
    ht = torch.from_numpy(ht)
    ht.to('cpu')
    ht.requires_grad = False

    criterion = MyLoss()

    opt = torch.optim.Rprop([theta], lr=0.2)
    old_loss = math.inf
    for epoch in range(max_epoch):
        # print(theta)
        opt.zero_grad()
        bl2 = complex_multi(bl, theta)

        loss = criterion(ht, bl2)
        epoch_loss = loss.item()
        print('epoch:{} loss:{}'.format(epoch, epoch_loss))
        # print('epoch:{} loss:{}'.format(epoch, loss.sum()))
        # print('epoch:{}'.format(epoch))
        loss.backward()
        # loss.backward(theta, retain_graph=True)
        # loss.backward(theta)  # torch.ones_like(theta)
        opt.step()
        if epoch > 0 and old_loss > epoch_loss and old_loss - epoch_loss < old_loss * 1e-2:
            break
        old_loss = epoch_loss

    th = theta.detach().numpy()
    print(th)
    th = np.nan_to_num(th, nan=0)

    tth = np.arccos(th[:, 0] / np.sqrt(np.power(th[:, 0], 2) + np.power(th[:, 1], 2)))
    tth = np.nan_to_num(tth, nan=0)
    tth[tth < 0] = 0
    tth = np.reshape(tth, (256, 256))
    tth = np.fliplr(tth)
    return tth
