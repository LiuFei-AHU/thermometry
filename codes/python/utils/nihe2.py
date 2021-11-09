"""
使用复数拟合的方式计算相位差
V1.0 BY LiuFei 2021.09.25
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from utils.util import merge_kspace_from_mrd, mrd_process, phase_correlation
from utils.model import train

# 0. set up data directory
rt = r'D:\personal\thermometry\codes\MRD_Parse\MRD\20211029\n7'
# pre-process
# mrd_process(rt, percent=60, show_image=True)  # 从MRD数据中提取幅度数据，相位数据/

# 1. start algorithm
mask_dir = os.path.join(rt, r'1\image_merged\1.npy')  # mask here
data = merge_kspace_from_mrd(rt, slice_num=0)  # get k-space data
baseline = data[:, :, 0]
heat = data[:, :, 7]
tth = train(heat, baseline)  # train the model, return phase diff
# np.save('../results/20211027/n2/20211027-n2-theat5.npy', tth)

# tth = phase_correlation(baseline, heat)
# tth = np.fliplr(tth)
# 1.1 get phase diff with itself ,as a reference
# baseline = data[:, :, 1]
# tth1 = train(heat, baseline)
# tth -= tth1

# 2. mask outer noise
# if mask_dir is not None:
#     mask = np.load(mask_dir)
#     tth[mask == 0] = 0
#     # tth1[mask == 0] = 0

# 3. rebuild temperature map,the constant should be different, here the value is just from other papers.
tth = tth * 20  # 7.781
# tth1 = tth1 * 7.781
# print(tth)

# 4. show result image
plt.figure(figsize=(10, 8))
plt.subplot(111)
h = plt.imshow(tth, vmin=0, vmax=np.max(np.max(tth)))
cb = plt.colorbar(h, label='Temperature / C')
# plt.subplot(122)
# h = plt.imshow(tth1, vmin=0, vmax=np.max(np.max(tth1)))
# plt.colorbar(h, label='Temperature / C')
plt.show()

# 5. save data ? you determine!
# savemat('theta-1009-v4.mat', {'theta': th, 'tth': tth})

