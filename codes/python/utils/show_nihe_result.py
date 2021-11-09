import matplotlib.pyplot as plt
from scipy.io import loadmat, savemat
import numpy as np
from math import pi
from PIL import Image


m = loadmat('./theta-1009-v4.mat')
th = m['theta']

tth = np.arccos(th[:, 0] / np.sqrt(np.power(th[:, 0], 2) + np.power(th[:, 1], 2)))
tth = np.nan_to_num(tth, nan=0)
print(tth)
print(np.min(np.min(tth)), np.max(np.max(tth)))
# for i in range(len(tth)):
#     if tth[i]>pi/2:
#         tth[i] = tth[i]-pi
# tth[tth < 0] = 0
tth = np.reshape(tth, (256, 256))
tth = np.fliplr(tth)
# # uncomment these codes to perform mask on theta
dir = r'D:\personal\thermometry\codes\MRD_Parse\MRD\20211009\vertical_4\6\image_merged\1.npy'
mask = np.load(dir)
tth[mask == 0] = 0
print(mask)

# tth = np.fliplr(tth)
# calculate the average angle drift in the heated area
# 125+-15,150+-15
# h_t = tth[135:165, 110:140]
# h_t = tth[110:140, 135:165]
# h_t = tth[90:110, 110:130]
# # print(h_t)
# print('min/max:', np.min(h_t), np.max(h_t))
# h_per = np.percentile(h_t, 30)
# print('percentile：', h_per)
# h_mean = np.mean(h_t)
# print('mean:', h_mean)
# h_t[h_t < h_per] = 0

# hh = np.array([x for x in np.nditer(h_t, order='C') if x > h_per])
# hh = np.mean(hh)
# print('mean2:', hh)

# print('mean3:', np.percentile(tth, 60))
# hh = np.percentile(tth, 60)

# tth[tth <= h_per] = 0
# tth[tth <= hh] = 0
tth = 7.781 * tth
# print(tth)

fig = plt.figure(figsize=(10, 8))
h = plt.imshow(tth, vmin=0, vmax=np.max(np.max(tth)))
cb = plt.colorbar(h, label='Temperature / C')
plt.show()

#显示中心区域
# h_t *= 7.781
# fig = plt.figure(figsize=(10, 8))
# h = plt.imshow(h_t, vmin=0, vmax=np.max(np.max(h_t)))
# cb = plt.colorbar(h, label='Temperature / C')
# plt.show()


# savemat('theta_mask.mat', {'theta': th, 'tth': tth})
# img_dst = Image.fromarray(h_t)
# plt.subplot(111)
# plt.imshow(img_dst, vmin=0, vmax=np.max(img_dst), cmap='gray')
# plt.show()