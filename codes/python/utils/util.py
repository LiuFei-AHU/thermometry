"""
核磁测温工具类代码，支持其他代码工作
V1.0 BY LiuFei 2021.09.25
---------------------------
数据文件组织方式必须为：
root-
     -(d)xx (one of MRI scans)
        -(d)image_merged (images fused from different channels)
            -(f)1.npy merged (image with noise removed)
        -(d)kspace_data_from_mrd (kspace_data attained from MRD in different channels)
        -(d)mrd_parsed_images (images parsed from MRD in different channels)
        -(f)Temp1-8.mrd (original MRD DATA)
    -xx ...
"""

from scipy.io import loadmat, savemat
import os
import numpy as np
from scipy.stats import gaussian_kde as kdf
from scipy import special as sp
from PIL import Image
import matplotlib.pyplot as plt
import cv2
import torch


def mrd_3d(filename, pathname):
    """
    解析MRD数据，从matlab将代码迁移到Python
    :param filename: MRD文件名
    :param pathname: 路径名
    :return:
    """
    reordering = 'seq'

    # 用于保存提取的图片
    curr = os.path.join(pathname, 'mrd_pared_images')
    if not os.path.exists(curr):
        os.mkdir(curr)

    # 用于保存提取的kspace数据
    kspace_filename = os.path.join(pathname, 'kspace_data_from_mrd')
    if not os.path.exists(kspace_filename):
        os.mkdir(kspace_filename)

    # disp([pathname, filename]
    f_name = os.path.join(pathname, filename)
    f = open(f_name, 'rb')  # Define the file id
    val = f.read(4)  # val为MRD文件前四个数据
    xdim = int(val[0])  # 128
    ydim = int(val[1])  # 128
    zdim = int(val[2])  # 1
    dim4 = int(val[3])  # 3

    print(xdim, ydim)

    f.seek(18)  # 从文件开始(Beginning of file)偏移18个字节的位置

    f.seek(48)  # 从文件开始偏移48个字节的位置
    f.seek(152)  # 从文件开始偏移152个字节的位置
    val = f.read(2)
    dim5 = int(val[0])  # 1
    dim6 = int(val[1])
    f.seek(256)  # 从文件开始偏移256个字节的位置
    no_samples = xdim  # NOS - X: 12
    no_views = ydim  # NOV - Y: 128
    no_views_2 = zdim  # NOV2 - Z: 1
    no_slices = dim4  # number of slices: 3
    no_echoes = dim5  # 回波: 1
    no_expts = dim6  # 输出: 1

    onlydatatype = 3  # onlydatatype = '3'
    iscomplex = 2  # 是否为复数：2为复数
    dataformat = 'int16'

    num2read = no_expts * no_echoes * no_slices * no_views_2 * no_views * no_samples * iscomplex  # *datasize;
    m_total = f.read(num2read)  # reading all the data at once
    print(m_total)


def merge_kspace_from_mrd(root_dir, slice_num=0):
    """
    从MRD中提取Kspace数据，直接叠加融合
    用作测温算法数据来源

    :param root_dir: 数据根目录
    :param slice_num: 切片位置,默认第一张切片
    :return: data 返回融合后的复数矩阵，尺寸为（n,n,s）
    """

    # r'E:\BaiduNetdiskDownload\20210903\{}\kspace_data_from_mrd'
    # dir = os.path.join(root_dir, r'{}\kspace_data_from_mrd')

    files = os.listdir(root_dir)
    files.sort(key=lambda x: int(x[:]))

    data = np.zeros((256, 256, len(files) + 1), np.complex_)

    for i, nm in enumerate(files):
        print(nm)
        d = os.path.join(root_dir, str(nm), 'kspace_data_from_mrd')
        one_data = None
        for idx, dd in enumerate(os.listdir(d)):
            # print(i, dd)
            f_name = os.path.join(d, dd)
            m = loadmat(f_name)
            if m['KspaceData'].ndim == 3:
                tmp = m['KspaceData'][:, :, slice_num]
            else:
                tmp = m['KspaceData'][:, :]
            if idx == 0:
                one_data = tmp
            else:
                one_data += tmp
            # if idx == 0:
            #     one_data = m['KspaceData'][:, :, slice_num]
            # else:
            #     one_data += m['KspaceData'][:, :, slice_num]
        data[:, :, i + 1] = one_data
    # print('data.shape',data.shape,data.shape[2])
    data[:, :, 0] = data[:, :, data.shape[2] - 1]
    return data


def merge_image_from_mrd(root_dir, percent=80):
    """
    从MRD数据中直接提取幅度数据，按比例融合，过滤数值，保存
    可用作mask
    :return: None
    """
    assert root_dir is not None, 'root_dir must not be none'

    files = os.listdir(root_dir)
    files.sort(key=lambda x: int(x[:]))

    for l in files:
        # d = dir.format(i)
        d = os.path.join(root_dir, str(l), 'kspace_data_from_mrd')
        imgs = []
        maxi = 0
        for idx, dd in enumerate(os.listdir(d)):
            f_name = os.path.join(d, dd)
            m = loadmat(f_name)
            if m['KspaceData'].ndim == 3:
                one_data = m['KspaceData'][:, :, 0]
            else:
                one_data = m['KspaceData'][:, :]
            one_data = np.fliplr(np.fft.ifftshift(np.fft.ifft2(one_data)))

            mag = np.absolute(one_data)
            # mag[mag < 0.008] = 0
            # print(mag)
            # print(np.max(np.max(mag)), np.min(np.min(mag)))
            imgs.append(mag)
            # maxj = np.max(np.max(mag))
            # if maxj > maxi:
            #     maxi = maxj
            # print(maxi)
        # do merge 8->4
        dst = []
        for i in range(0, len(imgs), 2):
            img1 = imgs[i]
            img2 = imgs[i + 1]
            dst.append(img1 * 0.5 + img2 * 0.5)
        # do merhge 4->2
        dst1 = []
        for i in range(0, len(dst), 2):
            img1 = dst[i]
            img2 = dst[i + 1]
            dst1.append(img1 * 0.5 + img2 * 0.5)

        img_dst = dst1[0] * 0.5 + dst1[1] * 0.5
        maxj = np.max(np.max(img_dst))
        if maxj > maxi:
            maxi = maxj
        im_per = np.percentile(img_dst, percent)
        print(im_per)
        # img_dst[img_dst < 0.01] = 0
        img_dst[img_dst < im_per] = 0
        # print(img_dst)
        # 保存
        s_path = os.path.join(root_dir, str(l), 'image_merged')
        if not os.path.exists(s_path):
            os.makedirs(s_path)
        with open(s_path + '/1.npy', mode='wb') as f:
            np.save(f, img_dst)
        # Image.fromarray(np.array(img_dst)).show()
        # plt.subplot(111)
        # plt.imshow(np.array(img_dst), vmin=0, vmax=im_per, cmap='gray')
        # plt.show()

    return None


def show_image_from_npy(root_dir):
    """
    显示过滤幅度图噪声后的幅度图像，9宫格形状
    :return: None
    """
    # root_dir = r'E:\BaiduNetdiskDownload\20210903'
    assert root_dir is not None, 'root_dir must not be none'

    files = os.listdir(root_dir)
    files.sort(key=lambda x: int(x[:]))

    fig = plt.figure()
    for l, nm in enumerate(files):
        if l == 9:
            break
        print(l)
        d = os.path.join(root_dir, str(nm), 'image_merged', '1.npy')
        img = np.load(d)
        maxi = np.max(np.max(img))
        plt.subplot(330 + l + 1)
        plt.imshow(np.array(img), vmin=0, vmax=maxi, cmap='gray')
        plt.axis('off')
        plt.subplots_adjust(wspace=0, hspace=0)
        # plt.savefig('merger_images.jpg', bbox_inches='tight')
    fig.tight_layout()
    plt.show()


def merge_images_from_parsed_img(root_dir):
    """
    将解析出来的8通道的核磁图片进行融合，即8->1
    :param root_dir: 数据存放的根目录
    :return: None
    """

    # root_dor = r'E:\BaiduNetdiskDownload\20210903'
    assert root_dir is not None, 'root_dir must not be none'

    files = os.listdir(root_dir)
    files.sort(key=lambda x: int(x[:]))

    # scans from 4-14
    for l in files:
        dir = os.path.join(root_dir, str(l), 'mrd_pared_images')
        # 4 slices
        for s in range(1, 2):  # 5
            # 8 channels
            imgs = []
            for i in range(1, 9):
                img = cv2.imread(os.path.join(dir, 'Temp' + str(i) + '.MRD' + str(s) + '.JPG'))
                imgs.append(img)

            # do merge 8->4
            dst = []
            for i in range(0, len(imgs), 2):
                img1 = imgs[i]
                img2 = imgs[i + 1]
                dst.append(cv2.addWeighted(img1, 0.5, img2, 0.5, 0))
            # do merhge 4->2
            dst1 = []
            for i in range(0, len(dst), 2):
                img1 = dst[i]
                img2 = dst[i + 1]
                dst1.append(cv2.addWeighted(img1, 0.5, img2, 0.5, 0))

            img_dst = cv2.addWeighted(dst1[0], 0.5, dst1[1], 0.5, 0)

            # save merged image
            s_path = os.path.join(root_dir, str(l), 'image_merged')
            if not os.path.exists(s_path):
                os.makedirs(s_path)
            # img_dst = img_dst[225:225+280,305:305+280]
            # # 增强对比度（图片变亮）
            # img_dst = exposure.adjust_gamma(img_dst, 1)
            # img_dst = exposure.rescale_intensity(img_dst)
            cv2.imwrite(s_path + '/{}-MRD{}.jpg'.format(l, s), img_dst)
            print(l, s, 'processed')

    # import numpy as np
    # import cv2
    # from PIL import Image, ExifTags
    #
    #
    # def calcMeanAndVariance(img):
    #     row = img.shape[0]
    #     col = img.shape[1]
    #     # channel=img.shape[2]
    #     total = row * col
    #     print(row, col, total)
    #     mean = np.zeros((3))
    #     variance = np.zeros((3))
    #     sum = np.zeros((3))
    #
    #     for i in range(row):
    #         for j in range(col):
    #             sum[0] += img[i][j][0]
    #             sum[1] += img[i][j][1]
    #             sum[2] += img[i][j][2]
    #
    #     mean[0] = sum[0] / total
    #     mean[1] = sum[1] / total
    #     mean[2] = sum[2] / total
    #     sum = np.zeros((3))
    #     for i in range(row):
    #         for j in range(col):
    #             sum[0] = np.square(img[i][j][0] - mean[0])
    #             sum[1] = np.square(img[i][j][1] - mean[1])
    #             sum[2] = np.square(img[i][j][2] - mean[2])
    #
    #     variance[0] = np.sqrt(sum[0] / total)
    #     variance[1] = np.sqrt(sum[1] / total)
    #     variance[2] = np.sqrt(sum[2] / total)
    #     print(mean, variance)
    #     return mean, variance
    #
    #
    # def cololTransit(img1, img2):
    #     # image1 = cv2.cvtColor(img1, cv2.CV_8U)
    #     # image2 = cv2.cvtColor(img2, cv2.CV_8U)
    #     image1, image2 = img1, img2
    #     mean1, variance1 = calcMeanAndVariance(image1)
    #     mean2, variance2 = calcMeanAndVariance(image2)
    #     # print (mean1,variance1)
    #     radio = np.zeros((3))
    #
    #     radio[0] = variance2[0] / variance1[0]
    #     radio[1] = variance2[1] / variance1[1]
    #     radio[2] = variance2[2] / variance1[2]
    #
    #     print('test', radio)
    #
    #     row = image1.shape[0]
    #     col = image1.shape[1]
    #     for i in range(row):
    #         for j in range(col):
    #             image1[i][j][0] = min(255, max(0, radio[0] * (image1[i][j][0] - mean1[0]) + mean2[0]))
    #             image1[i][j][1] = min(255, max(0, radio[1] * (image1[i][j][1] - mean1[1]) + mean2[1]))
    #             image1[i][j][2] = min(255, max(0, radio[2] * (image1[i][j][2] - mean1[2]) + mean2[2]))
    #     # image = cv2.cvtColor(image1, cv2.CV_8U)
    #     image = image1
    #     return image
    #
    #
    # if __name__ == '__main__':
    #     img1 = cv2.imread(r'E:\BaiduNetdiskDownload\20210903\5\mrd_pared_images\Temp1.MRD1.jpg')
    #     img2 = cv2.imread(r'E:\BaiduNetdiskDownload\20210903\5\mrd_pared_images\Temp2.MRD1.jpg')
    #     img3 = cv2.imread(r'E:\BaiduNetdiskDownload\20210903\5\mrd_pared_images\Temp3.MRD1.jpg')
    #     img4 = cv2.imread(r'E:\BaiduNetdiskDownload\20210903\5\mrd_pared_images\Temp4.MRD1.jpg')
    #     img5 = cv2.imread(r'E:\BaiduNetdiskDownload\20210903\5\mrd_pared_images\Temp5.MRD1.jpg')
    #     img6 = cv2.imread(r'E:\BaiduNetdiskDownload\20210903\5\mrd_pared_images\Temp6.MRD1.jpg')
    #     img7 = cv2.imread(r'E:\BaiduNetdiskDownload\20210903\5\mrd_pared_images\Temp7.MRD1.jpg')
    #     img8 = cv2.imread(r'E:\BaiduNetdiskDownload\20210903\5\mrd_pared_images\Temp8.MRD1.jpg')
    #
    #     # cv2.namedWindow('src')
    #     # cv2.namedWindow('dst')
    #     # cv2.resizeWindow('src',500,500)
    #     # cv2.resizeWindow('dst',500,500)
    #     # cv2.imshow('src', img1)
    #     # cv2.imshow('dst', img2)
    #     # cv2.waitKey()
    #     # cv2.destroyAllWindows()
    #
    #     img = cololTransit(img1, img2)
    #     img = cololTransit(img, img3)
    #     # img = cololTransit(img, img4)
    #     # img = cololTransit(img, img5)
    #     # img = cololTransit(img, img6)
    #     # img = cololTransit(img, img7)
    #     # img = cololTransit(img, img8)
    #     cv2.namedWindow('result')
    #     cv2.imshow('result', img)
    #     cv2.waitKey()
    #     cv2.destroyAllWindows()

    # img1 = cv2.imread(r'E:\BaiduNetdiskDownload\20210903\5\mrd_pared_images\Temp1.MRD1.jpg')
    # img2 = cv2.imread(r'E:\BaiduNetdiskDownload\20210903\5\mrd_pared_images\Temp2.MRD1.jpg')
    # img3 = cv2.imread(r'E:\BaiduNetdiskDownload\20210903\5\mrd_pared_images\Temp3.MRD1.jpg')
    # img4 = cv2.imread(r'E:\BaiduNetdiskDownload\20210903\5\mrd_pared_images\Temp4.MRD1.jpg')
    # img5 = cv2.imread(r'E:\BaiduNetdiskDownload\20210903\5\mrd_pared_images\Temp5.MRD1.jpg')
    # img6 = cv2.imread(r'E:\BaiduNetdiskDownload\20210903\5\mrd_pared_images\Temp6.MRD1.jpg')
    # img7 = cv2.imread(r'E:\BaiduNetdiskDownload\20210903\5\mrd_pared_images\Temp7.MRD1.jpg')
    # img8 = cv2.imread(r'E:\BaiduNetdiskDownload\20210903\5\mrd_pared_images\Temp8.MRD1.jpg')
    #
    # dst1 = cv2.addWeighted(img1, 0.5, img2, 0.5, 0)
    # dst2 = cv2.addWeighted(img3, 0.5, img4, 0.5, 0)
    # dst3 = cv2.addWeighted(img5, 0.5, img6, 0.5, 0)
    # dst4 = cv2.addWeighted(img7, 0.5, img8, 0.5, 0)
    #
    # dst5 = cv2.addWeighted(dst1, 0.5, dst2, 0.5, 0)
    # dst6 = cv2.addWeighted(dst3, 0.5, dst4, 0.5, 0)
    #
    # dst7 = cv2.addWeighted(dst5, 0.5, dst6, 0.5, 0)
    #
    # # cv2.imshow('dst1', dst1)
    # # cv2.imshow('dst2', dst2)
    # # cv2.imshow('dst3', dst3)
    # # cv2.imshow('dst4', dst4)
    #
    # cv2.imshow('dst5', dst5)
    # cv2.imshow('dst6', dst6)
    # cv2.imshow('dst7', dst7)
    # cv2.waitKey(0)
    # cv2.destroyAllWindows()


def mrd_process(dir, percent=80, show_image=False):
    merge_image_from_mrd(dir, percent)  # 从MRD数据中直接提取幅度数据，mask
    merge_kspace_from_mrd(dir)  # 从MRD中提取Kspace数据
    merge_images_from_parsed_img(dir)  # 将解析出来的8通道的核磁图片进行融合
    if show_image:
        show_image_from_npy(dir)  # 显示过滤幅度图噪声后的幅度图像


def phase_correlation(a, b):
    """
    计算两个复数的相位差
    :param a:
    :param b:
    :return:
    """
    G_a = np.fft.fft2(a)
    G_b = np.fft.fft2(b)
    G_b = np.fft.fftshift(G_b)
    G_a = np.fft.fftshift(G_a)

    conj_b = np.ma.conjugate(G_b)
    R = G_a * conj_b
    R /= np.absolute(R)
    R = np.angle(R)
    # r = np.fft.ifft2(R)
    # r = np.maximum(r, 0)
    # R = np.fft.fftshift(R)
    return R  # r.real


def ricernd_noise(v, s):
    """
    莱斯噪声
    :param v:
    :param s:
    :return: 噪声
    """
    dim = v.shape
    print(dim)
    x = s * np.random.randn(*dim) + v
    y = s * np.random.randn(*dim)
    r = np.sqrt(x ** 2 + y ** 2)
    return r


class Rice:
    """
    Rice 噪声生成器
    """
    # numSamples = 2*(10**6)  # the number of samples used in the simulation
    # numSamples = 64 #表示要产生的随机数的长度
    r = np.linspace(0, 6, 6000)  # theoretical envelope PDF x axes
    theta = np.linspace(-np.pi, np.pi, 6000)  # theoretical phase PDF x axes

    def __init__(self, K, r_hat_2, phi, numSamples):

        # # user input checks and assigns value
        self.K = self.input_Check(K, "K", 0, 50)
        self.r_hat_2 = self.input_Check(r_hat_2, "\hat{r}^2", 0.5, 2.5)
        self.phi = self.input_Check(phi, "\phi", -np.pi, np.pi)
        self.numSamples = numSamples
        # user input checks and assigns value

        # simulating and theri densities
        self.multipathFading = self.complex_Multipath_Fading()
        self.xdataEnv, self.ydataEnv = self.envelope_Density()
        self.xdataPh, self.ydataPh = self.phase_Density()

        # theoretical PDFs calculated
        self.envelopeProbability = self.envelope_PDF()
        self.phaseProbability = self.phase_PDF()

    # 输入格式检查
    def input_Check(self, data, inputName, lower, upper):
        # input_Check checks the user inputs

        # has a value been entered
        if data == "":
            raise ValueError(" ".join((inputName, "must have a numeric value")))

        # incase of an non-numeric input
        try:
            data = float(data)
        except:
            raise ValueError(" ".join((inputName, "must have a numeric value")))

        # data must be within the range
        if data < lower or data > upper:
            raise ValueError(" ".join((inputName, f"must be in the range [{lower:.2f}, {upper:.2f}]")))

        return data

    def calculate_Means(self):
        # calculate_means calculates the means of the complex Gaussians representing the
        # in-phase and quadrature components

        p = np.sqrt(self.K * self.r_hat_2 / (1 + self.K)) * np.cos(self.phi)
        q = np.sqrt(self.K * self.r_hat_2 / (1 + self.K)) * np.sin(self.phi)

        return p, q

    def scattered_Component(self):
        # scattered_Component calculates the power of the scattered signal component

        sigma = np.sqrt(self.r_hat_2 / (2 * (1 + self.K)))

        return sigma

    def generate_Gaussians(self, mean, sigma):
        # generate_Gaussians generates the Gaussian random variables

        gaussians = np.random.default_rng().normal(mean, sigma, self.numSamples)

        return gaussians

    def complex_Multipath_Fading(self):
        # complex_Multipath_Fading generates the complex fading random variables

        p, q = self.calculate_Means()
        sigma = self.scattered_Component()
        multipathFading = self.generate_Gaussians(p, sigma) + (1j * self.generate_Gaussians(q, sigma))

        return multipathFading

    def envelope_PDF(self):
        # envelope_PDF calculates the theoretical envelope PDF

        PDF = 2 * (1 + self.K) * self.r / self.r_hat_2 \
              * np.exp(- self.K - ((1 + self.K) * self.r ** 2) / self.r_hat_2) \
              * np.i0(2 * self.r * np.sqrt(self.K * (1 + self.K) / self.r_hat_2))

        return PDF

    def phase_PDF(self):
        # phase_PDF calculates the theoretical phase PDF

        def q_func(x):
            # Q-function

            return 0.5 - 0.5 * sp.erf(x / np.sqrt(2))

        PDF = (1 / (2 * np.pi)) * np.exp(- self.K) * (1 + (np.sqrt(4 * np.pi * self.K)
                                                           * np.exp(
                    self.K * (np.cos(self.theta - self.phi)) ** 2) * np.cos(self.theta - self.phi))
                                                      * (1 - q_func(
                    np.sqrt(2 * self.K) * np.cos(self.theta - self.phi))))

        return PDF

    def envelope_Density(self):
        # envelope_Density finds the envelope PDF of the simulated random variables

        R = np.sqrt((np.real(self.multipathFading)) ** 2 + (np.imag(self.multipathFading)) ** 2)
        kde = kdf(R)
        x = np.linspace(R.min(), R.max(), 100)
        p = kde(x)

        return x, p

    def phase_Density(self):
        # phase_Density finds the phase PDF of the simulated random variables

        R = np.angle(self.multipathFading)
        kde = kdf(R)
        x = np.linspace(R.min(), R.max(), 100)
        p = kde(x)

        return x, p


def complex_multi(a, b):
    """
    复数乘法
    :param a: 张量a[rel,img]
    :param b: 张量b[rel,img]
    :return: 复数乘法
    """
    r = a.shape[0]
    c = torch.zeros([r, 2])
    c[:, 0] = a[:, 0] * b[:, 0] - a[:, 1] * b[:, 1]
    c[:, 1] = a[:, 0] * b[:, 1] + a[:, 1] * b[:, 0]
    # for i in range(r):
    #     c[i, 0] = a[i, 0] * b[i, 0] - a[i, 1] * b[i, 1]
    #     c[i, 1] = a[i, 0] * b[i, 1] + a[i, 1] * b[i, 0]
    return c


def show_phase_temp_curve(x=None, y=None, plot='plot'):
    """
     画出相位差和温度差的关系曲线
     :param x: 相位差
     :param y: 温度差
     :param plot: 散点图 or 折现图 => ['plot', 'scatter']
     :return:
     """

    # X轴取值相位差，Y轴取值温度差
    if x is None:
        # 以2021.10.27 n2的数据为例，计算相位差和温度差的关系，画出关系曲线
        print('phase data is required! \n This is an example data from 2021.10.27 n2')
        x = [0.911758261,
             0.891536749,
             0.700557613,
             0.606077677,
             # 0.787883234,
             # 0.829426096,
             0.335970884,
             0.325331918,
             0.262157513
             ]
        y = [16.9,
             14.3,
             11.5,
             8.7,
             # 6.9,
             # 5.1,
             3.3,
             2.4,
             0.9
             ]

    # 坐标轴范围
    plt.xlim(0, 1)
    plt.ylim(0, 25)

    # 画出相位差-温度差曲线
    if plot == 'plot':
        plt.plot(x, y, c='red', marker='o', label='phase temp mapping')
    else:
        plt.scatter(x, y, c='red', marker='o', label='phase temp mapping')
    plt.ylabel('Temperature/℃')
    plt.xlabel('Phase/diff')
    plt.title("Curve of temperature variation")
    plt.legend(loc='best')

    plt.show()


if __name__ == '__main__':
    dir = r'D:\personal\thermometry\codes\MRD_Parse\MRD\20210903'
    merge_image_from_mrd(dir, percent=80)  # 从MRD数据中直接提取幅度数据，mask
    merge_kspace_from_mrd(dir)  # 从MRD中提取Kspace数据
    merge_images_from_parsed_img(dir)  # 将解析出来的8通道的核磁图片进行融合
    show_image_from_npy(dir)  # 显示过滤幅度图噪声后的幅度图像
    # mrd_3d('Temp1.MRD', r'D:\personal\thermometry\codes\MRD_Parse\MRD\20211014\hrizontal1\1')
