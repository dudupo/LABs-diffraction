import cv2 as cv
import numpy as np
from imageio import imread, imsave
from os import walk, system
from random import randint
import matplotlib.pyplot as plt
import pickle5 as pickle
from numpy.lib.type_check import imag
from scipy.stats import describe
from skimage.color import rgb2gray
import matplotlib.patches as mpatches
from skimage import data
from skimage.filters import threshold_otsu, sobel
# from skimage.segmentation import clear_border, watershed, expand_labels
from skimage.measure import label, regionprops
from skimage.morphology import closing, square, dilation
from skimage.color import label2rgb
import scipy.ndimage as ndimage
from random import choice, gammavariate
from copy import deepcopy
from re import S
import numpy as np
from numpy.lib import vectorize
import scipy.special
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.path.realpath(".."))

pkl_path = '/cs/usr/michael_cohen/Desktop/LAB/LABs-diffraction/bio/pkl/'


def colony_k(colonies, i, plot=False):
    colony = colonies[i]
    col0 = colony[0].pix_num
    k_arr = []
    bad = False
    for j, col_t in enumerate(colony):
        k = int((col_t.pix_num/col0) + 0.5)
        if j != 0:
            mult = k / k_arr[j-1]
            if mult != 2:
                bad = True
                break
        else:
            if col_t.pix_num > 550 or col_t.pix_num < 135:
                bad = True
                break
        k_arr.append(k)
    if plot and not bad:
        plt.plot(k_arr)
        # plt.show()
    return k_arr, bad

def get_multi_factor_pckl(t, col_idx=1):
    k = 1
    colonies = pickle.load( open("colonys-prob_test-2021-06-21_14-27-21.872035.pkl", "rb"))
    pkt = np.zeros((k, t))
    col_to_plot = []
    col_0 = 0
    for i, colony in enumerate(colonies):
        colony_zero = colony[0].pix_num
        start_idx = 0
        found = False
        for j, col_t in enumerate(colony):
            if col_t.pix_num // colony[0].pix_num < 2:
                continue
            elif not found:
                colony_zero = col_t.pix_num
                found = True
                start_idx = j
            if j < t and found:
                if i == col_idx:
                    col_to_plot.append(col_t.pix_num/ colony[0].pix_num)
                cur_k = int(col_t.pix_num // colony_zero)
                if cur_k <= k:
                    pkt[cur_k - 1][j - start_idx] += 1
                else:
                    pkt = np.r_[pkt, np.zeros((cur_k - k, t))]
                    k = cur_k
                    pkt[cur_k - 1][j - start_idx] += 1

    if len(col_to_plot) > 0:
        plt.plot(col_to_plot)
        plt.show()
    else:
        print("single colony plot failed, try different colony")
    return pkt


def get_multi_factor(colonies, start_time=0,  final_time=80):
    k=1
    t = final_time - start_time
    pkt = np.zeros((k, t))
    for colony in colonies:
        for j, col_t in enumerate(colony):
            if j < t:
                cur_k = int((col_t.pix_num/colony[0].pix_num) + 0.5)
                if cur_k <= k:
                    pkt[cur_k-1][j] += 1
                else:
                    pkt = np.r_[pkt, np.zeros((cur_k-k, t))]
                    k = cur_k
                    pkt[cur_k-1][j] += 1
    return pkt



def calc_avg_t0_size(colonies):
    outliers = 0
    outliers_idxs = []
    sum = 0
    l = len(colonies) + 1
    for i, col in enumerate(colonies):
        t0_size = col[0].pix_num
        sum += t0_size
        if t0_size > 550 or t0_size < 135:
            outliers += 1
            outliers_idxs.append(i)
    print("num of outliers: ", outliers)
    print("outliers indices: ", outliers_idxs)
    print(l)
    return sum/l



if __name__ == '__main__':
    with open("colonys-prob_test-2021-06-21_14-27-21.872035.pkl", "rb") as fh:
        colonies = pickle.load(fh)
    good_colonies_k = []
    good_colonies = []
    for i in range(len(colonies)):
        col, bad = colony_k(colonies, i, True)
        if not bad:
            good_colonies_k.append(col)
            good_colonies.append(colonies[i])

    plt.show()
    # print(calc_avg_t0_size(good_colonies))
    exit(0)