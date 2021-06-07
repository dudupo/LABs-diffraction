
import cv2 as cv
import numpy as np
from imageio import imread, imsave
from os import walk, system
import matplotlib.pyplot as plt
from scipy.stats import describe
from skimage.color import rgb2gray
import matplotlib.patches as mpatches
from skimage import data
from skimage.filters import threshold_otsu, sobel
from skimage.segmentation import clear_border, watershed, expand_labels

from skimage.measure import label, regionprops
from skimage.morphology import closing, square
from skimage.color import label2rgb



def read_image(path, gray_output):
    im = imread(path)
    im_float = im.astype(np.float64)/255
    if gray_output == 1:
        if im.ndim == 2:
            return im_float
        if im.ndim == 3:
            return rgb2gray(im_float)
    else:
        return im_float


def im_disp(im):
    plt.imshow(make_binary(im), cmap='gray')
    plt.show()




def make_binary(img):
    thresh = 1/255
    img[img <= thresh ] = 0
    img[img > thresh ] = 1
    return img


def find_regions(image):

    # apply threshold
    thresh = threshold_otsu(image)
    bw = closing(image > thresh, square(3))

    # remove artifacts connected to image border
    cleared = clear_border(bw)

    # label image regions
    label_image = label(cleared, connectivity=2)
    # to make the background transparent, pass the value of `bg_label`,
    # and leave `bg_color` as `None` and `kind` as `overlay`
    image_label_overlay = label2rgb(label_image, image=image, bg_label=0)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.imshow(image_label_overlay)

    for region in regionprops(label_image):
        # take regions with large enough areas
        if region.area >= 100:
            # draw rectangle around segmented coins
            minr, minc, maxr, maxc = region.bbox
            rect = mpatches.Rectangle((minc, minr), maxc - minc, maxr - minr,
                                      fill=False, edgecolor='red', linewidth=2)
            ax.add_patch(rect)

    ax.set_axis_off()
    plt.tight_layout()
    plt.show()

# def find_expand_regions(im):
#     coins = data.coins()
#
#     # Make segmentation using edge-detection and watershed.
#     edges = sobel(coins)
#
#     # Identify some background and foreground pixels from the intensity values.
#     # These pixels are used as seeds for watershed.
#     markers = np.zeros_like(coins)
#     foreground, background = 1, 2
#     markers[coins < 30.0] = background
#     markers[coins > 150.0] = foreground
#
#     ws = watershed(edges, markers)
#     seg1 = label(ws == foreground)
#
#     expanded = expand_labels(seg1, distance=10)
#
#     # Show the segmentations.
#     fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(9, 5),
#                              sharex=True, sharey=True)
#
#     color1 = label2rgb(seg1, image=coins, bg_label=0)
#     axes[0].imshow(color1)
#     axes[0].set_title('Sobel+Watershed')
#
#     color2 = label2rgb(expanded, image=coins, bg_label=0)
#     axes[1].imshow(color2)
#     axes[1].set_title('Expanded labels')
#
#     for a in axes:
#         a.axis('off')
#     fig.tight_layout()
#     plt.show()


if __name__ == '__main__':
    im = read_image('./test/img_000000100_YFPFast_000.tif', 1)
    im = make_binary(im)
    # im_disp(im)
    find_regions(im)





