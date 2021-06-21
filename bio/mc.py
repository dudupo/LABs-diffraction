
import cv2 as cv
import numpy as np
from imageio import imread, imsave
from os import walk, system
from random import randint
import matplotlib.pyplot as plt
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
    plt.imshow(make_binary_plot(im), cmap='gray')
    plt.show()


def make_binary(img):
    thresh = 20/255
    img = ndimage.gaussian_filter(img, sigma=2)
    img[img <= thresh ] = 0
    img[img > thresh ] = 1
    return img

def make_binary_plot(img):
    thresh = 20/255
    img = ndimage.gaussian_filter(img, sigma=2)
    img[img <= thresh ] = 0
    img[img > thresh ] = 1
    return dilation(img, square(7))


def get_middle(rect):
    x, y = rect.xy
    x += rect._width//2
    y += rect._height//2
    return x,y

def compute_dist(rect1, rect2):
    return (np.sum(\
        (np.array( get_middle(rect1)) -\
         np.array( get_middle(rect2)))**2))**0.5

def match_frames(f1, f2):
    if f1 is None:
        print(f"len: ~ , {len(f2)}")
        return [[_] for _ in f2]
    print(f"len: {len(f1)} , {len(f2)}")
    for colony in f1:
        minindex = np.argmin( np.vectorize( \
            lambda _f2 : compute_dist(colony[-1], _f2)) ( f2 ) )
        f2[minindex].id = colony[0].id
        colony.append(f2[minindex])
    return f1


class ColonyRect(mpatches.Rectangle):
    def __init__(self, maxc, minc, maxr, minr, pix_num, time):
        super(ColonyRect, self).__init__((minc, minr), maxc - minc, maxr - minr,
                           fill=False, edgecolor='red', linewidth=0.5)
        self.pix_num = pix_num
        self.id = randint(1,200) 
        self.time = time


def plotColonyTrace(ax, rects, _matrix):
    # plt.image_show( inital_iamge)
    for rect in rects:
        rect.plot(ax, _matrix)

def get_rect_arr(image, time):
    image = make_binary(image)
    rect_arr = list()
    # apply threshold
    thresh = threshold_otsu(image)
    bw = closing(image > thresh, square(3))

    # plt.imshow(bw)
    # plt.show()

    # remove artifacts connected to image border
    cleared = clear_border(bw)

    # label image regions
    label_image = label(cleared, connectivity=2)
    # to make the background transparent, pass the value of `bg_label`,
    # and leave `bg_color` as `None` and `kind` as `overlay`
    image_label_overlay = label2rgb(label_image, image=image, bg_label=0)
    for region in regionprops(label_image):
        # take regions with large enough areas
        if region.convex_area >= 10:
            # draw rectangle around segmented coins
            minr, minc, maxr, maxc = region.bbox
            pix_num = np.sum(image[minr:maxr, minc:maxc])
            rect_arr.append(ColonyRect( maxc, minc, maxr,
             minr, region.convex_area, time ))
    return rect_arr

def build_position_colonies(_path, first_frame = 9, t=32):
    _keyword = "YFP"  # "Phasefast"
    colonies_arr = None
    i=0
    j=0
    for dirpath, dirname, filenams in walk(_path):
        DIRname = dirpath.split("/")[-1]
        print(i)
        i += 1
        for _filename in sorted( filter(lambda s: ("00.tif" in s) and (_keyword in s), filenams) ):
            j += 1
            if j > t:
                break
            img_num = int(_filename.split("_")[1])
            __file = "{0}/{1}".format(dirpath, _filename)
            if first_frame < img_num < t:
                print(__file)
                cur_frame = get_rect_arr(read_image(__file, 1), img_num)
                colonies_arr = match_frames(colonies_arr, cur_frame)
                
    return colonies_arr if colonies_arr is not None else []


def run_on_all_positions(start_time,  final_time):
    all_colonies = list()
    for j in [ "" , "3"]:
        for i in range(27):
            all_colonies += build_position_colonies( f"tif{j}/Pos{i}", start_time, final_time)
    return all_colonies



def get_multi_factor(k, start_time,  final_time=32):
    colonies = run_on_all_positions(start_time,  final_time)
    
    return colonies
    
    pkt = None
    t = final_time - start_time
    pkt = np.zeros((k, t))
    for colony in colonies:
        for j, col_t in enumerate(colony):
            if j < t:
                cur_k = int(col_t.pix_num//colony[0].pix_num)
                if cur_k <= k:
                    pkt[cur_k-1][j] += 1
                else:
                    pkt = np.r_[pkt, np.zeros((cur_k-k, t))]
                    k = cur_k
                    pkt[cur_k-1][j] += 1
    return pkt



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
        if region.area >= 200:
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

def plot_regions(im):
    img = make_binary_plot(im)
    im_disp(im)
    im_disp(img)
    find_regions(img)





if __name__ == '__main__':
    # im = read_image('./img_000000080_YFPFast_000.tif', 1)
    # plot_regions(im)
    # im2 = read_image('./img_000000080_YFPFast_000.tif', 1)
    # build_position_colonies('tif-test')
    pkt = get_multi_factor(10, 10, 40)
    print(pkt)
    #
    # for k in range (5):
    #     plt.plot(pkt[k])
    # plt.show()
    # plot_colonies_graph()

