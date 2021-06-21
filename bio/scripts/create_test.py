import cv2 
from os import walk
import numpy as np
def loadf(_path):
    print(_path)
    return cv2.imread(f'{_path}', cv2.IMREAD_GRAYSCALE)


def create_test():
    _path = "./tif/Pos2"

    _keyword = "YFP" #"Phasefast"

    # _dict = {}
    for dirpath, dirname, filenams in walk(_path): 
        DIRname = dirpath.split("/")[-1]
        # _dict[DIRname] = [  ]
        for _filename in filter( lambda s : ("tif" in s ) and ( _keyword in s), filenams):               
            __file = "{0}/{1}".format( dirpath, _filename )
            print( loadf(__file)[ 200:300, 200:300])
            cv2.imwrite(  f"test/{_filename} "  , loadf(__file)[ 500:1000, 500:1000])
            print(__file)





if __name__ == "__main__":
    create_test()