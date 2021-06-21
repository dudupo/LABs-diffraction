

from os import walk, system
from copy import deepcopy
def rename( _path = "."):
	for dirpath, dirname, filenams in walk(_path): 
		for filename in filter( lambda s : "pkl" in s , filenams):   
			print(filename)
			q = deepcopy(filename).replace(":", "-").replace(" ", "_").replace("\'", "")
			system( f"mv '{filename}' ./{q}")

        #     DIRname = dirpath.split("/")[-1]

if __name__ == "__main__":
	rename()
