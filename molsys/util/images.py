# -*- coding: utf-8 -*-


import numpy as np

class NumpyDict(dict):
    """Allow for numpy array keys"""
    def __getitem__(self, key):
        return dict.__getitem__(self, tuple(key))

images = []
for x in range(-1,2):
    for y in range(-1,2):
        for z in range(-1,2):
            images.append([x,y,z])
images = np.array(images,"d")

imgidx, revimgidx = np.where((images == -images[:,np.newaxis]).all(axis=-1))

# INDEX TO ARRAY
idx2arr = dict(zip(imgidx,images))
# ARRAY TO INDEX
imgtuple = [tuple(i) for i in images.tolist()]
arr2idx_items = zip(imgtuple,imgidx)
arr2idx = NumpyDict(arr2idx_items)
