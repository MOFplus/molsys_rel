# -*- coding: utf-8 -*-


import numpy as np

images = []
for x in range(-1,2):
    for y in range(-1,2):
        for z in range(-1,2):
            images.append([x,y,z])
images = np.array(images,"d")
