#!/usr/bin/env python
import numpy as np
import time

# somehow the newest sde breaks with keras
#import pandas as pd
#from keras.models import Model
#from keras.layers import Dense, Input

import numpy.ctypeslib as npct
from ctypes import c_int
# load the library, using numpy mechanisms
libssc = npct.load_library("libssc", ".")
# setup the return types and argument types
libssc.ssc_mark_start.restype = c_int
libssc.ssc_mark_start.argtypes = [c_int, c_int]
libssc.ssc_mark_stop.restype = c_int
libssc.ssc_mark_stop.argtypes = [c_int, c_int]

_ = libssc.ssc_mark_stop(0,0)

a = np.random.rand(1000, 1000)
b = np.random.rand(1000, 1000)
t0 = time.time()

c = np.matmul(a,b)

_ = libssc.ssc_mark_start(1,0)
c = np.matmul(c,b)
_ = libssc.ssc_mark_stop(1,0)

t1 = time.time()
print("Walltime of the main kernel: %s sec" % (t1 - t0))
print(c.shape)

_ = libssc.ssc_mark_start(0,0)
