#!/usr/bin/env python
import numpy as np
import time

#from ctypes import c_int
import numpy.ctypeslib as npct
# load the library, using numpy mechanisms
libssc = npct.load_library("libssc", ".")
# setup the return types and argument types
libssc.ssc_mark_start.restype = None
libssc.ssc_mark_start.argtypes = None
libssc.ssc_mark_stop.restype = None
libssc.ssc_mark_stop.argtypes = None

libssc.ssc_mark_stop()

a = np.random.rand(1000, 1000)
b = np.random.rand(1000, 1000)
t0 = time.time()

c=np.matmul(a,b)

libssc.ssc_mark_start()
c=np.matmul(c,b)
libssc.ssc_mark_stop()

t1 = time.time()
print("Walltime of the main kernel: %s sec" % (t1 - t0))
print(c.shape)

libssc.ssc_mark_start()
