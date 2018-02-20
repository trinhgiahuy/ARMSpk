#!/usr/bin/env python
import numpy as np
import time
#a = np.random.rand(4000, 4000)
#b = np.random.rand(4000, 3800)
t0 = time.time()
#c=np.matmul(a,b)
t1 = time.time()
print("Walltime of the main kernel: %s sec" % (t1 - t0))
#print(c.shape)
