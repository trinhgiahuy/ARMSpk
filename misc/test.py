import numpy as np
a = np.random.rand(1000, 1000)
b = np.random.rand(1000, 10)
c=np.matmul(a,b)
print(c.shape)
