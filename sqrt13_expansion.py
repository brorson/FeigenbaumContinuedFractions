import gmpy2
from gmpy2 import mpfr
import numpy as np
import math

import matplotlib.pyplot as plt

gmpy2.get_context().precision=1700

delta = gmpy2.sqrt(13.0)

N = 500

a = np.zeros((N,), dtype=object)

a[0] = np.floor(delta)
xi = delta-a[0]
x = mpfr("1.0")/xi

for i in range(1,N):
    a[i] = np.floor(x)
    xi = x - a[i]
    x = mpfr("1.0")/xi

y = np.zeros(len(a), dtype=object)
for i in range(len(a)):
    y[i] = a[i]
    print(y[i])

plt.plot(y, 'r.')
plt.yscale('log')
plt.show()

    
dcf = a[N-1]
for i in range((N-2), -1, -1):
    dcf = a[i] + mpfr("1.0")/dcf

print("-------------------------------------")    
print(delta)
print("-------------------------------------")
print(dcf)
print("-------------------------------------")

K = np.power(np.prod(a),(mpfr("1.0")/N)) 
print(K)

