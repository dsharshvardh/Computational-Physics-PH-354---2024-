# Question 6

import numpy as np
import matplotlib.pyplot as plt


l=5e-7
a=2e-4
b=10*a
f=1
N=1000


def func(n):
    y=n*b/N-a/2
    if y<=1e-4:
        return abs(np.sin(np.pi*y/20e-6))
    else:
        return 0


y=np.zeros(int(N))
for i in range(int(N)):
    y[i]=func(i)

c_k=np.fft.fft(y)

Intensity=np.zeros(400)
for i in range(200):
    Intensity[i]=abs(c_k[N+i-200])**2
    Intensity[i+200]=abs(c_k[i])**2
Intensity=Intensity*(b/N)**2
plt.plot(Intensity,c="red")
plt.ylabel("Intensity")
plt.xlabel("x(cm)")
plt.show()





