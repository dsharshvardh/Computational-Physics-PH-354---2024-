import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft
import dcst

def num_ck(percent,lst): # you tell the percent it will find number of coefficients from list
    return (len(lst)*percent)//100


# Question 4(a)

dow2_data=np.loadtxt("dow2.txt")
plt.plot(dow2_data,c="blue")
plt.xlabel("Days",fontname='Times New Roman',fontsize=18)
plt.ylabel("Dow_2",fontname='Times New Roman',fontsize=18)
plt.show()
c_k=fft.rfft(dow2_data)
c_k_2_percent=np.zeros(len(c_k),dtype=complex)
limit_2=num_ck(2,c_k)
for i in range(limit_2):
    c_k_2_percent[i]=c_k[i]
dow2_data_2=fft.irfft(c_k_2_percent)
plt.plot(dow2_data,c="blue")
plt.plot(dow2_data_2,c="red",linewidth=2)
plt.legend(["Dow_2 Data","First 2% Fourier Coefficients"])
plt.show()




# Question 4(b)

def dct(y1):
    N = len(y1)
    y2 = np.empty(2*N,float)
    y2[:N] = y1[:]
    y2[N:] = y1[::-1]
    c = fft.rfft(y2)
    phi = np.exp(-1j*np.pi*np.arange(N)/(2*N))
    return np.real(phi*c[:N])

def idct(a):
    N = len(a)
    c = np.empty(N+1,complex)
    phi = np.exp(1j*np.pi*np.arange(N)/(2*N))
    c[:N] = phi*a
    c[N] = 0.0
    return fft.irfft(c)[:N]


c_k=dct(dow2_data)
c_k_1=np.zeros(len(c_k),dtype=complex)
limit_again=num_ck(2,c_k)
for i in range(limit_again):
    c_k_1[i]=c_k[i]
dow2_data_3=idct(c_k_1)
plt.plot(dow2_data,c="blue")
plt.plot(dow2_data_3,c="red")
plt.legend(["Dow Data","First 2% Cosine Transform Coefficients"])
plt.show()