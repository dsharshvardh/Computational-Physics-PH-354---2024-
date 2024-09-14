
import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft

# Question 3(a)
dow_data=np.loadtxt("dow.txt")
plt.plot(dow_data,color='blue')
plt.xlabel("Days", fontname='Times New Roman',fontsize=18)
plt.ylabel("Dow Jones", fontname='Times New Roman',fontsize=18)
plt.show()



# Question 3(b)
def num_ck(percent,lst): # you tell the percent it will find number of coefficients from list
    return (len(lst)*percent)//100
c_k = fft.rfft(dow_data)
c_k_10_percnt=np.zeros(len(c_k),dtype=complex)



# Question 3(c)
limit_10 = num_ck(10,c_k)
for i in range(limit_10):
    c_k_10_percnt[i]=c_k[i]



# Question 3(d)
dow_data_10=fft.irfft(c_k_10_percnt)
plt.plot(dow_data,c="blue")
plt.plot(dow_data_10,c="red")
plt.legend(["Dow Data","First 10% Fourier Coefficients"])
plt.show()


'''
When we are using only first 10 percent coefficients, the plot is doing enough justice to the origial data. This tell us that first 10 percent coefficients represent frequencies
which has major contribution to the function representing the data and we get good graph by using first 10 percent coefficients to represent actual data.
'''



# Question 3(e)
c_k_2_percnt = np.zeros(len(c_k),dtype=complex)
limit_2=num_ck(2,c_k)

for i in range(limit_2):
    c_k_2_percnt[i]=c_k[i]

dow_data_50=fft.irfft(c_k_2_percnt)
plt.plot(dow_data,c="blue")
plt.plot(dow_data_50,c="red")
plt.legend(["Dow Data","First 2% Fourier Coefficients"])
plt.show()


'''
We can see that with the first 2% coefficients, it is not giving that much accurate graph with data points.
 Therefore, we are neglecting some frequencies which are important in the representation of the data.
'''


# Together showing all three

plt.plot(dow_data,c="blue")
plt.plot(dow_data_10,c="red")
plt.plot(dow_data_50,c="green")
plt.legend(["Dow Data","First 10% Fourier Coefficients","First 2% Fourier Coefficients"])
plt.show()



# Square wave function smoothing using fourier transform

lst=[1 if i<500 else -1 for i in range(1000)]
lst=np.array(lst)
c_k_smooth=fft.rfft(lst)
c_k_smooth_new=np.zeros(len(c_k_smooth),dtype=complex)
for i in range(10):
    c_k_smooth_new[i]=c_k_smooth[i]
inv_ft=fft.irfft(c_k_smooth_new)
plt.plot(lst,c="blue")
plt.plot(inv_ft,c="red")
plt.xlabel("Number of sampled points(1ms)", fontname='Times New Roman',fontsize=18)
plt.legend(["Original Data","First 10 Fourier Coefficients"])
plt.show()


'''
We see that when we use first 10 fourier coefficients in the square wave signal it has the basic structure of the square wave, but near the
discontinuity(including the end points), there is a high amplitude of oscillations of the fourier terms. This happens because to represent a sudden jump(from +1 to -1 or
vice versa), we need all the fourier terms(each of which has a continuous form) giving a substantial contribution. 
Since the entire process of Fourier Transformation happens assuming periodic functions, this also means that there will be large discontinuities at the end points where
the function tries to force the waveform to be periodic since the function is a propagating square wave train with inequal values at it's start and end points.
'''