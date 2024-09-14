import numpy as np
import matplotlib.pyplot as plt



# Question 7(a)

data=np.loadtxt("blur.txt")
N=len(data)
print(N)
plt.imshow(data,cmap="gray")
plt.show()




# Question 7(b)

sgma=25
x=len(data)
y=len(data[0])

def func(x,y):
    return np.exp(-(x**2+y**2)/(2*sgma**2))

lst=[[func((j+x/2)%x-x/ 2,(i+y/2)%y-y/2) for i in range(x)]for j in range(y)]

plt.imshow(lst,cmap="gray")
plt.show()





# Question 7(c)

blurr_photo=np.fft.rfft2(data)
point_spread=np.fft.rfft2(lst)

unblurr_photo=np.zeros([y,x//2+1],complex)
epsilon = 1e-3
for i in range(x//2+1):
    for j in range(y):
        if abs(point_spread[j][i])<epsilon:
            unblurr_photo[j][i] = blurr_photo[j][i]
        else:
            unblurr_photo[j][i]=blurr_photo[j][i]/(point_spread[j][i])

plt.imshow(np.fft.irfft2(unblurr_photo),cmap="gist_gray")
plt.show()



# Question 7(d)

'''
We are dividing with a small number at points when the Fourier coefficients of the point spread function are very small which causes machine precision errors on a computer
and thus, we cannot completely recover the image despite knowing the point spread function whenever we come across poles in the Fourier coefficients.
'''