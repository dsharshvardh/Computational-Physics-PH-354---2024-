import numpy as np
import matplotlib.pyplot as plt


pitch=np.loadtxt("pitch.txt",dtype=complex)
c_k=np.zeros(len(pitch),dtype=complex)



def fastft(lst,a):
    
    if len(lst)==1:
        return lst[0]
    else:
        lst_even=[lst[2*i] for i in range(len(lst)//2)]
        lst_odd=[lst[2*j+1] for j in range(len(lst)//2)]
        lst_even_1=fastft(lst_even,a)
        lst_odd_1=fastft(lst_odd,a)
    c_k_1 = lst_even_1 + np.exp(-(2j)*np.pi*a/len(lst))*lst_odd_1
    if len(lst)==len(pitch):
        c_k[a] = c_k_1
        c_k[a+len(pitch)//2]=lst_even_1 - np.exp(-(2j)*np.pi*a/len(lst))*lst_odd_1
    return c_k_1



c_np=np.fft.fft(pitch)

for i in range(len(pitch)//2):
    fastft(pitch,i)

plt.plot(abs(c_k),c="red")
plt.plot(abs(c_np),c="blue")
plt.ylabel("Fourier coefficient")
plt.xlabel("Index")
plt.legend(["My Code","Numpy.fft"])
plt.show()


'''
Results of my code and np.fft.fft is shown in same graph and both of them have overlapping graph thats why we are seeing only one colored plot. 
So, they are giving similiar results.
'''