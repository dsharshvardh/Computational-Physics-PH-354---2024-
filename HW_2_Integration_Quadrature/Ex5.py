# Question (a)

# Transmission function, q(u) = sin(alpha*u)**2
# This function should have same periodicity as spacings in grating.
# Also alpha must have dimension of [L]**-1 as u has the dimension [L]. 
# So these arguments motivate a plausible choice of alpha as
#          alpha = math.pi/d where d is the spacing between slits     


# Question (b)

import math
import numpy as np


d=20*(10**-6) # slit sepration in metres
alpha = math.pi/d # relation between alpha and d


def q_1(u): #  Intensity transmission function for question 5(a) to 5(d)
    return np.sin(alpha*u)**2


# Question (c)

nslit=10
w = (nslit)*d    

# Simpon's rule for integration
def simps(q,a,b,i,x): # 'q' is the funtion we want to integrate, 'a' is lower limt, 'b' is upper limit and 'i' is the number of slices we want, 'x' is distance from centre of screen  
    lmbda = 500*(10**-9) # in metre
    f = 1 # in metre
    h = (b-a)/i
    u_lst = [a+j*h for j in range(i+1)]
    integ_real = 0
    integ_imag = 0
    for i in range(len(u_lst)):
        if i==0 or i==(len(u_lst)-1):
            integ_real = integ_real + math.sqrt(q(u_lst[i]))*math.cos(2*math.pi*x*u_lst[i]/(lmbda*f))
            integ_imag = integ_imag + math.sqrt(q(u_lst[i]))*math.sin(2*math.pi*x*u_lst[i]/(lmbda*f))
        
        elif i!=0 and i!=(len(u_lst)-1) and i%2==1:
            integ_real = integ_real + 4*math.sqrt(q(u_lst[i]))*math.cos(2*math.pi*x*u_lst[i]/(lmbda*f))
            integ_imag = integ_imag + 4*math.sqrt(q(u_lst[i]))*math.sin(2*math.pi*x*u_lst[i]/(lmbda*f))
        
        elif i!=0 and i!=(len(u_lst)-1) and i%2==0:
            integ_real = integ_real + 2*math.sqrt(q(u_lst[i]))*math.cos(2*math.pi*x*u_lst[i]/(lmbda*f))
            integ_imag = integ_imag + 2*math.sqrt(q(u_lst[i]))*math.sin(2*math.pi*x*u_lst[i]/(lmbda*f))

    I_real = (h/3)*(integ_real)
    I_imag = (h/3)*(integ_imag)
    
    return I_real**2 + I_imag**2
# Using Simpson method just because i found it easy and also it approach accuracy with less number of slices as compared to trapezoidal

x_lst = np.linspace(-0.05,0.05,1000)
Intensity = []
for i in x_lst:
    Intensity.append(simps(q_1,-w/2,w/2,1000,i))

import matplotlib.pyplot as plt

plt.plot(x_lst, Intensity, color='red', linestyle='-', linewidth=1)
plt.xlabel("x (metres)") 
plt.ylabel("Intensity")  
plt.grid() 
plt.show()  


# Quesion (d)

sc=np.zeros([len(Intensity),len(Intensity)])

for i in range(0,len(Intensity),1):
    for j in range(0,len(Intensity),1):
        sc[i,j]=Intensity[j]


# Plotting the intensity values with a perceptually uniform colormap
plt.imshow(sc, cmap='viridis', extent=[-0.5, 0.5, -0.1, 0.1], vmin=np.min(sc), vmax=np.max(sc))
plt.colorbar(label='Intensity', cmap='viridis')  
plt.xlabel("x (metres)")
plt.yticks([])
plt.show()




# Question(e) (i)


beta = alpha/2
def q_2(u): #  Intensity transmission function for question 5(e)
    return (np.sin(alpha*u)**2)*(np.sin(beta*u)**2)

Intensity_1 = []
for i in x_lst:
    Intensity_1.append(simps(q_2,-w/2,w/2,1000,i))

import matplotlib.pyplot as plt

plt.plot(x_lst, Intensity_1, color='red', linestyle='-', linewidth=1)
plt.xlabel("x (metres)") 
plt.ylabel("Intensity")  
plt.grid() 
plt.show()  



sc=np.zeros([len(Intensity_1),len(Intensity_1)])

for i in range(0,len(Intensity_1),1):
    for j in range(0,len(Intensity_1),1):
        sc[i,j]=Intensity_1[j]


# Plotting the intensity values with a perceptually uniform colormap
plt.imshow(sc, cmap='viridis', extent=[-0.5, 0.5, -0.1, 0.1], vmin=np.min(sc), vmax=np.max(sc))
plt.colorbar(label='Intensity', cmap='viridis')  
plt.xlabel("x (metres)")
plt.yticks([])
plt.show()



# Question(e) (ii)


w_1 = 90*(10**-6)

def q_3(u):
    if u >= 25*(10**-6) and u <= 45*(10**-6):
        return 1
    elif u <= -35*(10**-6) and u >= -45*(10**-6):
        return 1
    else:
        return 0
    

Intensity_2 = []
for i in x_lst:
    Intensity_2.append(simps(q_3,-w_1/2,w_1/2,1000,i))

import matplotlib.pyplot as plt

plt.plot(x_lst, Intensity_2, color='red', linestyle='-', linewidth=1)
plt.xlabel("x (metres)") 
plt.ylabel("Intensity")  
plt.grid() 
plt.show()  



sc=np.zeros([len(Intensity_2),len(Intensity_2)])

for i in range(0,len(Intensity_2),1):
    for j in range(0,len(Intensity_2),1):
        sc[i,j]=Intensity_2[j]


# Plotting the intensity values with a perceptually uniform colormap
plt.imshow(sc, cmap='viridis', extent=[-0.5, 0.5, -0.1, 0.1], vmin=np.min(sc), vmax=np.max(sc))
plt.colorbar(label='Intensity', cmap='viridis')  
plt.xlabel("x (metres)")
plt.yticks([])
plt.show()