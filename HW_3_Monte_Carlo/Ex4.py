# Question 4(a)


import numpy as np
import math
import matplotlib.pyplot as plt

def f(x):
    return x**2 - math.cos(4*math.pi*x)

def annealing_min(f,x,T,Tmin): # Here f(x) is the function,  x is starting point of code, T = initial temprature, Tmin is last cooled temprature where we want to break the itteration    
    xlst=[]
    Tmax=T
    time=0
    while T>Tmin:
        time=time+1
        x_new = x + np.random.normal(0, 1) # to take a random number with gaussian distribution having locus 0 and standard deviration 1
        E_0=f(x)
        E_n=f(x_new)

        if E_n<=E_0:
            xlst.append(x_new)
            x=x_new

        else:
            if np.random.random() < math.exp(-(E_n-E_0)/T):
                xlst.append(x_new)
                x=x_new
            else:               
                xlst.append(x)

        T=Tmax*np.exp(-time/1e4) # exponential cooling schedule
    return xlst

xlst=(annealing_min(f,2,100,1e-5))



plt.plot(xlst, marker='.', markersize=3, color='red', linestyle='None')
plt.xlabel('Time')
plt.ylabel('x values')
plt.title('For minima of $f(x) = x^2 - \cos(4\pi x)$')
plt.show()

print('minima of this function is:',xlst[-1])






# Question 4(b)

def f_2(x):
    return math.cos(x) + math.cos(math.sqrt(2)*x) + math.cos(math.sqrt(3)*x)


def annealing_min_2(f,x,T,Tmin): 
    xlst=[]
    Tmax=T
    time=0
    while T>Tmin:
        time=time+1
        x_new = x + np.random.normal(0, 1) 
        
        while x_new>50 or x_new<0:
            x_new=x + np.random.normal(0,1)

        E_0=f(x)
        E_n=f(x_new)

        if E_n<=E_0:
            xlst.append(x_new)
            x=x_new

        else:
            if np.random.random() < math.exp(-(E_n-E_0)/T):
                xlst.append(x_new)
                x=x_new
            else:               
                xlst.append(x)

        T=Tmax*np.exp(-time/1e5)  # exponential cooling schedule
    return xlst


x_2_lst=(annealing_min_2(f_2,30,100,1e-5))



plt.plot(x_2_lst, marker='.', markersize=3, color='red', linestyle='None')
plt.xlabel('Time')
plt.ylabel('x values')
plt.title('For minima of $f(x) = \cos(x) + \cos(\sqrt{2}x) + \cos(\sqrt{3}x)$')
plt.show()

print('minima of this function is:',x_2_lst[-1])


# When we run this code in Question 2(b) for f_2 function, then it is not only showing the minimas which are given in the question but it was showing other also like x=34
# So to confirm I checked the plot of this equation and I found out that 34 is also a minima and it also has some significant depth. 