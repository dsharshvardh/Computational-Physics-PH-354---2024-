# Explaination of Question (a) and (c) is given in image named as HW_2_Ans_Q3.a),c)_Q4.a)

# Question (b)

import numpy as np

def gaussxw(N):

    # Initial approximation to roots of the Legendre polynomial
    a = np.linspace(3,4*N-1,N)/(4*N+2)
    x = np.cos(np.pi*a+1/(8*N*N*np.tan(a)))

    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta>epsilon:
        p0 = np.ones(N,float)
        p1 = np.copy(x)
        for k in range(1,N):
            p0,p1 = p1,((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))

    # Calculate the weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)

    return x,w


def gaussxwab(N,a,b):
    x,w = gaussxw(N)
    return 0.5*(b-a)*x+0.5*(b+a),0.5*(b-a)*w   


def f(m,x,a):
        return 4/(np.sqrt((2/m)*(a**4-x**4)))

# Finding Time period
N=20
x_1,w_1=gaussxwab(N,0,2)

t=0  # starting time
for i in range(N):
    t = t + w_1[i]*f(1,x_1[i],2)

print("Time period is =",t)



import matplotlib.pyplot as plt
def time_period(b):
    N=20
    a=0
    x_2,w_2=gaussxwab(N,a,b)
    t=0

    for i in range(N):
        t = t+ w_2[i]*f(1,x_2[i],b)
    return t

amp_lst = []
n=200
h = (2-0.1)/n
for i in range(n+1):
    amp_lst.append(0.1+i*h)


time=[]
for i in amp_lst:
    time.append(time_period(i))



plt.scatter(amp_lst, time, marker='o', color='red') 
plt.plot(amp_lst, time, color='black', linestyle='-', linewidth=1)
plt.title("Anharmonic Oscillator")  
plt.xlabel("Amplitude") 
plt.ylabel("Time period")  
plt.grid() 
plt.savefig("HW_2_Fig_3")  
plt.show()  
