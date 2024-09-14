# Question (a) The derivation asked in this question is shared in a photo named as HW_2_Ans_Q3.a),c)_Q4.a)


# Question (b)

import numpy as np
import matplotlib.pyplot as plt

G=6.674*1e-11
sigma=100


z=np.linspace(0,10,400)


def gaussxw(N):

    a = np.linspace(3,4*N-1,N)/(4*N+2)
    x = np.cos(np.pi*a+1/(8*N*N*np.tan(a)))

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

    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)
    return x,w

def gaussxwab(N,a,b):
    x,w = gaussxw(N)
    return 0.5*(b-a)*x+0.5*(b+a),0.5*(b-a)*w



def f(x,y,z):
    return ((G*sigma*z)/(x**2+y**2+z**2)**1.5)

def force(z):
    N=100
    x,w=gaussxwab(N,-10/2,10/2)
    intgrtn=0
    for i in range(N):
        for j in range(N):
           intgrtn = intgrtn + w[i]*(w[j]*f(x[i],x[j],z))   # Here j loop is acting like integration over one of the variable 
                                                            # (doesn't matter whichever you think) for each value of the other variable which took care by i loop.  
    return intgrtn

frc=[]
for i in (z):
    frc.append(force(i))


plt.scatter(z, frc, marker='D',s=50, color='red') 
plt.plot(z, frc, color='black', linestyle='-', linewidth=1)
plt.title("Gravitational pull")  
plt.xlabel("z") 
plt.ylabel("Force")  
plt.grid()  
plt.show()  



# Question (c)

# Fz has discontinuity at z=0, Fz increases for smaller values of z but then suddenly drops to zero telling
# that force Fz has falling trend for small vaues of z which is wrong. A way to solve this is by decreasing 
# the spacing between consecutive z values. If the spacing is small than the desired z resolution, this numercal 
# artifact can be avoided upto desired z resolution. However, accurate integral for smaller z values require higher 
# number of Quadrature points. 
# I am attacing added a photo of the plot with higher number of sampling points = 800 to get more accurate curve to show what i just said.