
import numpy as np
'''
d2x/dt2 = -x
if y =dx/dt
dy/dt = -x
'''
import matplotlib.pyplot as plt



# Question 5(a) and 5(b)


def harmonic(x,y):
    w=1
    a=y
    b=-w**2*x
    return a,b



def func(N,amplitude,Time,f): # time is the highest time
    xlst=[amplitude]
    ylst=[0]
    time=[0]
    # N=100000
    h=Time/N
    h1=50/N

    for i in range(N):
        x1=xlst[i]
        y1 = ylst[i]
        t=time[i]
        k1x,k1y=f(x1,y1)
        k2x,k2y=f(x1+ h1 *k1x/2,y1+ h1 *k1y/2)
        k3x,k3y=f(x1+ h1 *k2x/2,y1+ h1 *k2y/2)
        k4x, k4y = f(x1 +  h1 *k3x, y1 + h1 *k3y)
        time.append(t+h)
        xlst.append(x1+ h1 *(k1x+2*(k2x+k3x)+k4x)/6)
        ylst.append(y1 + h1 * (k1y + 2 * (k2y + k3y) + k4y) / 6)

    return xlst,ylst,time



xlst_1,ylst_1,time_1 = func(100000,1,50,harmonic)
xlst_2,ylst_2,time_2 = func(100000,2,50,harmonic)
plt.xlabel("Time")
plt.ylabel("x")
plt.title("Harmonic Oscillator")
plt.plot(time_1,xlst_1,label="Amplitude = 1",c='red')
plt.plot(time_2,xlst_2,label="Amplitude = 2",c='blue')
plt.legend()
plt.show()


'''
We can see in the plot that graph of both the amplitudes 1 and 2 are overlapping with different amplitude.
'''


# Question 5(c)

# Anharmonic Oscillator
def anharmonic(x,y):
    a=y
    b=-x**3
    return a,b


xlst_an_1,ylst_an_1,time_an_1 = func(100000,1,50,anharmonic)
xlst_an_2,ylst_an_2,time_an_2 = func(100000,2,50,anharmonic)
xlst_an_3,ylst_an_3,time_an_3 = func(100000,0.5,50,anharmonic)

plt.xlabel("Time")
plt.ylabel("x")
plt.title("Anharmonic Oscillator")
plt.plot(time_an_3,xlst_an_3,label="Amplitude = 0.5",c='green')
plt.plot(time_an_1,xlst_an_1,label="Amplitude = 1",c='red')
plt.plot(time_an_2,xlst_an_2,label="Amplitude = 2",c='blue')
plt.legend()
plt.show()




# Question 5(d)
plt.xlabel("x")
plt.ylabel("dx/dt")
plt.title("Phase space")
plt.plot(xlst_1,ylst_1,label="Harmonic Oscillator",c='red')
plt.plot(xlst_an_1,ylst_an_1,label="Anharmonic Oscillator",c='blue')
plt.legend()
plt.show() 





# Question 5(e)

# Van Der Pol Oscillator
def van_der(x,y,mu):
    w=1
    a=y
    b=mu*(1-x**2)*y-w**2*x
    return a,b


def func_2(N,amplitude,mu,f):
    xlst=[amplitude]
    ylst=[0]
    time=[0]    
    h1=20/N

    for i in range(N):
        x1=xlst[i]
        y1 = ylst[i]
        t=time[i]
        k1x,k1y=f(x1,y1,mu)
        k2x,k2y=f(x1+ h1 *k1x/2,y1+ h1 *k1y/2,mu)
        k3x,k3y=f(x1+ h1 *k2x/2,y1+ h1 *k2y/2,mu)
        k4x, k4y = f(x1 +  h1 *k3x, y1 + h1 *k3y,mu)
        time.append(t+h1)
        xlst.append(x1+ h1 *(k1x+2*(k2x+k3x)+k4x)/6)
        ylst.append(y1 + h1 * (k1y + 2 * (k2y + k3y) + k4y) / 6)
    return xlst,ylst,time

xlst_van_mu1, ylst_van_mu1, time = func_2(100000,1,1,van_der)
xlst_van_mu2, ylst_van_mu2, time = func_2(100000,1,2,van_der)
xlst_van_mu3, ylst_van_mu3, time = func_2(100000,1,4,van_der)


plt.title("Van-Der-Pol Oscillator")
plt.xlabel("x")
plt.ylabel("dx/dt")
plt.title("Phase space")
plt.plot(xlst_van_mu1,ylst_van_mu1,label="\u03BC=1",c='green')
plt.plot(xlst_van_mu2,ylst_van_mu2,label="\u03BC=2",c='blue')
plt.plot(xlst_van_mu3,ylst_van_mu3,label="\u03BC=4",c='red')
plt.grid()
plt.legend()
plt.show()






    


