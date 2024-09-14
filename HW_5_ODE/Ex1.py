# Question 1(a)


import matplotlib.pyplot as plt

def f(Vout,t,RC):
    
    if int(2*t)%2==0:
        return (1-Vout)/RC
    else:
        return (-1-Vout)/RC



def vout(time,N,RC): # N = number of grid, time = end time
    h=(time/N)
    Time = [0+i*h for i in range(N+1)]
    Vout=[]
    v=0
    for i in range(len(Time)):
       
       k1=h*f(v,Time[i],RC)
       k2=h*f(v+k1/2,Time[i]+h/2,RC)
       k3=h*f(v+k2/2,Time[i],RC)
       k4=h*f(v+k3,Time[i]+h,RC)
       Vout.append(v+(k1+2*(k2+k3)+k4)/6)
       v=Vout[i]
    
    return Vout


Vout_1 = vout(10,1000000,1)

plt.plot(Vout_1,c='red')
plt.title("RC=1")
plt.xlabel("Time(s)")
plt.ylabel("Vout(V)")
plt.show()


Vout_2 = vout(10,1000000,0.1)

plt.plot(Vout_2,c='red')
plt.title("RC=0.1")
plt.xlabel("Time(s)")
plt.ylabel("Vout(V)")
plt.show()


Vout_3 = vout(10,1000000,0.01)

plt.plot(Vout_3,c='red')
plt.title("RC=0.01")
plt.xlabel("Time(s)")
plt.ylabel("Vout(V)")
plt.show()




'''
We get less accurate plots for lower values of h. We start getting nice continuous plots from h=0.0001. We will take h=1e-5 for this case.

When adjusting the value of RC, which dictates the time constant and controls the rate of charge/discharge, we observe a transition when Vin = +1, Vout begins to increase until it 
reaches either Vout = 1 or Vin switches to -1. Subsequently, the potential undergoes a decay, unless Vout = 1 is maintained.
This behavior leads to a cyclic pattern of charging and discharging due to the oscilliating applied potential. 
The time taken for Vout to approach 1 is determined by RC(following V as exp(-t/RC)). A larger RC extends this time, making it longer than 0.5 (the switching time of Vin),
resulting in Vout decreasing during the discharging phase. This causes an overall decay in oscillation, 
deviating from the uniform charging and discharging cycles due to the non-uniform initial conditions at each stage.

As a result, for RC = 1, we observe a decaying envelope, while for RC = 0.1 and RC = 0.01, we witness charging and discharging cycles with peaks nearing +1 or -1.

'''