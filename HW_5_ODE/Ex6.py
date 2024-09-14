import numpy as np
import matplotlib.pyplot as plt
import math



# Question 6(a)

def func(r,t): # r = [x,y]
    grav_comp = []
    R = math.sqrt(r[0]**2+r[1]**2)
    GM = 6.6738*1.9891e19
    for i in range(2):
        grav_comp.append(-GM*r[i]/R**3)
    return np.array(grav_comp)


time=[0]
T=0
h=3600
rlst=[np.array([1.47e11,0])]
vlst=[np.array([0,30287])]

a=vlst[0]+0.5*h*func(rlst[0],T)

while T<2e8:
    r = rlst[len(rlst)-1]
    r_1 = r + h*a
    k = h*func(r_1,T+h)
    v_1 = a + 0.5*k
    rlst.append(r_1)
    vlst.append(v_1)
    a=a + k
    T+=h
    time.append(T)
xlst=[]
ylst=[]
for i in range(len(rlst)):
    xlst.append(rlst[i][0])
    ylst.append(rlst[i][1])

plt.plot(xlst,ylst,c='red')
plt.grid()
plt.xlabel("x")
plt.ylabel("y")
plt.show()


# Question 6(b)

m=5.9722e24
tot_en=[]
pot_en=[]
kin_en=[]

for i in range(len(rlst)):
    r=rlst[i]
    v=vlst[i]
    net_v = 0
    net_r = 0
    for i in range(len(v)):
        net_v = net_v + v[i]**2
        net_r = net_r + r[i]**2

    k_e = 0.5*m*net_v
    p_e=-((6.6738*1.9891e19)/math.sqrt(net_r))*m
    kin_en.append(k_e)
    pot_en.append(p_e)
    tot_en.append(k_e + p_e)

plt.plot(time,kin_en,label="Kinetic Energy",c='red')
plt.plot(time,pot_en,label="Potential Energy",c='blue')
plt.plot(time,tot_en,label="Total Energy",c='green')
plt.legend()
plt.xlabel("Time")
plt.ylabel("Energy")
plt.grid()
plt.show()


# Question 6(c)

time_new=[]
for i in range(len(time)): 
    time_new.append(time[i]/31536000)   # Converting from seconds to years
plt.plot(time_new,tot_en,c='red')
plt.title("Total Energy vs Time")
plt.xlabel("Time(in years)")
plt.ylabel("Energy")
plt.show()

