# Question 2(a)

import matplotlib.pyplot as plt

alpha=1
beta=0.5
gamma=0.5
delta=2

def f_rab(x,y):
    return alpha*x-beta*(x*y)

def f_fox(x,y):
    return gamma*(x*y)-delta*y


def lotka_voltera(time,N):
    rab=[2]
    fox=[2]
    h=(time/N)
    Time = [0+i*h for i in range(N+1)]
    
    for i in range(len(Time)-1):
        k1_r = h*f_rab(rab[i],fox[i])
        k1_f = h*f_fox(rab[i],fox[i])
        
        k2_r = h*f_rab(rab[i] + k1_r/2, fox[i] + k1_f/2)
        k2_f = h*f_fox(rab[i] + k1_r/2, fox[i] + k1_f/2)
        
        k3_r = h*f_rab(rab[i] + k2_r/2, fox[i] + k2_f/2)
        k3_f = h*f_fox(rab[i] + k2_r/2, fox[i] + k2_f/2)

        k4_r = h*f_rab(rab[i] + k3_r, fox[i] + k3_f)
        k4_f = h*f_fox(rab[i] + k3_r, fox[i] + k3_f)

        r = rab[i] + (k1_r + 2*k2_r + 2*k3_r +k4_r)/6
        f = fox[i] + (k1_f + 2*k2_f + 2*k3_f +k4_f)/6
        rab.append(r)
        fox.append(f)

    return rab,fox,Time

rab,fox,time = lotka_voltera(30,10000)

plt.plot(time,rab,label="Rabbit Population",c='blue')
plt.plot(time,fox,label="Fox Population",c='red')
plt.xlabel("Time")
plt.ylabel("Population in thousands")
plt.ylim(0,8.4)
plt.legend()
plt.show()


# Question 2(b)

'''
In HW5_Fig_2, we can see oscillations in the prey-predator system. Initially, the prey population rises due to a faster reproduction rate than the rate of predation,
while the predator population declines due to a higher death rate compared to sustenance. As the prey population grows sufficiently, each predator finds more sustenance 
than the death rate, leading to a shift where predator numbers start increasing while prey numbers decline. This cycle continues until competition for prey intensifies, 
resetting the cycle.

The system's equilibrium point occurs at x=4, y=2, and starting from this condition results in a stable equilibrium. However, since our initial conditions differ, 
we observe periodic oscillations between phases dominated by prey and phases dominated by predators.
'''
        

