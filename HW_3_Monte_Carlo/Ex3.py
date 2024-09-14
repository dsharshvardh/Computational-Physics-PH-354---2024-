
# Question 3(a)

import numpy as np
import math


def ising_energy(J,lattice_spin):
    E=0
    E = E + np.sum(np.multiply(np.roll(lattice_spin,1,axis=0),lattice_spin))
    E+= np.sum(np.multiply(np.roll(lattice_spin, 1, axis=1), lattice_spin))
    Energy=(-J)*E
    return Energy

lattice_spin=[]
h_rng=1000
for i in range(20):
    lst_row=[]
    for j in range(20):
        p=np.random.random()
        if p<0.5:
            lst_row.append(-1)
        else:
            lst_row.append(+1)
    lattice_spin.append(lst_row)






# Question 3(b),(c) combined

def ising_magnetisation(J,steps,T,lattice_spin):
    time=[]   # For steps
    t=0
    magnetisation=[] # Will store magentisation on each step
    Energy=[] # Will store energy on each step

    for i in range(steps):
        mag=0 # for magnetisation of lattice on a particular step

        mag=np.sum(lattice_spin)
        magnetisation.append(mag)
        time.append(t)
        t=t+1
        
        old_energy=ising_energy(J,lattice_spin)
        Energy.append(old_energy)  
        spin_switch=[np.random.randint(0,20),np.random.randint(0,20)]  # Randomly choosing which lattice site spin to switch
        
        lattice_spin[spin_switch[0]][spin_switch[1]]=lattice_spin[spin_switch[0]][spin_switch[1]]*(-1)
        new_energy=ising_energy(1,lattice_spin)
        energy_change = new_energy - old_energy  # Lattice after and before spin switch
        if energy_change <= 0 or np.random.rand()<math.exp(-(new_energy-old_energy)/T):
                pass
        else:
            lattice_spin[spin_switch[0]][spin_switch[1]]=lattice_spin[spin_switch[0]][spin_switch[1]]*(-1)
    
    return magnetisation, Energy, time, lattice_spin

magnetisation,Energy,time,last_lattice_spin=ising_magnetisation(1,1000000,1,lattice_spin)           

'''
Please wait, it will take some time  :)
'''



import matplotlib.pyplot as plt

# Plotting Magnetization with a specific color (e.g., blue)
plt.plot(time, magnetisation, label='Magnetization', color='blue')
plt.xlabel('Time')
plt.ylabel('Total Magnetization')
plt.title('Magnetization vs Time for T=1')
plt.legend()
plt.show()

# Plotting Energy with a different color (e.g., red)
plt.plot(time, Energy, label='Energy', color='red')
plt.xlabel('Time')
plt.ylabel('Total Energy')
plt.title('Energy vs Time for T=1')
plt.legend()
plt.show()



# Question 3(e)

plt.imshow(last_lattice_spin,vmin=-1,vmax=1)
plt.title(f'Final orientation of spins T={1}')
plt.axis('off')
plt.show()





# Question 3(d)


# Part D
# At T=1 we see that lattice want to stay in ferromagnetic state and after many steps the magnetisation attained the value of 400 (total spin) which is maximum possible magnetisation 
# a 20 x 20 lattice can have although the sign of the magnetisation is random (+ or -) and it depends on direction of spins in initial steps when we randomly switch the spins. So we can get magnetisation
# in both the direction. And value will saturate to 400.




'''
# Question 3(e)

At T=1, we see a saturation of the magnetisation, consistent with the ferromagnetic phase.

At T=2 as well we see a similar saturation near the maximum possible magnetisation(400), but we see some oscilliations, due to the higher thermal energy permitting
excitations from the ground state.

At higher T=3, we start seeing oscillations near the magnetisation = 0, without any saturation to 400. That means we get the paramagnetic phase, with
a net 0 magnetisation.

For a 2D Ising model, we know that the critical temperature is almost equal to 2.269*J/k, so for T>2.269, we start seeing the paramagnetic phase.

From the plots, we see the formation of large domains in the ferromagnetic phase having same spins. These contribute to a net magnetisation in a given direction. 
In the T=2 plot, due to the thermal fluctuations, we see smaller domains with spins aligned in the opposite directions which reduce the total magnetisation from the maximum magnetisation.

In the paramagnetic phase however, we see multiple such domains with approximately similar numbers of spins opposite either way. The domains grow smaller with temperature
and at high enough temperature, we expect to see a random arrangement of spins.

The behaviour of the system is always such that it minimises it's free energy, which consists of the energy as well as entropy. At low temperatures, we consider only the energy term,
which is minimised if all the spins are aligned. Meanwhile, as T increases, we get contributions from the entropy which tends to randomise the system, resulting in a randomly oriented system
which has spins not aligned with each other.

'''