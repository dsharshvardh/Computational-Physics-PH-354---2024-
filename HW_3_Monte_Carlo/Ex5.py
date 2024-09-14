
import numpy as np
import os
import random
import math
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


# Function for energy, In the lattice where will be an atom of a dimer, we will insert 1 in the list and when there is no atom, we will insert 0

def Energy(lattice):
    sum=0
    for i in range(len(lattice)):
        for j in range(len(lattice[i])):
            if lattice[i][j]>0:
                sum=sum+1
            else:
                continue
    Energy=-sum/2
    return Energy



def random_dimer(L): # for randomly putting a dimer in a lattice size of L

    x=random.randint(0,L-1) # x and y position of first atom of dimer
    y=random.randint(0,L-1)


    # We will have to decide neighbour of first atom depending on where was first atom - bulk, corner or edges
    if x==0 and y==0: # For corners
        toss=random.randint(0,1)
        if toss==0:
            x_1=x+1
            y_1=y
        else:
            x_1=x
            y_1=y+1

    elif x==L-1 and y==0: # For corners
        toss=random.randint(0,1)
        if toss==0:
            x_1=x-1
            y_1=y
        else:
            x_1=x
            y_1=y+1

    elif x==0 and y==L-1: # For corners
        toss=random.randint(0,1)
        if toss==0:
            x_1=x+1
            y_1=y
        else:
            x_1=x
            y_1=y-1

    elif x==L-1 and y==L-1: # For corners
        toss=random.randint(0,1)
        if toss==0:
            x_1=x-1
            y_1=y
        else:
            x_1=x
            y_1=y-1

    elif x==0 and (y!=0 and y!=L-1): # For edges
        toss=random.randint(0,2)
        if toss==0:
            x_1=x+1
            y_1=y
        elif toss==1:
            x_1=x
            y_1=y+1
        else:
            x_1=x
            y_1=y-1

    elif x==L-1 and (y!=0 and y!=L-1):  # For edges
        toss=random.randint(0,2)
        if toss==0:
            x_1=x-1
            y_1=y
        elif toss==1:
            x_1=x
            y_1=y+1
        else:
            x_1=x
            y_1=y-1

    elif y==0 and (x!=0 and x!=L-1):  # For edges
        toss=random.randint(0,2)
        if toss==0:
            x_1=x+1
            y_1=y
        elif toss==1:
            x_1=x-1
            y_1=y
        else:
            x_1=x
            y_1=y+1
    
    elif y==L-1 and (x!=0 and x!=L-1):  # For edges
        toss=random.randint(0,2)
        if toss==0:
            x_1=x+1
            y_1=y
        elif toss==1:
            x_1=x-1
            y_1=y
        else:
            x_1=x
            y_1=y-1
    
    else: # For Bulk
        toss=random.randint(0,3)
        if toss==0:
            x_1=x+1
            y_1=y
        elif toss==1:
            x_1=x-1
            y_1=y
        elif toss==2:
            x_1=x
            y_1=y+1
        else:
            x_1=x
            y_1=y-1

   
    return x,y,x_1,y_1
        



# Question 5(a)


def annealing_a(L,steps): # L is size of the lattice      This function is without any cooling
     
    lattice = [[0 for i in range(L)]for j in range(L)]

    
    for i in range(steps):
        value = np.max(lattice)
        x,y,x_1,y_1=random_dimer(L)
        if lattice[x][y]==0 and lattice[x_1][y_1]==0:
            
            lattice[x][y] = value+1
            lattice[x_1][y_1] = value+1
            
        
        elif lattice[x][y]!=0 and lattice[x_1][y_1]!=0 and lattice[x][y]==lattice[x_1][y_1]:
            lattice[x][y]=0
            lattice[x_1][y_1]=0
            
    return lattice

lattice_dimer_a=annealing_a(50,10000) 
no_dimer = np.ma.masked_equal(lattice_dimer_a, 0)
cmap = plt.cm.viridis
cmap.set_bad(alpha=0)  
plt.imshow(no_dimer, cmap=cmap, interpolation='none', norm=Normalize(vmin=no_dimer.min(), vmax=no_dimer.max()))
lattice_energy=Energy(lattice_dimer_a)
plt.title(f'tau = {500} and Energy = {lattice_energy}')
cbar = plt.colorbar(label='Dark color = oldest dimer and Light color = youngest dimer', orientation='vertical', fraction=0.046, pad=0.1)
cbar.set_ticks([])

plt.show()






# Question 5(b)

'''
This function includes cooling. With cooling, T is decreasing so probability of removal of dimer is decreasing, so energy will decrease faster with cooling.
'''

def annealing(L,tau): # L is size of the lattice and tau is the step size
    
    T=1
    t=0
    Tmax=T
    Tmin=1e-1
    lattice = [[0 for i in range(L)]for j in range(L)]

    
    while   t < tau:  #T>Tmin    # This is for seeing the snapshots for animation of dimer arrangement. If you want to see how speed of cooling affects,
                       # we can run the complete process for long period of time to see that how tau(cooling) affects  dimer arrangment. To do it 
                       # comment on T > Tmin condition and comment off t < tau condition and also comment on the commented off code in the end for 
                       # plotting the last picture of dimers after many steps (till T<Tmin) and comment off code for snapshots which is commented 
                       # on now and then you will see that even if we try to put dimers for a many steps (till T < 10e-1), still we get lesser dimers at tau=500 then tau 10000.  
        t=t+1
        
        value = np.max(lattice)
        x,y,x_1,y_1=random_dimer(L)
        if lattice[x][y]==0 and lattice[x_1][y_1]==0:
            
            lattice[x][y] = value+1
            lattice[x_1][y_1] = value+1
            
        
        elif lattice[x][y]!=0 and lattice[x_1][y_1]!=0 and lattice[x][y]==lattice[x_1][y_1]:
            
            new_lattice = np.copy(lattice[x][y])
            lattice[x][y]=0
            lattice[x_1][y_1]=0
            if np.random.uniform(0,1)< np.exp(-1/T):
                l=1
                
                
            else:
                lattice[x][y] = new_lattice
                lattice[x_1][y_1] = new_lattice
                
                
                        
        T=Tmax*math.exp(-t/tau) 

    
    return lattice



'''
To check the effect of tau put 500 also in place of 10000 to see the effect of tau.
'''

# lattice_dimer=annealing(50,10000) 
# no_dimer = np.ma.masked_equal(lattice_dimer, 0)
# cmap = plt.cm.viridis
# cmap.set_bad(alpha=0)  
# plt.imshow(no_dimer, cmap=cmap, interpolation='none', norm=Normalize(vmin=no_dimer.min(), vmax=no_dimer.max()))
# lattice_energy=Energy(lattice_dimer)
# plt.title(f'tau = {500} and Energy = {lattice_energy}')
# cbar = plt.colorbar(label='Dark color = oldest dimer and Light color = youngest dimer', orientation='vertical', fraction=0.046, pad=0.1)
# cbar.set_ticks([])

# plt.show()



'''
For capturing snapshots. If you want to see lattice on any specific step then you can uncomment the above commented code 
'''
lattice_dimer=annealing(50,10000) 
save_folder = 'HW_3_Fig_5_snapshots'
os.makedirs(save_folder, exist_ok=True)
cmap = plt.cm.viridis
cmap.set_bad(alpha=0)

fig, ax = plt.subplots()
im = ax.imshow(np.zeros((1, 1)), cmap=cmap, interpolation='none', norm=Normalize(vmin=0, vmax=1))
cbar = plt.colorbar(im, label='Dimer age(Dark = oldest)', orientation='vertical', fraction=0.046, pad=0.1)
cbar.set_ticks([])  

for i in range(0, 10001, 500):
    lattice_dimer= annealing(50, i)
    no_dimer = np.ma.masked_equal(lattice_dimer, 0)
    ax.clear()
    im = ax.imshow(no_dimer, cmap=cmap, interpolation='none', norm=Normalize(vmin=no_dimer.min(), vmax=no_dimer.max()))
    lattice_energy=Energy(lattice_dimer)
    plt.title(f'Dimer step = {i} and Energy = {lattice_energy}')
    plt.savefig(os.path.join(save_folder,f'lattice_dimer_step_{i}.png'))

plt.close()