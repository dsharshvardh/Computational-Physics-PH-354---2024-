# Question (a)

# This process involves guess of a root and then keep improving the guess.
# In it we conver the from of function h(x)=0 to x=f(x) and then do itteration for approximate value of x
# For the given equation, one root will be 0 and other we will find by this code 
import numpy as np
import matplotlib.pyplot as plt
def f(c,x):
    return (1-np.exp(-c*x))
c=2
#initial guesses of root

x=0.5
iter = 0
while(abs(x-f(c,x))>=1e-6):
        x=f(c,x)
        iter +=1


print("The root found by fixed point iteration is",x,"with iterations",iter)
    


# Question (b) for values of c from 0 to 3 in steps of 0.01 
n=int(3/0.01)
c_lst = [i*0.01 for i in range(0,n+1)]

import math

root_lst = []
for c in c_lst:   
    x0 = 0.5
    while(abs(x0-f(c,x0))>=1e-6):
         x0=f(c,x0)
    
    root_lst.append(x0)





# Question (c) fixed point iteration with acceleration

# Here instead taking next value of x as x_n+1=f(x_n), we will take x_n=1 = x_n - (f(x_n)-x_n)/g(x_n) where g(x_n) = (f(f(x_n))-f(x_n))/(f(x_n)-x_n)-1

def g(f,c,x):
     return (f(c,f(c,x))-f(c,x))/(f(c,x)-x) - 1

c=2
x=0.5
iter=0
while(abs(x-f(c,x))>=1e-6):
    x = x-(f(c,x)-x)/g(f,c,x)
    iter+=1

print("The root found by fixed point iteration with acceleration is",x,"with iterations",iter)

print("Completed")

# Yes, my code is giving my almost same answer with very less itterations


# Graph for Question (b)
print("The graph of percolation threshold is")  

# Plot
plt.plot(c_lst, root_lst, marker='o', linestyle='-', color='r',markerfacecolor='b')
plt.xlabel('c')
plt.ylabel('x')
plt.title('Percolation threshold')

plt.show()
     



