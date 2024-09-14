
import math
import scipy.constants as const
import numpy as np


w = 1 #nm
V = 20 #eV
R2 = (const.m_e*(w*1e-9)**2/(2*const.hbar**2))*(V*const.e)



def sym(a):
    if a**2>R2:
        print(" E>V for after this")
    else:
        return  -a*np.sin(a) + np.sqrt(R2-a**2)* np.cos(a)
      

def antisym(a):
    return a*np.cos(a) + np.sqrt(R2-a**2)* np.sin(a)


def roots(func,a,b):
    c = (a+b)/2.
    while abs(func(c))>=10**-8:
        c = (a*func(b)-b*func(a))/(func(b)-func(a))
        if (func(a)*func(c)>=0):
            a=c
        else:
            b=c
    return c
        
        

# Now lets find energies
print("The 6 first energy levels are")
for n in range(3):
    a = (1+3*n)-0.3*3
    b = (1+3*n)+0.3*3
    Energy = ((1/(2*const.m_e))*(2*const.hbar*roots(sym,a,b)/(w*1e-9))**2)/const.e
    if Energy<V:
        print("{:10.5f}".format(Energy),"eV - symmetric")
    
    a = (2.7+3*n)-0.3*3
    b = (2.7+3*n)+0.3*3
    Energy = ((1/(2*const.m_e))*(2*const.hbar*roots(antisym,a,b)/(w*1e-9))**2)/const.e
    if Energy<V:
        print("{:10.5f}".format(Energy),"eV - antisymmetric")



    


      
