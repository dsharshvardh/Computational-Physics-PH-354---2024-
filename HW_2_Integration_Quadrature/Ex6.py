import math
import cmath
import numpy as np

def f(z):
    return cmath.exp(2*z)

def deriv(m,N):
    sum = 0
    for k in range(N):
        z = cmath.exp(complex(0,2*math.pi*k/N))
        sum = sum + f(z)*cmath.exp(complex(0,-2*math.pi*k*m/N))
    
    mth_deriv = math.factorial(m)*sum/N
    return mth_deriv

N = 10000

for m in range(21):
    derivt = deriv(m,N)
    round(np.real(derivt))
    print("The",m,"th derivative is ",round(np.real(derivt))," + (",(np.imag(derivt)),"j)")
    


