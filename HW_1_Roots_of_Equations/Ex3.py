import math
import scipy.constants as const

# Let lambda = b/T (b is displacement constant)

# Question (a)

h=const.h
c=const.c
k=const.k
cnst = h*c/k

def f(b):
    return ((cnst*math.exp(cnst/b))/(b*(math.exp(cnst/b)-1)**2)) - (5/(math.exp(cnst/b)-1))   # I got this equation after putting d[I(lambda)]/d(lambda) = 0 and putting T = b/lambda, we will find the value of b now when put this equation equals to zero..


a=0.001
b=0.005
c=(a+b)/2
# Finding roots by Bisection method
while abs(f(c))>(10**-6):
    if f(a)*f(c)<0:
        b=c
    else:
        a=c
    c=(a+b)/2

def lmbda(T):  # T is Temprature
    print("The maximum wavelength in nm and the displacement constant are")
    return (c*10**9)/T, c

print(lmbda(5000))


# Question (b)

def Temp(lmbda):
    print("The surface temprature of the sun is")
    return (c*10**9)/lmbda

print(Temp(502))