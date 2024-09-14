import math
#Constants in S.I. unit
G=6.674*(10**(-11))
M=5.974*(10**24)
m=7.348*(10**22)
R=3.844*(10**8)
w=2.662*(10**-6)

# Two initial points
r0=0.01*(10**8)
r1=3.75*(10**8)

def f(r):
	return ((G*M)/r**2)-((G*m)/(R-r)**2)-(w**2)*r
	

#Applying Seacent method for root finding

f0=f(r0)
f1=f(r1)


while abs(f1)>(10**-6):
	r2 = (r0*f(r1)-r1*f(r0))/(f(r1)-f(r0))
	f0=f1
	f1=f(r2)
	r0=r1
	r1=r2
	
	

print("The lagrange point between earth and moon is=",r1)





