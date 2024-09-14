# Question(a)
# Two solve this integration by adaptive trapezoidal method, we will divide the range of x in two parts, 0 to 0.6 and 0.6 to 1 and do the itteration for both of them 
# seprately because from 0 to 0.6 curve of this function are having sharp curvature
print("Question(a)")
# define function
import math
def f(x):
    return (math.sin(math.sqrt(100*x)))**2 

print()
print("Integration by adaptive Trapezoidal method")
print()

def trap(a,b,i): # 'a' is lower limt, 'b' is upper limit and 'i' is the number of slices we want. 
    h = (b-a)/i
    x_lst = [a+i*h for i in range(i+1)]
    integ = 0
    for i in range(len(x_lst)):
        if i==0 or i==(len(x_lst)-1):
            integ = integ + f(x_lst[i])
        elif i!=0 and i!=(len(x_lst)-1):
            integ = integ + 2*f(x_lst[i])

    I = (h/2)*(integ)

    return I

I_old = 0
I_new = trap(0,0.6,1) # giving just random values
n=1
print("For the range of 0 to 0.6")
print(1,"slice    Integration = ",I_new)
while abs(I_new - I_old)>(10**-6):
    I_old = I_new
    n=n*2
    I_new = trap(0,0.6,n)
    error = abs(I_new - I_old)/I_new
    print(n,"slice    Integration = ",I_new,"    Error =",error)
    


print()
I_old_1 = 0
I_new_1 = trap(0.6,1,1) # giving just random values
n=1
print("For the range of 0.6 to 1")
print(1,"slice    Integration = ",I_new_1)
while abs(I_new_1 - I_old_1)>(10**-6):
    I_old_1 = I_new_1
    n=n*2
    I_new_1 = trap(0.6,1,n)
    error_1 = abs(I_new_1 - I_old_1)/I_new
    print(n,"slice    Integration = ",I_new_1,"    Error =",error_1)

print()
print("The integration of the function from 0 to 1 is",I_new+I_new_1)
print() 
print() 
print() 


# Question (b)
print("Question(b)")
# By Romberg method
print("Integration by romberg method")
print()
def romb(N): # N is till what step we want to go
    n=1
    I_rom = [[trap(0,1,n)]]
    while n != 2**(N-1):
        clmns=[]
        n=n*2
        clmns.append(trap(0,1,n))
        for i in range(len(I_rom[len(I_rom)-1])):
            elemnt = ((4**(i+1))*clmns[i] - I_rom[len(I_rom)-1][i])/((4**(i+1))-1)
            clmns.append(elemnt)
        I_rom.append(clmns)
    
    
    for i in range(len(I_rom)):   
        for j in range(len(I_rom[i])):
            print("{:10.10f}".format(I_rom[i][j]),end="    ")
        print()

    return I_rom[N-1][N-1]

I_old_rom = 0
I_new_rom = romb(1) 
N=1
print(1,"slice    Integration = ",I_new_rom)
print()
while abs(I_new_rom - I_old_rom)>(10**-6):
    I_old_rom = I_new_rom
    N=N+1
    I_new_rom = romb(N)
    error_rom = abs(I_new_rom - I_old_rom)/I_new_rom
    
    print((2**(N-1)),"slice    Integration = ",I_new_rom,"    Error =",error_rom)
    print()
print()
print()
print()



# Question (c) By adaptive simpson's rule
print("Question(c)")
print("Integration by adaptive Simpson's method")
print()
def simps(a,b,i): # 'a' is lower limt, 'b' is upper limit and 'i' is the number of slices we want. 
    h = (b-a)/i
    xs_lst = [a+i*h for i in range(i+1)]
    integ_s = 0
    for i in range(len(xs_lst)):
        if i==0 or i==(len(xs_lst)-1):
            integ_s = integ_s + f(xs_lst[i])
        
        elif i!=0 and i!=(len(xs_lst)-1) and i%2==1:
            integ_s = integ_s + 4*f(xs_lst[i])
        
        elif i!=0 and i!=(len(xs_lst)-1) and i%2==0:
            integ_s = integ_s + 2*f(xs_lst[i])

    I_s = (h/3)*(integ_s)
    return I_s

I_old_s = 0
I_new_s = simps(0,0.6,1) # giving just random values
n=1
print("For the range of 0 to 0.6")

while abs(I_new_s - I_old_s)>(10**-6):
    I_old_s = I_new_s
    n=n*2
    I_new_s = simps(0,0.6,n)
    error = abs(I_new_s - I_old)/I_new_s
    print(n,"slice    Integration = ",I_new_s,"    Error =",error)


print()
I_old_s_1 = 0
I_new_s_1 = trap(0.6,1,1) # giving just random values
n=1
print("For the range of 0.6 to 1")

while abs(I_new_s_1 - I_old_s_1)>(10**-6):
    I_old_s_1 = I_new_s_1
    n=n*2
    I_new_s_1 = simps(0.6,1,n)
    error_1 = abs(I_new_1 - I_old_1)/I_new
    print(n,"slice    Integration = ",I_new_s_1,"    Error =",error_1)

print()
print("The integration of the function from 0 to 1 is",I_new_s + I_new_s_1)


# The results confirm here that the order of reuired number of slices is as : -
#            Adaptive Trapezoidal < Adaptive Simpson < Romberg method