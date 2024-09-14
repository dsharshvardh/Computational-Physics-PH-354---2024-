# Question (a)

import numpy as np
print("Question (a)")
# defining function
def f(x):
    return x**4 - 2*x + 1

# Writing trapezoidal integration

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
    frac_error = (I-4.4)/4.4
    print("Then integration by trapezoidal rule with", i ,"slices and its fractional error with the correct value 4.4 is")
    return [I,frac_error]

print(trap(0,2,10))  
print()


# Question (b)
print("Question (b)")

# Writing simpson's integration
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
    frac_s_error = (I_s-4.4)/4.4
    print("Then integration by simpson's rule with", i, "slices and its fractional error with the correct value 4.4 is")
    return [I_s,frac_s_error]

print(simps(0,2,10))
print()



# Question (c)

print("Question (c)")

print(trap(0,2,100))  
print()
print(simps(0,2,100))
print()
print(trap(0,2,1000)) 
print() 
print(simps(0,2,1000))
print() 


# So the value of Integration calculated by Simpson's method is always having more accurate values than Trapezoidal method which is expected and with increase in number of slices, the accuracy is increasing. 



# Question (d)
print("Question (d)")
trap_10 = trap(0,2,10) # 10 slices
print(trap_10[0])
print()
trap_20 = trap(0,2,20) # 20 slices
print(trap_20[0])
print()
print("Error = ",np.abs(-((trap_20[0]-trap_10[0])/(1-(20/10)**2))))
print("Analytical error = ",trap_20[0]-4.4)

# They do not agree perfectly because in method 1.d) error calculation has been done using Taylor series approximation and neglecting terms of order h^4  
