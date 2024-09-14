# Question (a)
# T0 find the solution of ax^2 + bx + c = 0 
import math
def roots_1(a,b,c):
    roots = []
    x1 = (-b+math.sqrt((b**2)-4*a*c))/(2*a)
    x2 = (-b-math.sqrt((b**2)-4*a*c))/(2*a)
    roots.append(x1)
    roots.append(x2)
    return roots

print("Roots are")
print(roots_1(0.001,1000,0.001))


# Question (b)

# x = [{-b +- (b^2 - 4ac)^1/2}/2ac] *  [-b ∓ (b^2−4ac)^1/2]
# x = [(-b)^2 - (b^2 - 4ac)] / 2a[-b -+ (b^2 - 4ac)^1/2]
# x = (4ac) / 2a[-b -+ (b^2 - 4ac)^1/2]
# x = (2c) / [-b -+ (b^2 - 4ac)^1/2]

import math
def roots_rearr(a,b,c):
    roots = []
    x3 = 2*c/(-b-math.sqrt((b**2)-4*a*c))
    x4 = 2*c/(-b+math.sqrt((b**2)-4*a*c))
    roots.append(x3)
    roots.append(x4)
    return roots

print("Rearranged equation roots are")
print(roots_rearr(0.001,1000,0.001))

# The two functions differ in their results. This happens when b^2>>|4ac|. In either of the two
# formulae, finding one of the roots involve subtracting two nearly equal quantities (Discriminant is almost b^2).
# This results in loss of numerical precision corresponding to that root. So in one formula if Root 1 is more accurate
# than Root 2, then in the other formula, Root2 is more accurate than Root 1. Now the question is which root is more accurate
# in which formula?
# Well that depends on the sign of b in ax^2+bx+c=0
# When b>0, -b-(b^2-4ac)^1/2 doesn't involve subtracting two nearly equal quantities. So Root 2 of roots_1 function and
# Root 1 of roots_rearr function gives accurate results.
# Converse happens for b<0




# Question (c)

def roots_acc(a,b,c):
    if b>0:
        return [roots_rearr(a,b,c)[0],roots_1(a,b,c)[1]]
    else:
        return [roots_1(a,b,c)[10],roots_rearr(a,b,c)[1]]
    
print('Accurate Roots are:')
print(roots_acc(.001, 1000, .001))



# Solving the equation by bisection method also to confirm the answer
def func_quad(a,b,c,r):
     return a*(r**2)+b*r+c


def roots_num(a,b,c,lst_range):  # lst range is 2D list of all ranges in which root is present
    roots=[]
    for i in (lst_range):
        m = i[0]
        n = i[1]
        r=(m+n)/2
        y=func_quad(a,b,c,r)
        
        while abs(y) > (10**-6):
            if func_quad(a,b,c,r)*func_quad(a,b,c,m) < 0:
                
                n=r
            else:
                
                m=r
        
            r = (m+n)/2
            y=func_quad(a,b,c,r)
        roots.append(r)        
    
    return roots

print("Roots obtained from bisection method")
print(roots_num(0.001,1000,0.001,[[-0.5*(10**6),-1.5*(10**6)],[-0.5*(10**-6),-1.5*(10**-6)]]))
print("completed")