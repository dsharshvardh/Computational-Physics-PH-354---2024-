import numpy as np

# Lets generate a function of d dimensions

def monte_integ(d,N): # d = dimension and N = Number of points
    rand_points = []
    for i in range(N):
        varble_points=[]
        for j in range(d):
            varble_points.append(np.random.uniform(-1,1))
        rand_points.append(varble_points)

    fx=0
    for i in range(len(rand_points)):
        sum=0        
        for j in range(len(rand_points[0])):
            sum = sum + rand_points[i][j]**2
        if sum<=1:
            fx = fx+1
        else:
            continue
    
    I = ((2**d)/N)*fx
    return I



print("The area of sphere in 10 dimension with unit radius is" ,monte_integ(10,1000000))
print("The area of a circle(2D) with unit radius is", monte_integ(2,1000000))

