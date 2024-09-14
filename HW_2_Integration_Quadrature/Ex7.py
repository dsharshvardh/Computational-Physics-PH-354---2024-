import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps



# Question (a) 
# Here I am making a function which will read a list having data and make a list of dw/dx, dw/dy and Intensity of the given list 
# To find differentiation, I will have to use all F = froward difference method, B = Backward difference method and C = Central difference method (giving shorthand notation as F,B,C )
# Because on edges, we will not be able to use central difference method as it needs two before and after points of the edge point which is not possible. So I have to use F and B method
# for them according to requirements. And those data points of altitude which are not in edges can be done by C method. 

def Intensity(lst,phi,h):  # F = froward difference, B = Backward difference, C = Central difference
    r = len(lst)
    c = len(lst[0])
    intensity=[[0 for i in range (c)]for j in range(r)]
    dw_dx=[[0 for i in range (c)]for j in range(r)]
    dw_dy=[[0 for i in range (c)]for j in range(r)]
    for i in range(len(intensity)):
        for j in range(len(intensity[0])):
            if i==j==0:
                dwdx = (lst[i][j+1] - lst[i][j])/h  # F
                dwdy = (lst[i+1][j] - lst[i][j])/h  # F
            elif i==0 and j==(len(intensity[0])-1):
                dwdx = (lst[i][j] - lst[i][j-1])/h  # B
                dwdy = (lst[i+1][j] - lst[i][j])/h  # F
            elif i==0 and (j!=0 or j!=(len(intensity[0])-1)):
                dwdx = (lst[i][j+1] - lst[i][j-1])/(2*h)  # C
                dwdy = (lst[i+1][j] - lst[i][j])/h        # F  
            
            elif i==(len(intensity)-1) and j==0:
                dwdx = (lst[i][j+1] - lst[i][j])/h   # F
                dwdy = (lst[i][j] - lst[i-1][j])/h   # B
            elif i==(len(intensity)-1) and j==(len(intensity[0])-1):
                dwdx = (lst[i][j] - lst[i][j-1])/h   # B
                dwdy = (lst[i][j] - lst[i-1][j])/h   # B
            elif i==(len(intensity)-1) and (j!=(len(intensity[0])-1)):
                dwdx = (lst[i][j+1] - lst[i][j-1])/(2*h)  # C
                dwdy = (lst[i][j] - lst[i-1][j])/h        # B

            elif (i!=0 or i!=(len(intensity)-1)) and j==0:
                dwdx = (lst[i][j+1] - lst[i][j])/h        # F
                dwdy = (lst[i+1][j] - lst[i-1][j])/(2*h)  # C
            elif (i!=0 or i!=(len(intensity)-1)) and j==(len(intensity[0])-1):
                dwdx = (lst[i][j] - lst[i][j-1])/h  # B
                dwdy = (lst[i+1][j] - lst[i-1][j])/(2*h)  # F

            else:  
                dwdx = (lst[i][j+1] - lst[i][j-1])/(2*h)  # C
                dwdy = (lst[i+1][j] - lst[i-1][j])/(2*h)  # C
                
            dw_dx[i][j]=(dwdx) # it will make derivative list with respect to x
            dw_dy[i][j]=(dwdy) # it will make derivative list with respect to y
            intensity[i][j] = (math.cos(phi)*dwdx + math.sin(phi)*dwdy)/(math.sqrt(dwdx**2 + dwdy**2 +1)) # it will give me intensity 2D array
    return dw_dx,dw_dy,intensity




# Question (b)

alt_lst=np.loadtxt("altitude.txt") # loading data from file
dw_dx,dw_dy,intensity = Intensity(alt_lst,np.deg2rad(45),30000)

plt.imshow(intensity,cmap="twilight")
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.title("Density_plot")
plt.show()



# Question (c)

stm_lst=np.loadtxt("stm.txt") # loading data from file
dw_dx,dw_dy,intensity = Intensity(stm_lst,np.deg2rad(45),2.5)

plt.imshow(intensity,cmap="cubehelix")
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.title("Density_plot")
plt.show()