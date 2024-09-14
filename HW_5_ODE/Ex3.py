# Question 3(a)

import matplotlib.pyplot as plt

sigma = 10
r = 28
b = 8/3

def deriv_x(x,y,z):
    return sigma*(y-x)

def deriv_y(x,y,z):
    return r*x-y-x*z

def deriv_z(x,y,z):
    return x*y-b*z


def lorentz(time,N):
    x=[0]
    y=[1]
    z=[0]
    h=(time/N)
    Time = [0+i*h for i in range(N+1)]
    
    for i in range(len(Time)-1):
        k1_x = h*deriv_x(x[i], y[i], z[i])
        k1_y = h*deriv_y(x[i], y[i], z[i])
        k1_z = h*deriv_z(x[i], y[i], z[i])

        k2_x = h*deriv_x(x[i] + k1_x/2, y[i] + k1_y/2, z[i] + k1_z/2)
        k2_y = h*deriv_y(x[i] + k1_x/2, y[i] + k1_y/2, z[i] + k1_z/2)
        k2_z = h*deriv_z(x[i] + k1_x/2, y[i] + k1_y/2, z[i] + k1_z/2)

        k3_x = h*deriv_x(x[i] + k2_x/2, y[i] + k2_y/2, z[i] + k2_z/2)
        k3_y = h*deriv_y(x[i] + k2_x/2, y[i] + k2_y/2, z[i] + k2_z/2)
        k3_z = h*deriv_z(x[i] + k2_x/2, y[i] + k2_y/2, z[i] + k2_z/2)

        k4_x = h*deriv_x(x[i] + k3_x, y[i] + k3_y, z[i] + k3_z)
        k4_y = h*deriv_y(x[i] + k3_x, y[i] + k3_y, z[i] + k3_z)
        k4_z = h*deriv_z(x[i] + k3_x, y[i] + k3_y, z[i] + k3_z)

        x_1 = x[i] + (k1_x + 2*k2_x + 2*k3_x +k4_x)/6
        y_1 = y[i] + (k1_y + 2*k2_y + 2*k3_y +k4_y)/6
        z_1 = z[i] + (k1_z + 2*k2_z + 2*k3_z +k4_z)/6
        

        x.append(x_1)
        y.append(y_1)
        z.append(z_1)

    return x,y,z,Time

x,y,z,time = lorentz(50,10000)


plt.plot(time,y,c='red')
plt.xlabel("Time")
plt.ylabel("y")
plt.show()



# Question 3(b)

plt.plot(x,z,c='red')
plt.xlabel("x")
plt.ylabel("z")
plt.show()