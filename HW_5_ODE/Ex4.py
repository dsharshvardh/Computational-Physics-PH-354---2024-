import numpy as np
import matplotlib.pyplot as plt


# Question 4(a)
# I have added an image named as HW5_Explaination_4(a)

# Question 4 (b) and (c)

def grav(i,j,r): # We take two bodies i and j here and calculate the gravitational acceleration on body i due to body j
    mass=[150,200,250]
    return mass[j]*(r[j]-r[i])/abs(r[j]-r[i])**3



def f(r,t):
    dr=np.zeros(6,dtype=complex)
    for i in range(len(dr)):
        obj=[]
        if i < 3:
            dr[i] = r[3+i]
        else:
            for j in range(3):
                if j != (i-3):
                    obj.append(j)            
            dr[i] = grav(i-3,obj[0],r)+grav(i-3,obj[1],r)

    return dr


# Considering 2d vectors as complex numbers in the argand plane (z=x+iy).

def three_body(h,rlist,t):
    time=[0]
    tolerance=1e-3

    while t<=2:
        r=rlist[-1]
        r1 = r.copy()
        r2=r.copy()
        for i in range(2):            
            k1 = h*f(r1, t)
            k2 = h*f(r1 + k1/2, t + h/2)
            k3 = h*f(r1 + k2/2, t + h/2)
            k4 = h*f(r1 + k3, t + h)
            r1 = r1+(k1+2*(k2+k3)+k4)/6

        k1 = 2*h*f(r2, t)
        k2 = 2*h*f(r2 + k1/2, t + h)
        k3 = 2*h*f(r2 + k2/2, t)
        k4 = 2*h*f(r2 + k3, t + 2*h)
        r2 = r2+(k1+2*(k2+k3)+k4)/6

        diff=[abs(r2[i]-r1[i]) for i in range(3)]
        l=max(diff)
        h_1=h*(30*h*tolerance/l)**(1/4)
        if h_1>2*h:
            h_1=2*h
        h=h_1
        # rnew=rk(hopt,r,T)
        k1 = h_1*f(r, t)
        k2 = h_1*f(r + k1/2, t + h)
        k3 = h_1*f(r + k2/2, t)
        k4 = h_1*f(r + k3, t + 2*h)
        r = r+(k1+2*(k2+k3)+k4)/6
        rnew = r

        rlist.append(rnew)
        t+=h_1
        time.append(t)
    return rlist


rlst=[np.array([3+1j,-1-2j,-1+1j,0,0,0])]
rlst1 = three_body(1e-6,rlst,0)
        

        

x1lst=[ele[0].real for ele in rlst1]
y1lst=[ele[0].imag for ele in rlst1]
x2lst=[ele[1].real for ele in rlst1]
y2lst=[ele[1].imag for ele in rlst1]
x3lst=[ele[2].real for ele in rlst1]
y3lst=[ele[2].imag for ele in rlst1]
plt.plot(x1lst,y1lst,label="First Body",c='red')
plt.plot(x2lst,y2lst,label="Second Body",c='blue')
plt.plot(x3lst,y3lst,label="Third Body",c='green')
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.show()





