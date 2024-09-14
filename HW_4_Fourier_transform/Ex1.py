# Question 1(a)

import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft

sunspots_months=np.loadtxt("sunspots.txt",float)
months=sunspots_months[:,0]
sunspots=sunspots_months[:,1]
plt.plot(months,sunspots, color='red', linewidth=0.5)
plt.xlabel("Months", fontname='Times New Roman',fontsize=18)
plt.ylabel("Number of sunspots", fontname='Times New Roman',fontsize=18)
plt.xticks(fontname='Times New Roman',fontsize=16)
plt.yticks(fontname='Times New Roman',fontsize=16)
plt.show()

# We hae 24 cycles in an overall span of 3000 months, that means a period of about 125 months for each cycle.



# Question 1(b)

ck=fft.fft(sunspots)
k=np.arange(0,len(ck),1)
plt.plot(k,abs(ck)**2,color='red', linewidth=1)
plt.title("Power Spectrum of the Sunspot Signal", fontname='Times New Roman',fontsize=18)
plt.xlabel("k", fontname='Times New Roman',fontsize=18)
plt.ylabel("$|C_{\\mathrm{k}}^{\\mathrm{2}}|$", fontname='Times New Roman', fontsize=18)
plt.show()



# Question 1(c)

t=[abs(ck[i])**2 for i in range(1,len(ck))]
tmax=max(t)
kmax=t.index(tmax) # maximum fourier coefficient occur at this k
print(" Maximum fourier coefficient occur at",kmax)
print("Period of the sine wave with the value of kmax is",len(ck)/kmax)
x=[1+kmax]
y=[tmax]
plt.plot(k[1:100],abs(ck[1:100])**2,color='red', linewidth=1)
plt.scatter(x,y,c="red",s=10)
plt.xlabel("k", fontname='Times New Roman',fontsize=18)
plt.ylabel("$|C_{\\mathrm{k}}^{\\mathrm{2}}|$", fontname='Times New Roman', fontsize=18)
plt.show()


'''
The maximum amplitude out of fourier coefficients occurs at k=23. k/N will be the frequency at which this sinusoidal wave 
oscillates and period of this will be N/k=136.65 months, which is quite close to what we obtained from eyeballing the graph.
'''





