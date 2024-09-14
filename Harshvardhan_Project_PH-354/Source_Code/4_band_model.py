import numpy as np
import math
import cmath
import matplotlib.pyplot as plt


'''
In this file, we are plotting band structure for 4-band model
'''

def find_eigen(matrix):
    eigenvalues, eigenvectors = np.linalg.eig(matrix)
    unsorted_eigenvalues = [value.real for value in eigenvalues]
    real_eigenvalues = sorted(unsorted_eigenvalues)
    return real_eigenvalues , eigenvectors


# for eigen values of a four band matrix - ti = hoppings , ka,kb = kx*ax, ky*ay 

def fourband(t1,t2,t3,t4,t5,ka,kb):
    Ak = t2 + t5*(cmath.exp(complex(0,-ka)))
    Bk = 4*t4*(cmath.exp(complex(0,-(ka-kb)/2)))*math.cos(ka/2)*math.cos(kb/2)
    Ck = 2*(cmath.exp(complex(0,kb/2)))*math.cos(kb/2)*(t3+t1*(cmath.exp(complex(0,-ka))))
    Dk = 2*(cmath.exp(complex(0,kb/2)))*math.cos(kb/2)*(t1+t3*(cmath.exp(complex(0,-ka))))
    matrix = [[0,Ak,Bk,Ck],[np.conj(Ak),0,Dk,Bk],[np.conj(Bk),np.conj(Dk),0,Ak],[np.conj(Ck),np.conj(Bk),np.conj(Ak),0]]
    eigen_values,eigen_vec = find_eigen(matrix)
    eigen_values_new = sorted(eigen_values)
    return eigen_values_new


kalst = np.linspace(0,math.pi,100)
kblst = np.linspace(0,math.pi,100)

Energy = []
sym_point = [0]
p=-1
for i in range(len(kalst)-1,-1,-1):
    Energy.append(fourband(-1.220,3.665,-0.205,-0.105,-0.055,kalst[i],kblst[len(kblst)-1]))
    p=p+1
sym_point.append(p)  

for i in range(len(kblst)-1,-1,-1):
    Energy.append(fourband(-1.220,3.665,-0.205,-0.105,-0.055,0,kblst[i]))
    p=p+1
sym_point.append(p)  

for i in kalst:
    Energy.append(fourband(-1.220,3.665,-0.205,-0.105,-0.055,i,0))
    p=p+1
sym_point.append(p)  

for i in kblst:
    Energy.append(fourband(-1.220,3.665,-0.205,-0.105,-0.055,kalst[len(kalst)-1],i))
    p=p+1
sym_point.append(p)  


E_plot=[] # Modifying Eigenvalue list (Energies to make it plotted)
for i in range(len(Energy[0])):
    E=[]
    for j in range(len(Energy)):
        E.append(Energy[j][i])
    E_plot.append(E)



symbols = ['S','Y',r'$\Gamma$','X','S']
x_tick_labels = [''] * 401  # Create enough empty slots for 0 to 400
for x, symbol in zip(sym_point, symbols):
    x_tick_labels[x] = symbol

for sublist in E_plot:
    plt.plot(sublist,linewidth = 2)

plt.xticks(sym_point, [x for x in x_tick_labels if x != ''])
for x in sym_point:
    plt.axvline(x=x, color='black', linestyle='--', linewidth=1)  # Add vertical lines

plt.rcParams['font.family'] = 'Times New Roman'

plt.xticks(sym_point, [x for x in x_tick_labels if x != ''],fontsize=14)
plt.yticks(fontsize=14)
plt.xlim(0,399)
plt.xlabel('wave vector',fontsize=16)
plt.ylabel('Energy(eV)',fontsize=16)
plt.show()










