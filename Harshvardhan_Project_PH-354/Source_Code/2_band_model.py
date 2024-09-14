import numpy as np
import math
import cmath
import matplotlib.pyplot as plt

def find_eigen(matrix):
    eigenvalues, eigenvectors = np.linalg.eig(matrix)
    unsorted_eigenvalues = [value.real for value in eigenvalues]
    real_eigenvalues = sorted(unsorted_eigenvalues)
    return real_eigenvalues #, eigenvectors


# for eigen values of a two band matrix :- ti = hoppings , ka,kb = kx*ax, ky*ay 
def twoband(t1,t2,t3,t4,t5,ka,kb):
    Ak = t2 + t5*(cmath.exp(complex(0,-ka)))
    Bk = 4*t4*(cmath.exp(complex(0,-(ka-kb)/2)))*math.cos(ka/2)*math.cos(kb/2)
    Ck = 2*(cmath.exp(complex(0,kb/2)))*math.cos(kb/2)*(t3+t1*(cmath.exp(complex(0,-ka))))
    Dk = 2*(cmath.exp(complex(0,kb/2)))*math.cos(kb/2)*(t1+t3*(cmath.exp(complex(0,-ka))))
    mat00 = Bk*(cmath.exp(complex(0,(ka-kb)/2)))
    mat01 = Ak+Ck*(cmath.exp(complex(0,(ka-kb)/2)))
    mat10 = np.conj(Ak) + np.conj(Ck)*(cmath.exp(complex(0,-(ka-kb)/2)))
    mat11 = Bk*(cmath.exp(complex(0,(ka-kb)/2)))
    matrix = [[mat00,mat01],[mat10,mat11]]
    eigen_values = find_eigen(matrix)
    return eigen_values

# K path list
kalst = np.linspace(0,math.pi,100)
kblst = np.linspace(0,math.pi,100)


Energy = []
sym_point = [0]
p=-1
for i in range(len(kalst)-1,-1,-1):
    Energy.append(twoband(-1.220,3.665,-0.205,-0.105,-0.055,kalst[i],kblst[len(kblst)-1]))
    p=p+1
sym_point.append(p)  

for i in range(len(kblst)-1,-1,-1):
    Energy.append(twoband(-1.220,3.665,-0.205,-0.105,-0.055,0,kblst[i]))
    p=p+1
sym_point.append(p)  

for i in kalst:
    Energy.append(twoband(-1.220,3.665,-0.205,-0.105,-0.055,i,0))
    p=p+1
sym_point.append(p)  

for i in kblst:
    Energy.append(twoband(-1.220,3.665,-0.205,-0.105,-0.055,kalst[len(kalst)-1],i))
    p=p+1
sym_point.append(p)  


E_plot=[] # Modifying Eigenvalue list (Energies to make it plotted)
for i in range(len(Energy[0])):
    E=[]
    for j in range(len(Energy)):
        E.append(Energy[j][i])
    E_plot.append(E)


symbols = ['S','Y',r'$\Gamma$','X','S']
x_tick_labels = [''] * 401 

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











