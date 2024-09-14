import numpy as np
import math
import cmath
import matplotlib.pyplot as plt


'''
For eigen values of a matrix (For clarity of my code, I will use this function just for eigen values for ploting of bands 
and construct another function for both eigen values and eigen vectors)
'''

def find_eigen(matrix):
    eigenvalues, eigenvectors = np.linalg.eig(matrix)
    real_eigenvalues = [value.real for value in eigenvalues]
    new_eigenvalues = sorted(real_eigenvalues)
    return new_eigenvalues #, eigenvectors


def find_eigen_vec(matrix):
    eigenvalues, eigenvec = np.linalg.eig(matrix)
    real_eigenvalues = [value.real for value in eigenvalues]
    # Convert eigenvectors to a list
    eigenvectors = [vector.flatten().tolist() for vector in eigenvec.T]

    # Combine elements from both lists into pairs
    combined_list = list(zip(real_eigenvalues, eigenvectors))

    # Sort based on the first list
    combined_list.sort(key=lambda x: x[0])

    # Unzip the sorted pairs back into the original lists
    new_eigenvalues, new_eigenvectors = zip(*combined_list)

    # Convert to lists
    new_eigenvalues = list(new_eigenvalues)
    new_eigenvectors = list(new_eigenvectors)
    return new_eigenvalues , new_eigenvectors


# To find a modulus square of an imaginary number
def modulus_square(z):
    # z = a + bi
    a = z.real
    b = z.imag
    
    # |z|^2 = (a^2 + b^2)
    modulus_square_value = a**2 + b**2
    
    return modulus_square_value


'''
Lets construct a momentum space hamiltonian for z-PNR having width N. To construct these Hamiltonians, we will make Hamiltonians for different hoppings seprately and in the end, we will add them all
and then will run the calculations that we want.
'''

def Hmltn(t1,t2,t3,t4,t5,kb,N):

    # First lets make matrix having only t1 and t2 hoppings
    H0 = [] # Hamiltonian for hopping of atoms within the unit cell
    row_0=[]
    for i in range(2*N):
        if i!=1:
            row_0.append(0)
        else:
            row_0.append(-t1)
    
    H0.append(row_0)
    
    for j in range(1,2*N-1):
        rows = []
        for k in range(2*N):
            if j%2==1:
                if k == j-1:
                    rows.append(-t1)
                elif k == j+1:
                    rows.append(-t2)
                else:
                    rows.append(0)
            
            else: 
                if k == j-1:
                    rows.append(-t2)
                elif k == j+1:
                    rows.append(-t1)
                else:
                    rows.append(0)
        
        H0.append(rows)
    
    row_2N=[]
    for i in range(2*N):
        if i!=2*N-2:
            row_2N.append(0)
        else:
            row_2N.append(-t1)
    
    H0.append(row_2N)
        
    H1=[[0 for i in range(2*N)]for j in range(2*N)]

    row_lst = [0,]
    column_lst = []
    for i in range(0,2*N-1):
        a=3*(i+1)+i
        if a>(2*N-1):
            break
        else:
           row_lst.append(a)
        
        if (a+1)>(2*N-1):
            break
        else:
           row_lst.append(a+1)
    
    for i in range(len(row_lst)):
        a=row_lst[i]+1
        b=row_lst[i]-1
        if i%2==0 and a<=(2*N-1):
            column_lst.append(a)
        elif i%2==1 and b<=(2*N-1):
            column_lst.append(b)
        
    
    for i in range(len(row_lst)):
        H1[row_lst[i]][column_lst[i]] = -t1*cmath.exp(complex(0,-kb))
    
    
    H_1 = [[0 for i in range(2*N)]for j in range(2*N)]
    
    for i in range(len(H1)):
        for j in range(len(H1)):
            H_1[i][j] = np.conj(H1[j][i])
            

    H_t1_t2 = [[0 for i in range(len(H_1))]for j in range(len(H_1))] # matrix haivng all t1 and t2 hoppings
    for i in range(len(H_1)):
        for j in range(len(H_1)):
            H_t1_t2[i][j]=(H0[i][j]+H1[i][j]+H_1[i][j])

    H_net=H_t1_t2

    # t3 hopping

    H_t3=[[0 for i in range(2*N)]for j in range(2*N)]  # matrix haivng all t3 hoppings 
    H_t3_f=[[0 for i in range(2*N)]for j in range(2*N)]
    H_t3_b=[[0 for i in range(2*N)]for j in range(2*N)]
    columns_t3_lst=[]
    for i in range(2*N):
        if i%2==1:
            columns_t3_lst.append(i)
        else:
            continue
    
    rows_t3_lst=[]
    for i in columns_t3_lst:
        if i+3<=(2*N-1):
            rows_t3_lst.append(i+3)
        else:
            break
    

    for i in range(len(rows_t3_lst)):
        if i%2==0:
            H_t3_f[rows_t3_lst[i]][columns_t3_lst[i]] = (-t3*(1+cmath.exp(complex(0,-kb))))

        else:
            H_t3_f[rows_t3_lst[i]][columns_t3_lst[i]] = (-t3*(1+cmath.exp(complex(0,kb))))

    for i in range(len(H_t3_b)):
        for j in range(len(H_t3_b)):
            H_t3_b[i][j]=np.conj(H_t3_f[j][i])

    for i in range(len(H_t3_b)):
        for j in range(len(H_t3_b)):
            H_t3[i][j] = H_t3_f[i][j]+H_t3_b[i][j]

    

    # t4 hopping 
        
    H_t4 =[[0 for i in range(2*N)]for j in range(2*N)]

    # lets do forward hoppings first and then backward hoppings (hopping from i to i+2)
    H_t4_f = [[0 for i in range(2*N)]for j in range(2*N)]
    H_t4_b = [[0 for i in range(2*N)]for j in range(2*N)]
    rows_t4_lst_d = [] # rows which will do hopping -ky direction unit cell 
    rows_t4_lst_u = [] # rows which will do hopping +ky direction unit cell 
    columns__t4_lst_d = [] # Columns which will do hopping -ky direction unit cell 
    columns__t4_lst_u = [] # Columns which will do hopping +ky direction unit cell
    
    i=0
    for j in range(3,2*N):
        if i==0 or i==1:
           rows_t4_lst_u.append(j)
        i = i+1
        if i==4:
           i=i-4
        else:
           continue
    
    for i in range(2,2*N):
        k=0
        for j in rows_t4_lst_u:
            if i==j:
                k=k+1
                break
            elif i!=j:
                continue
        if k==0:
            rows_t4_lst_d.append(i)
        else:
            continue
    
    # Now columns
    for i in rows_t4_lst_u:
        columns__t4_lst_u.append(i-2)
    
    for j in rows_t4_lst_d:
        columns__t4_lst_d.append(j-2)

    for i in range(len(rows_t4_lst_u)):
         H_t4_f[rows_t4_lst_u[i]][columns__t4_lst_u[i]]=-t4*(1+cmath.exp(complex(0,-kb)))

    for i in range(len(rows_t4_lst_d)):
         H_t4_f[rows_t4_lst_d[i]][columns__t4_lst_d[i]]=-t4*(1+cmath.exp(complex(0,kb)))
    
    for i in range(2*N):
        for j in range(2*N):
            H_t4_b[i][j] = np.conj(H_t4_f[j][i])

    
    for i in range(2*N):
        for j in range(2*N):
            H_t4[i][j] = H_t4_f[i][j] + H_t4_b[i][j] # overall H by t4 hopping
    

    # t5 hopping
    H_t5 =[[0 for i in range(2*N)]for j in range(2*N)]
    H_t5_f =[[0 for i in range(2*N)]for j in range(2*N)]
    H_t5_b =[[0 for i in range(2*N)]for j in range(2*N)]
    column_t5_lst=[]
    row_t5_lst=[]
    for i in range(2*N-3):
        if i%2==0:
            column_t5_lst.append(i)
        else:
            continue
    for j in column_t5_lst:
        row_t5_lst.append(j+3)

    for i in range(len(row_t5_lst)):
        H_t5_f[row_t5_lst[i]][column_t5_lst[i]]=-t5

    for i in range(2*N):
        for j in range(2*N):
            H_t5_b[i][j] = np.conj(H_t5_f[j][i])


    H_t5_f_2 =[[0 for i in range(2*N)]for j in range(2*N)]
    H_t5_b_2 =[[0 for i in range(2*N)]for j in range(2*N)]
    column_t5_lst_2=[]
    row_t5_lst_2=[]
    for i in range(2*N-1):
        if i%2==1:
            column_t5_lst_2.append(i)
        else:
            continue
    
    for j in column_t5_lst_2:
        row_t5_lst_2.append(j+1)

    for i in range(len(row_t5_lst_2)):    
            H_t5_f_2[row_t5_lst_2[i]][column_t5_lst_2[i]] = -2*t5*(math.cos(kb))



    for i in range(2*N):
        for j in range(2*N):
            H_t5_b_2[i][j] = np.conj(H_t5_f_2[j][i])

    
    for i in range(2*N):
        for j in range(2*N):
            H_t5[i][j] = H_t5_f[i][j] + H_t5_b[i][j] + H_t5_f_2[i][j] + H_t5_b_2[i][j] # overall H by t5 hopping
        
    
    '''
    Overall tight binding Hamiltonian with all hoppings 
    '''
    H_net=[[0 for i in range(2*N)]for j in range(2*N)]
    for i in range(len(H_t1_t2)):
        for j in range(len(H_t1_t2)):
            H_net[i][j] = H_t1_t2[i][j] + H_t3[i][j] + H_t4[i][j]  + H_t5[i][j]
        
    return H_net 


'''
Lets plot bandstructures for 50-zPNR with three different ratios |t2/t1| = 1,2,3
'''

kblst=np.linspace(-math.pi,math.pi,200) # k-path

Energies_t2_t1=[] # Collecting eigen values for k-path grid
Energies_t2_2t1 = []
Energies_t2_3t1 = []

for i in kblst:
    eigen_lst,eigen_vec = find_eigen_vec(Hmltn(1.22,-1.22,0.205,0.105,0.055,i,50))
    Energies_t2_t1.append(eigen_lst)

    eigen_lst_t2_2t1,eigen_vec_t2_2t1 = find_eigen_vec(Hmltn(1.22,-2.44,0.205,0.105,0.055,i,50))
    Energies_t2_2t1.append(eigen_lst_t2_2t1)

    eigen_lst_t2_3t1,eigen_vec_t2_3t1 = find_eigen_vec(Hmltn(1.22,-3.66,0.205,0.105,0.055,i,50))
    Energies_t2_3t1.append(eigen_lst_t2_3t1)

E_plot_t2_t1=[] # Modifying Eigenvalue list (Energies to make it plotted)
E_plot_t2_2t1=[]
E_plot_t2_3t1=[]
for i in range(len(Energies_t2_3t1[0])):
    E = []
    E2 = []
    E3 = []
    for j in range(len(Energies_t2_t1)):
        E.append(Energies_t2_t1[j][i])
    E_plot_t2_t1.append(E)

    for j in range(len(Energies_t2_2t1)):
        E2.append(Energies_t2_2t1[j][i])
    E_plot_t2_2t1.append(E2)
    
    for j in range(len(Energies_t2_3t1)):
        E3.append(Energies_t2_3t1[j][i])
    E_plot_t2_3t1.append(E3)
 


'''
Plotting Band structure for 50 width (100 bands)
'''

# t2 = t1

N = len(E_plot_t2_t1)
plt.figure(figsize=(2, 6)) 
for i, y_values in enumerate(E_plot_t2_t1):
    if i == N // 2 or i == (N // 2) - 1:  # Color N/2 and (N/2)+1 bands in red
        plt.plot(kblst, y_values, color='red', label=f'List {i+1}',linewidth=1)
    else:
        plt.plot(kblst, y_values, color='blue', label=f'List {i+1}',linewidth=1)


plt.xlim(-math.pi, math.pi)
plt.ylim(-3, 4)
plt.xlabel('k',fontsize=16)
plt.ylabel('Energy (eV)',fontsize=16)
plt.title('$t_2$ = |$t_1$|',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.show()


# t2 = 2|t1|

N = len(E_plot_t2_2t1)
plt.figure(figsize=(2, 6)) 
for i, y_values in enumerate(E_plot_t2_2t1):
    if i == N // 2 or i == (N // 2) - 1:  # Color N/2 and (N/2)+1 bands in red
        plt.plot(kblst, y_values, color='red', label=f'List {i+1}',linewidth=1)
    else:
        plt.plot(kblst, y_values, color='blue', label=f'List {i+1}',linewidth=1)

plt.xlim(-math.pi, math.pi)
plt.ylim(-4, 4)
plt.xlabel('k',fontsize=16)
plt.ylabel('Energy (eV)',fontsize=16)
plt.title('$t_2$ = 2|$t_1$|',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.show()


# t2 = 3|t1|

N = len(E_plot_t2_3t1)
plt.figure(figsize=(2, 6)) 
for i, y_values in enumerate(E_plot_t2_3t1):
    if i == N // 2 or i == (N // 2) - 1:  # Color N/2 and (N/2)+1 bands in red
        plt.plot(kblst, y_values, color='red', label=f'List {i+1}',linewidth=1)
    else:
        plt.plot(kblst, y_values, color='blue', label=f'List {i+1}',linewidth=1)

plt.xlim(-math.pi, math.pi)
plt.ylim(-5, 5)
plt.xlabel('k',fontsize=16)
plt.ylabel('Energy (eV)',fontsize=16)
plt.title('$t_2$ = 3|$t_1$|',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.show()



'''
Now we will find the probability amplitude of different sites in the quasi-flat band for k = 0 of a zigzag phosphorene nanoribbon for different ratios of |t2/t1|. 
'''

import matplotlib.pyplot as plt

# t1 = |t2|  50-zPNR
eigen_val_lst_1,eigen_vec_lst_1 = find_eigen_vec(Hmltn(1.22,-1.22,0,0,0,0,50)) 

prob_amp_1=[]
for i in eigen_vec_lst_1[49]:
    a=modulus_square(i)
    prob_amp_1.append(a)

xi = [i for i in range(100)]


# t1 = 2|t2|  50-zPNR
eigen_val_lst_2,eigen_vec_lst_2 = find_eigen_vec(Hmltn(1.22,-2.44,0,0,0,0,50))

prob_amp_2=[]
for i in eigen_vec_lst_2[49]:
    a=modulus_square(i)
    prob_amp_2.append(a)


# t1 = 3|t2|  50-zPNR
eigen_val_lst_3,eigen_vec_lst_3 = find_eigen_vec(Hmltn(1.22,-3.66,0,0,0,0,50))

prob_amp_3=[]
for i in eigen_vec_lst_3[49]:
    a=modulus_square(i)
    prob_amp_3.append(a)



import matplotlib.pyplot as plt

'''
For plotting Probability amplitude for different hoppigs.
It demands two differents scaling in y axis from 0 to 0.05 and 0.05 to 0.3 to compare the result with reported graph in Fig 4 of paper scaling laws of band gap. 
So, we manually transform the data for nonlinear scaling.
'''
def nonlinear_scale(value):
    if value <= 0.05:
        return value / 0.05 * 0.5  # Map [0, 0.05] to [0, 0.5]
    else:
        return 0.5 + (value - 0.05) / 0.25 * 0.5  # Map (0.05, 0.3] to (0.5, 1]

transformed_prob_amp_1 = [nonlinear_scale(value) for value in prob_amp_1]
transformed_prob_amp_2 = [nonlinear_scale(value) for value in prob_amp_2]
transformed_prob_amp_3 = [nonlinear_scale(value) for value in prob_amp_3]

# Set the font family to Times New Roman
plt.rcParams['font.family'] = 'Times New Roman'

# Plotting the graphs with transformed data
plt.plot(xi, transformed_prob_amp_1, marker='o', markersize=2, linestyle='-', linewidth=1, label='$t_2$ = |$t_1$|')
plt.plot(xi, transformed_prob_amp_2, marker='s', markersize=2, linestyle='--', linewidth=1, label='$t_2$ = 2|$t_1$|')
plt.plot(xi, transformed_prob_amp_3, marker='D', markersize=2, linestyle='-', linewidth=1, label='$t_2$ = 3|$t_1$|')

plt.xlabel('X$_i$', fontsize=18)
plt.ylabel(r'$|\Psi_i|^2$', fontsize=18)
plt.title('Probability amplitude of flat band of 50-zPNR at k=0 ', fontsize=16)
plt.legend(fontsize=18)
plt.xticks(fontsize=16)
plt.yticks([0, 0.5, 1], ['0', '0.05', '0.3'], fontsize=16)
plt.show()



'''
For probability amplitude of the upper valence band eigenstate for k = 0 of a zigzag phosphorene nanoribbon for |t2/t1|=3 using different width. 
'''

# 14-z-PNR 
eigen_val_lst_14,eigen_vec_lst_14 = find_eigen_vec(Hmltn(1.22,-3.66,0.205,0.105,0,0,14))

prob_amp_14=[]
for i in eigen_vec_lst_14[13]:
    a=modulus_square(i)
    prob_amp_14.append(a)


xi_14 = [i for i in range(-14,14,1)]


# 6-z-PNR
eigen_val_lst_6,eigen_vec_lst_6 = find_eigen_vec(Hmltn(1.22,-3.66,0.205,0.105,0,0,6))

prob_amp_6=[]
for i in eigen_vec_lst_6[5]:
    a=modulus_square(i)
    prob_amp_6.append(a)

xi_6 = [i for i in range(-6,6,1)]


import matplotlib.pyplot as plt

# For plotting Probability amplitude for different width
plt.rcParams['font.family'] = 'Times New Roman'

plt.plot(xi_14, prob_amp_14, marker='.', linestyle='--',linewidth=1, label='14-zPNR')
plt.plot(xi_6, prob_amp_6, marker='.', linestyle='-.', linewidth=1, label='6-zPNR')
plt.xlabel('X$_i$',fontsize=18)
plt.ylabel(r'$|\Psi_i|^2$',fontsize=18)
plt.title('Probability amplitude of flat band at k=0',fontsize=16)

plt.legend(fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.show()



'''
For effect of in plane external electric field in z-PNR
'''

def Eext_pert_mat(theta,a,b,Eext,N): # a = t1 hopping bond distance, b = t2 hopping bond distance, theta is angle between two t1 bonds, N = width
    mat=[[0 for i in range(2*N)]for j in range(2*N)] 
    dist=0
    e = 1 # Finding energy in eV
    for i in range(1,2*N): # Taking 0th atom as reference
        if i%2!=0:
            dist=dist + (a*math.cos(math.radians(theta)))
        else:
            dist=dist + 0.797284  # distance parallel with electric field between 2 atoms connected by t2
        
        mat[i][i] = -e*Eext*dist
    
    return mat  # mat is the perturbation matrix that we have to add with over tight binding Hamiltonian


# t1 = 3|t2|,  10-zPNR
Energies_10z_unpert=[]
Energies_10z_pert=[] # Collecting eigen values for k-path grid
for i in kblst:
    Hmltn_10zPNR = Hmltn(1.22,-3.66,0.205,0.105,0.055,i,10)    # 0.205,0.105,0.055   
    Eext_10zPNR_unpert = Eext_pert_mat(98.15/2,2.164,2.207,0,10) 
    Eext_10zPNR_pert = Eext_pert_mat(98.15/2,2.164,2.207,0.016,10)
    n=len(Hmltn_10zPNR)
    overall_mat_unpert=[[Hmltn_10zPNR[i][j]+Eext_10zPNR_unpert[i][j] for j in range(n)]for i in range(n)] 
    overall_mat_pert=[[Hmltn_10zPNR[i][j]+Eext_10zPNR_pert[i][j] for j in range(n)]for i in range(n)]   
    eigen_lst_10z_unpert,eigen_vec10z_unpert = find_eigen_vec(overall_mat_unpert) 
    eigen_lst_10z_pert,eigen_vec10z_pert = find_eigen_vec(overall_mat_pert)
    Energies_10z_unpert.append(eigen_lst_10z_unpert)
    Energies_10z_pert.append(eigen_lst_10z_pert)


E_plot_10z_unpert=[] # Modifying Eigenvalue list (Energies to make it plotted)
for i in range(len(Energies_10z_unpert[0])):
    E=[]
    for j in range(len(Energies_10z_unpert)):
        E.append(Energies_10z_unpert[j][i])
    E_plot_10z_unpert.append(E)

E_plot_10z_pert=[] # Modifying Eigenvalue list (Energies to make it plotted)
for i in range(len(Energies_10z_pert[0])):
    E=[]
    for j in range(len(Energies_10z_pert)):
        E.append(Energies_10z_pert[j][i])
    E_plot_10z_pert.append(E)


import matplotlib.pyplot as plt

'''
Unperturbed
'''
N = len(E_plot_10z_unpert)
plt.figure(figsize=(2, 6)) 
for i, y_values in enumerate(E_plot_10z_unpert):
    if i == N // 2 or i == (N // 2) - 1:  # Color N/2 and (N/2)+1 bands in red
        plt.plot(kblst, y_values, color='red', label=f'List {i+1}',linewidth=1)
    else:
        plt.plot(kblst, y_values, color='blue', label=f'List {i+1}',linewidth=1)

plt.xlim(0, math.pi)
plt.ylim(-6, 6) 
plt.xlabel('k',fontsize=18)
plt.ylabel('Energy (eV)',fontsize=18)
plt.title('$E_{ext}$ = 0',fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.show()


'''
Perturbed E = 0.016 
'''

N = len(E_plot_10z_pert)
plt.figure(figsize=(2, 6)) 
for i, y_values in enumerate(E_plot_10z_pert):
    if i == N // 2 or i == (N // 2) - 1:  # Color N/2 and (N/2)+1 bands in red
        plt.plot(kblst, y_values, color='red', label=f'List {i+1}',linewidth=1)
    else:
        plt.plot(kblst, y_values, color='blue', label=f'List {i+1}',linewidth=1)

plt.xlim(0, math.pi)
plt.ylim(-6, 6) 
plt.xlabel('k',fontsize=18)
plt.ylabel('Energy (eV)',fontsize=18)
plt.title('$E_{ext}$ = 0.016',fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.show()


'''
Probability amplitude of perturbed and unperturbed Hamiltonian
'''

Hmltn_10zPNR = Hmltn(1.22,-3.66,0.205,0.105,0.055,0,10)    # 0.205,0.105,0.055 
Eext_10zPNR = Eext_pert_mat(30,2.164,2.207,0.016,10)
n=len(Hmltn_10zPNR)
overall_mat=[[Hmltn_10zPNR[i][j]+Eext_10zPNR[i][j] for j in range(n)]for i in range(n)]    
eigen_lst_10z,eigen_vec10z = find_eigen_vec(overall_mat)
eigen_lst_10z_unpert,eigen_vec10z_unpert = find_eigen_vec(Hmltn_10zPNR)

'''
Perturbed
'''

prob_amp_10zPNR_cond=[]
for i in eigen_vec10z[10]:
    a=modulus_square(i)
    prob_amp_10zPNR_cond.append(a)

prob_amp_10zPNR_val=[]
for i in eigen_vec10z[9]:
    a=modulus_square(i)
    prob_amp_10zPNR_val.append(a)

'''
Unperturbed
'''

prob_amp_10zPNR_cond_unpert=[]
for i in eigen_vec10z_unpert[10]:
    a=modulus_square(i)
    prob_amp_10zPNR_cond_unpert.append(a)

prob_amp_10zPNR_val_unpert=[]
for i in eigen_vec10z_unpert[9]:
    a=modulus_square(i)
    prob_amp_10zPNR_val_unpert.append(a)

'''
We want to make a bubble plot in which we will make our unit cell for flat bands at k=0 and size of red and blue color bubble will tell us
probability amplitude of the sites. Here, lst and lst2 make unit cell for upper flat band at k=0 and lst1 and lst3 make unit cell for lower flat band at k=0.
'''

lst=[[0,0.5]]   
k=0
for i in range(1,len(prob_amp_10zPNR_cond)):
    if k==0:
        if i%2!=0:
            a=lst[i-1][0]+0.28
            b=lst[i-1][1]+0.28
            k=k+1
        else:
            a=lst[i-1][0]+0.1
            b=lst[i-1][1]
        
    elif k==1:
        if i%2!=0:
            a=lst[i-1][0]+0.28
            b=lst[i-1][1]-0.28
            k=0
        else:
            a=lst[i-1][0]+0.1
            b=lst[i-1][1]

    lst.append([a,b]) 

lst2=[[0,0.5]]   
k=0
for i in range(1,len(prob_amp_10zPNR_cond)):
    if k==0:
        if i%2!=0:
            a=lst2[i-1][0]+0.28
            b=lst2[i-1][1]-0.28
            k=k+1
        else:
            a=lst2[i-1][0]+0.1
            b=lst2[i-1][1]
        
    elif k==1:
        if i%2!=0:
            a=lst2[i-1][0]+0.28
            b=lst2[i-1][1]+0.28
            k=0
        else:
            a=lst2[i-1][0]+0.1
            b=lst2[i-1][1]

    lst2.append([a,b]) 



# lst1 is unit cell for below flat band at k=0
    
lst1=[[0,-0.5]] 
k=0
for i in range(1,len(prob_amp_10zPNR_val)):
    if k==0:
        if i%2!=0:
            a=lst1[i-1][0]+0.28
            b=lst1[i-1][1]+0.28
            k=k+1
        else:
            a=lst1[i-1][0]+0.1
            b=lst1[i-1][1]
        
    elif k==1:
        if i%2!=0:
            a=lst1[i-1][0]+0.28
            b=lst1[i-1][1]-0.28
            k=0
        else:
            a=lst1[i-1][0]+0.1
            b=lst1[i-1][1]

    lst1.append([a,b]) 


# lst3 is unit cell for below flat band at k=0
    
lst3=[[0,-0.5]] 
k=0
for i in range(1,len(prob_amp_10zPNR_val)):
    if k==0:
        if i%2!=0:
            a=lst3[i-1][0]+0.28
            b=lst3[i-1][1]-0.28
            k=k+1
        else:
            a=lst3[i-1][0]+0.1
            b=lst3[i-1][1]
        
    elif k==1:
        if i%2!=0:
            a=lst3[i-1][0]+0.28
            b=lst3[i-1][1]+0.28
            k=0
        else:
            a=lst3[i-1][0]+0.1
            b=lst3[i-1][1]

    lst3.append([a,b]) 



'''
Plotting probability amplitude for perturbed system
'''

x_coords, y_coords = zip(*lst)
marker_sizes_cond = [p * 2000 for p in prob_amp_10zPNR_cond]  # Adjust the multiplier for appropriate bubble sizes

x_coords2, y_coords2 = zip(*lst2)
marker_sizes_cond_2 = [p * 2000 for p in prob_amp_10zPNR_cond]  # Adjust the multiplier for appropriate bubble sizes

x_coords1, y_coords1 = zip(*lst1)
marker_sizes_val = [p * 2000 for p in prob_amp_10zPNR_val]  # Adjust the multiplier for appropriate bubble sizes

x_coords3, y_coords3 = zip(*lst3)
marker_sizes_val_3 = [p * 2000 for p in prob_amp_10zPNR_val]  # Adjust the multiplier for appropriate bubble sizes


plt.figure(figsize=(6, 6))
plt.scatter(x_coords, y_coords, s=marker_sizes_cond, alpha=0.5, c='red', label='k=0 at upper flat band')
for i in range(len(lst) - 1):
    plt.plot([lst[i][0], lst[i + 1][0]], [lst[i][1], lst[i + 1][1]], color='black')

plt.scatter(x_coords2, y_coords2, s=marker_sizes_cond_2, alpha=0.5, c='red')
for i in range(len(lst2) - 1):
    plt.plot([lst2[i][0], lst2[i + 1][0]], [lst2[i][1], lst2[i + 1][1]], color='black')
    
plt.scatter(x_coords1, y_coords1, s=marker_sizes_val, alpha=0.5, c='blue', label='k=0 at lower flat band')
for i in range(len(lst1) - 1):
    plt.plot([lst1[i][0], lst1[i + 1][0]], [lst1[i][1], lst1[i + 1][1]], color='black')

plt.scatter(x_coords3, y_coords3, s=marker_sizes_val, alpha=0.5, c='blue')
for i in range(len(lst3) - 1):
    plt.plot([lst3[i][0], lst3[i + 1][0]], [lst3[i][1], lst3[i + 1][1]], color='black')

plt.xlabel('X',fontsize = 18)
plt.ylabel('Y',fontsize = 18)
plt.ylim(-2,2)
plt.title('Probability ampltiude for $E_{ext}$ = 0.016',fontsize = 16)
plt.legend(fontsize = 16)
plt.show()



'''
Plotting probability amplitude for unperturbed system
'''

marker_sizes_cond = [p * 2000 for p in prob_amp_10zPNR_cond_unpert]  # Adjust the multiplier for appropriate bubble sizes

marker_sizes_cond_2 = [p * 2000 for p in prob_amp_10zPNR_cond_unpert]  # Adjust the multiplier for appropriate bubble sizes

marker_sizes_val = [p * 2000 for p in prob_amp_10zPNR_val_unpert]  # Adjust the multiplier for appropriate bubble sizes

marker_sizes_val_3 = [p * 2000 for p in prob_amp_10zPNR_val_unpert]  # Adjust the multiplier for appropriate bubble sizes

plt.figure(figsize=(6, 6))
plt.scatter(x_coords, y_coords, s=marker_sizes_cond, alpha=0.5, c='red', label='k=0 at upper flat band')
for i in range(len(lst) - 1):
    plt.plot([lst[i][0], lst[i + 1][0]], [lst[i][1], lst[i + 1][1]], color='black')

plt.scatter(x_coords2, y_coords2, s=marker_sizes_cond_2, alpha=0.5, c='red')
for i in range(len(lst2) - 1):
    plt.plot([lst2[i][0], lst2[i + 1][0]], [lst2[i][1], lst2[i + 1][1]], color='black')
    
plt.scatter(x_coords1, y_coords1, s=marker_sizes_val, alpha=0.5, c='blue', label='k=0 at lower flat band')
for i in range(len(lst1) - 1):
    plt.plot([lst1[i][0], lst1[i + 1][0]], [lst1[i][1], lst1[i + 1][1]], color='black')

plt.scatter(x_coords3, y_coords3, s=marker_sizes_val, alpha=0.5, c='blue')
for i in range(len(lst3) - 1):
    plt.plot([lst3[i][0], lst3[i + 1][0]], [lst3[i][1], lst3[i + 1][1]], color='black')

plt.xlabel('X',fontsize = 18)
plt.ylabel('Y',fontsize = 18)
plt.title('Probability ampltiude for $E_{ext}$ = 0',fontsize = 16)
plt.ylim(-2,2)
plt.legend(fontsize = 16)
plt.show()


'''
In the end, we are plotting 2 different width of zPNR with some specific electric field
'''

color=['blue','red']

Energies_20z_pert=[]
Energies_20z_unpert=[] 
for i in kblst:
    Hmltn_20zPNR = Hmltn(1.22,-3.66,0.205,0.105,0.055,i,20)    
    Eext_20zPNR_pert = Eext_pert_mat(98.15/2,2.164,2.207,0.008,20)  
    n=len(Hmltn_20zPNR)
  
    eigen_lst_20z_unpert,eigen_vec20z_unpert = find_eigen_vec(Hmltn_20zPNR)
    Energies_20z_unpert.append(eigen_lst_20z_unpert)
    
    overall_mat_pert=[[Hmltn_20zPNR[i][j]+Eext_20zPNR_pert[i][j] for j in range(n)]for i in range(n)]    
    eigen_lst_20z_pert,eigen_vec20z_pert = find_eigen_vec(overall_mat_pert)
    Energies_20z_pert.append(eigen_lst_20z_pert)

    E_plot_20z_pert=[] # Modifying Eigenvalue list (Energies to make it plotted)
    for i in range(len(Energies_20z_pert[0])):
        E=[]
        for j in range(len(Energies_20z_pert)):
            E.append(Energies_20z_pert[j][i])
        E_plot_20z_pert.append(E)

    E_plot_20z_unpert=[] # Modifying Eigenvalue list (Energies to make it plotted)
    for i in range(len(Energies_20z_unpert[0])):
        E=[]
        for j in range(len(Energies_20z_unpert)):
            E.append(Energies_20z_unpert[j][i])
        E_plot_20z_unpert.append(E)

plt.figure(figsize=(4, 2))
plt.plot(kblst, E_plot_20z_unpert[19],c=color[0],label='E = 0')
plt.plot(kblst, E_plot_20z_unpert[20],c=color[1])
plt.xlabel('wave vector',fontsize = 18)
plt.ylabel('Energy (eV)',fontsize = 18)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.title('20-zPNR')
plt.xlim(-math.pi,math.pi)
plt.ylim(-0.6,0.1)
plt.legend(fontsize = 12)
plt.show()

plt.figure(figsize=(4, 2))
plt.plot(kblst, E_plot_20z_pert[19],c=color[0],label='E = 0.007')
plt.plot(kblst, E_plot_20z_pert[20],c=color[1])
plt.xlabel('wave vector',fontsize = 18)
plt.ylabel('Energy (eV)',fontsize = 18)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.xlim(-math.pi,math.pi)
plt.ylim(-0.85,0.05)
plt.title('20-zPNR')
plt.legend(fontsize = 12)
plt.show()


Energies_6z_pert=[]
Energies_6z_unpert=[] 
for i in kblst:
    Hmltn_6zPNR = Hmltn(1.22,-3.66,0.205,0.105,0.055,i,6)    
    Eext_6zPNR_pert = Eext_pert_mat(98.15/2,2.164,2.207,0.012,6)  
    n=len(Hmltn_6zPNR)
  
    eigen_lst_6z_unpert,eigen_vec6z_unpert = find_eigen_vec(Hmltn_6zPNR)
    Energies_6z_unpert.append(eigen_lst_6z_unpert)
    
    overall_mat_pert=[[Hmltn_6zPNR[i][j]+Eext_6zPNR_pert[i][j] for j in range(n)]for i in range(n)]    
    eigen_lst_6z_pert,eigen_vec6z_pert = find_eigen_vec(overall_mat_pert)
    Energies_6z_pert.append(eigen_lst_6z_pert)

    E_plot_6z_pert=[] # Modifying Eigenvalue list (Energies to make it plotted)
    for i in range(len(Energies_6z_pert[0])):
        E=[]
        for j in range(len(Energies_6z_pert)):
            E.append(Energies_6z_pert[j][i])
        E_plot_6z_pert.append(E)

    E_plot_6z_unpert=[] # Modifying Eigenvalue list (Energies to make it plotted)
    for i in range(len(Energies_6z_unpert[0])):
        E=[]
        for j in range(len(Energies_6z_unpert)):
            E.append(Energies_6z_unpert[j][i])
        E_plot_6z_unpert.append(E)

plt.figure(figsize=(4, 2))
plt.plot(kblst, E_plot_6z_unpert[5],c=color[0],label='E = 0')
plt.plot(kblst, E_plot_6z_unpert[6],c=color[1])
plt.xlabel('wave vector',fontsize = 18)
plt.ylabel('Energy (eV)',fontsize = 18)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.xlim(-math.pi,math.pi)
plt.ylim(-0.65,0.1)
plt.title('6-zPNR')
plt.legend(fontsize = 12)
plt.show()

plt.figure(figsize=(4, 2))
plt.plot(kblst, E_plot_6z_pert[5],c=color[0],label='E = 0.012')
plt.plot(kblst, E_plot_6z_pert[6],c=color[1])
plt.xlabel('wave vector',fontsize = 18)
plt.ylabel('Energy (eV)',fontsize = 18)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.xlim(-math.pi,math.pi)
plt.ylim(-0.9,0.3)
plt.title('6-zPNR')
plt.legend(fontsize = 12)
plt.show()






    

    
    

    


    

            
    

    
            
    
    


            




    

            
    

    
    


        
    


