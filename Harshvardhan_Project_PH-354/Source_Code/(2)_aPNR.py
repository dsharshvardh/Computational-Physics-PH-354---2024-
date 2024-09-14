import numpy as np
import sympy as sp
import math
import cmath
import matplotlib.pyplot as plt


'''
In this file, at first we will find band structure of 10-aPNR then we will see the effect of an in-plane external electric field in the band structure
and the change in probability amplitude because of this in-plane electric field. Then we will see an interesting phenomena in which band gap opened a little after 
critical electric field Eext = 0.339 and then again closed.  
'''

'''
For eigen values of a matrix (For clarity of my code, I will use this function just for eigen values for ploting 
of bands and construct another function for both eigen values and eigen vectors).
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
Lets construct a momentum space hamiltonian for a-PNR having width N. To construct these Hamiltonians, we will make Hamiltonians for different hoppings seprately and in the end, we will add them all
and then will run the calculations that we want.
'''

def Hmltn(t1,t2,t3,t4,t5,ka,N):

    H_net=[[0 for i in range(2*N)]for j in range(2*N)]
    
    '''
    t1 hopping
    '''

    H_t1=[]
    H_0=[0 if i<(2*N-1) else -t1 for i in range(2*N)] # defining edge rows in the start
    H_1=[-t1 if i==2 else 0 for i in range(2*N)]
    H_t1.append(H_0)
    H_t1.append(H_1)
    for i in range(2,2*N): # For Bulk
        H_i=[]
        if i==N:
            for j in range(2*N):
                if j==i-1:
                    H_i.append(-t1)
                else:
                    H_i.append(0)
        
        elif i==N+1:
            for j in range(2*N):
                if j==i+1:
                    H_i.append(-t1)
                else:
                    H_i.append(0)
        
        elif i==2*N-1:
            for j in range(2*N):
                if j==i-1:
                    H_i.append(-t1)
                elif j==0:
                    H_i.append(-t1)
                else:
                    H_i.append(0)
        
        else:
            for j in range(2*N):
                if j==i-1 or j==i+1:
                    H_i.append(-t1)
                else:
                    H_i.append(0)
        
        H_t1.append(H_i)
    

    '''
    t2 hopping
    '''
    H_t2=[[0 for i in range(2*N)]for j in range(2*N)]
    H_t2_1=[[0 for i in range(2*N)]for j in range(2*N)]
    H_t2_2=[[0 for i in range(2*N)]for j in range(2*N)]
    if N%2==0:
        H_t2_1[1][0]=-t2
        
        # for t2 inside unit cell
        row=[]
        column=[]
        for i in range(3,N,2):
            row.append(i)
        
        for j in range(2*N-2,N+1,-2):
            column.append(j)
        
        for i in range(len(row)):
            H_t2_1[row[i]][column[i]]=-t2

        #for t2 outside unit cell
        
        row_out=[]
        column_out=[]
        for i in range(2,N+1,2):
            row_out.append(i)
        
        for j in range(2*N-1,N,-2):
            column_out.append(j)

        for i in range(len(row_out)):
            H_t2_1[row_out[i]][column_out[i]]=-t2*cmath.exp(complex(0,ka))
        

        # Using Hermiticity
        for i in range(len(H_t2_1)):
           for j in range(len(H_t2_1)):
                H_t2_2[i][j] = np.conj(H_t2_1[j][i])

        for i in range(len(H_t2)):
            for j in range(len(H_t2)):
                H_t2[i][j]=(H_t2_1[i][j]+H_t2_2[i][j])

    else:
        H_t2_1[1][0]=-t2
        
        # for t2 inside unit cell
        row=[]
        column=[]
        for i in range(3,N+1,2):
            row.append(i)
        
        for j in range(2*N-2,N,-2):
            column.append(j)
        
        for i in range(len(row)):
            H_t2_1[row[i]][column[i]]=-t2

        #for t2 outside unit cell
        
        row_out=[]
        column_out=[]
        for i in range(2,N,2):
            row_out.append(i)
        
        for j in range(2*N-1,N+1,-2):
            column_out.append(j)

        for i in range(len(row_out)):
            H_t2_1[row_out[i]][column_out[i]]=-t2*cmath.exp(complex(0,ka))
        

        # Using Hermiticity
        for i in range(len(H_t2_1)):
           for j in range(len(H_t2_1)):
                H_t2_2[i][j] = np.conj(H_t2_1[j][i])

        for i in range(len(H_t2)):
            for j in range(len(H_t2)):
                H_t2[i][j]=(H_t2_1[i][j]+H_t2_2[i][j])
   
    
    '''
    t3 hopping
    '''
    H_t3=[[0 for i in range(2*N)]for j in range(2*N)]
    H_t3_1=[[0 for i in range(2*N)]for j in range(2*N)]
    H_t3_2=[[0 for i in range(2*N)]for j in range(2*N)]

    for i in range(1,2*N-1):
        if i<N:
            if i%2==0:
                H_t3_1[i][i+1] = -t3*cmath.exp(complex(0,ka))
            else:
                H_t3_1[i][i+1] = -t3*cmath.exp(complex(0,-ka))
            
        elif i==N:
            continue

        elif i>N:
            if i%2==0:
                H_t3_1[i][i+1] = -t3*cmath.exp(complex(0,ka))
            else:
                H_t3_1[i][i+1] = -t3*cmath.exp(complex(0,-ka))

    H_t3_1[0][2*N-1] = -t3*cmath.exp(complex(0,ka))
        
    # Using Hermiticity
    for i in range(len(H_t3_1)):
        for j in range(len(H_t3_1)):
            H_t3_2[i][j] = np.conj(H_t3_1[j][i])

    for i in range(len(H_t3)):
        for j in range(len(H_t3)):
            H_t3[i][j]=(H_t3_1[i][j]+H_t3_2[i][j])
        


    '''
    t4 hopping
    '''
    H_t4=[[0 for i in range(2*N)]for j in range(2*N)]
    H_t4_1=[[0 for i in range(2*N)]for j in range(2*N)]
    H_t4_2=[[0 for i in range(2*N)]for j in range(2*N)]
    row_1=[]  
    column_1=[]
    row_2=[0,]
    column_2=[]
    for i in range(1,N):
        row_1.append(i)
    for j in range(2*N-1,N,-1):
        column_1.append(j)
    
    for i in range(2*N-1,N+1,-1):
        row_2.append(i)
    for j in range(2,N+1):
        column_2.append(j)

    for i in range(len(row_1)):
        H_t4_1[row_1[i]][column_1[i]] = -t4*(1+cmath.exp(complex(0,ka)))
        H_t4_1[row_2[i]][column_2[i]] = -t4*(1+cmath.exp(complex(0,-ka)))
    
    for i in range(len(H_t4_1)):
        for j in range(len(H_t4_1)):
            H_t4_2[i][j] = np.conj(H_t4_1[j][i])

    for i in range(len(H_t4)):
        for j in range(len(H_t4)):
            H_t4[i][j]=(H_t4_1[i][j]+H_t4_2[i][j])
    

    
    '''
    t5 hopping
    '''
    H_t5=[[0 for i in range(2*N)]for j in range(2*N)]
    H_t5_1=[[0 for i in range(2*N)]for j in range(2*N)]
    H_t5_2=[[0 for i in range(2*N)]for j in range(2*N)]
    row_1=[1]
    column_1=[0]
    for i in range(2,N+1):
            row_1.append(i)
    for j in range(2*N-1,N,-1):
            column_1.append(j)
        
    for i in range(len(row_1)):
            if i%2==0:
                H_t5_1[row_1[i]][column_1[i]] = -t5*(cmath.exp(complex(0,ka)))
            else:
                H_t5_1[row_1[i]][column_1[i]] = -t5
        

    row_2=[]
    column_2=[]
    for i in range(1,N-1):
            row_2.append(i)
    for j in range(2*N-2,N,-1):
            column_2.append(j)
        
    for i in range(len(row_2)):
            if i%2==0:
                H_t5_1[row_2[i]][column_2[i]] = -t5
            else:
                H_t5_1[row_2[i]][column_2[i]] = -t5*(cmath.exp(complex(0,ka)))

    row_3=[0]
    column_3=[3]
    for i in range(2*N-1,N+2,-1):
            row_3.append(i)
    for j in range(4,N+1):
            column_3.append(j)

    for i in range(len(row_3)):
            if i%2!=0:
                H_t5_1[row_3[i]][column_3[i]] = -t5*(cmath.exp(complex(0,-ka)))
            else:
                H_t5_1[row_3[i]][column_3[i]] = -t5
        
    for i in range(len(H_t5_1)):
            for j in range(len(H_t5_1)):
                H_t5_2[i][j] = np.conj(H_t5_1[j][i])

    for i in range(len(H_t5)):
            for j in range(len(H_t5)):
                H_t5[i][j]=(H_t5_1[i][j]+H_t5_2[i][j])
    
    '''
    Net matrix
    '''
    for i in range(len(H_net)):
        for j in range(len(H_net)):
            H_net[i][j] = H_t1[i][j] + H_t2[i][j] + H_t3[i][j] + H_t4[i][j] + H_t5[i][j]

    return H_net


'''
Plotting band structure for 8-aPNR. You can see the results for 50-aPNR by changing 8 to 50 below and different hopping ratios which I talked about in my report.
'''

kblst=np.linspace(-np.pi,np.pi,201) # k-path

Energies=[] # Collecting eigen values for k-path grid
for i in kblst:
    eigen_lst,eigen_vec = find_eigen_vec(Hmltn(1.22,-3.66,0.205,0.105,0.055,i,8)) 
    Energies.append(eigen_lst)


E_plot=[] # Modifying Eigenvalue list (Energies to make it plotted)
for i in range(len(Energies[0])):
    E=[]
    for j in range(len(Energies)):
        E.append(Energies[j][i])
    E_plot.append(E)

k_0=int((len(E_plot[0])-1)/2)
N=int(len(E_plot)/2)
N_1=int((len(E_plot)/2)-1)
BG = E_plot[N][k_0] - E_plot[N_1][k_0]
frac = abs(E_plot[N][k_0])/BG


import matplotlib.pyplot as plt
plt.figure(figsize=(2, 6))

for i, y_values in enumerate(E_plot):
    plt.plot(kblst, y_values,color='blue', label=f'List {i+1}')

plt.xlim(-math.pi, math.pi)
plt.xlabel('k',fontsize = 18)
plt.ylabel('Energy (eV)', fontsize = 18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.title('$t_2$ = 3|$t_1$|',fontsize = 16)
plt.show()
    


'''
For effect of in plane external electric field in 8-aPNR
'''

def Eext_pert_mat(Eext,N): # N = width
    mat=[[0 for i in range(2*N)]for j in range(2*N)] 
    dist=0
    e = 1 # Finding energy in eV
    for i in range(2,2*N): # Taking 0th atom as reference
        if i <= N:
            dist = dist + 3.27/2
        
        elif i == N+1:
            dist = dist

        elif i > N+1:
            dist = dist - 3.27/2
        
        mat[i][i] = -e*Eext*dist
    
    return mat  # mat is the perturbation matrix that we have to add with over tight binding Hamiltonian


# t1 = 3|t2|,  8-aPNR

Energies_8a=[] # Collecting eigen values for k-path grid
for i in kblst:
    Hmltn_8aPNR = Hmltn(1.22,-3.66,0.205,0.105,0.025,i,8)     
    Eext_8aPNR = Eext_pert_mat(0.339,8) 
    n=len(Hmltn_8aPNR)
    overall_mat=[[Hmltn_8aPNR[i][j]+Eext_8aPNR[i][j] for j in range(n)]for i in range(n)]    
    eigen_lst_8a,eigen_vec8a = find_eigen_vec(overall_mat)
    Energies_8a.append(eigen_lst_8a)


E_plot_8a=[] # Modifying Eigenvalue list (Energies to make it plotted)
for i in range(len(Energies_8a[0])):
    E=[]
    for j in range(len(Energies_8a)):
        E.append(Energies_8a[j][i])
    E_plot_8a.append(E)

k_0=int((len(E_plot_8a[0])-1)/2)
N=int(len(E_plot_8a)/2)
N_1=int((len(E_plot_8a)/2)-1)
BG_1 = E_plot_8a[N][k_0] - E_plot_8a[N_1][k_0]
frac_1 = frac*BG_1
fac = E_plot_8a[N][k_0] - frac_1

for i in range(len(E_plot_8a)):
    for j in range(len(E_plot_8a[0])):
        E_plot_8a[i][j] = E_plot_8a[i][j] - fac 



# Plotting Band structure for perturbed 8a-PNR 
plt.figure(figsize=(2, 6))
for i, y_values in enumerate(E_plot_8a):
    plt.plot(kblst, y_values,color='blue', label=f'List {i+1}')

plt.xlim(-math.pi, math.pi)
plt.ylim(-8, 8) 
plt.xlabel('k',fontsize = 18)
plt.ylabel('Energy (eV)', fontsize = 18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.title('$E_{ext} = 0.339$',fontsize = 16)
plt.show()



'''
Probability amplitude of perturbed and unperturbed Hamiltonian
'''

Hmltn_8aPNR = Hmltn(1.22,-3.66,0.205,0.105,0.055,0,8)   
Eext_8aPNR = Eext_pert_mat(0.339,8)
n=len(Hmltn_8aPNR)
overall_mat=[[Hmltn_8aPNR[i][j]+Eext_8aPNR[i][j] for j in range(n)]for i in range(n)]    
eigen_lst_8a,eigen_vec8a = find_eigen_vec(overall_mat)
eigen_lst_8a_unpert,eigen_vec8a_unpert = find_eigen_vec(Hmltn_8aPNR)


'''
Perturbed
'''
prob_amp_8aPNR_cond=[]
for i in eigen_vec8a[8]:
    a=modulus_square(i)
    prob_amp_8aPNR_cond.append(a)

prob_amp_8aPNR_val=[]
for i in eigen_vec8a[7]:
    a=modulus_square(i)
    prob_amp_8aPNR_val.append(a)


'''
Unperturbed
'''
prob_amp_8aPNR_cond_unpert=[]
for i in eigen_vec8a_unpert[8]:
    a=modulus_square(i)
    prob_amp_8aPNR_cond_unpert.append(a)

prob_amp_8aPNR_val_unpert=[]
for i in eigen_vec8a_unpert[7]:
    a=modulus_square(i)
    prob_amp_8aPNR_val_unpert.append(a)



# We want to make a bubble plot in which we will make our unit cell for CBM and VBM at k=0 and size of red and blue color bubble will tell us
# probability amplitude of the sites. Here, lst is unit cell for CBM (at k=0) and lst1 is the unit cell for VBM (at k=0).

# lst is unit cell for CBM (at k=0)
lst=[[0,0.9],[0,1.1]]   

for i in range(2,9):
        a=lst[i-1][0] + 0.4
        if i%2!=0:
            b=lst[i-1][1] - 0.3
        else:
            b=lst[i-1][1] + 0.3

        lst.append([a,b])
lst.append([lst[8][0],lst[8][1]-0.8]) 
for i in range(10,16):
        a=lst[i-1][0] - 0.4
        if i%2!=0:
            b=lst[i-1][1] - 0.3
        else:
            b=lst[i-1][1] + 0.3

        lst.append([a,b]) 

lst.append([lst[0][0],lst[0][1]])
prob_amp_8aPNR_cond.append(prob_amp_8aPNR_cond[0])
prob_amp_8aPNR_cond_unpert.append(prob_amp_8aPNR_cond_unpert[0])


# lst1 is unit cell for VBM (at k=0)
    
lst1=[[0,-1.1],[0,-0.9]]   

for i in range(2,9):
        a=lst1[i-1][0] + 0.4
        if i%2!=0:
            b=lst1[i-1][1] - 0.3
        else:
            b=lst1[i-1][1] + 0.3

        lst1.append([a,b])
lst1.append([lst1[8][0],lst1[8][1]-0.8]) 
for i in range(10,16):
        a=lst1[i-1][0] - 0.4
        if i%2!=0:
            b=lst1[i-1][1] - 0.3
        else:
            b=lst1[i-1][1] + 0.3

        lst1.append([a,b]) 

lst1.append([lst1[0][0],lst1[0][1]])
prob_amp_8aPNR_val.append(prob_amp_8aPNR_val[0])
prob_amp_8aPNR_val_unpert.append(prob_amp_8aPNR_val_unpert[0])


'''
Plotting probability amplitude for perturbed system
'''

x_coords, y_coords = zip(*lst)
marker_sizes_cond = [p * 2000 for p in prob_amp_8aPNR_cond]  # Adjust the multiplier for appropriate bubble sizes


x_coords1, y_coords1 = zip(*lst1)
marker_sizes_val = [p * 2000 for p in prob_amp_8aPNR_val]  # Adjust the multiplier for appropriate bubble sizes



plt.figure(figsize=(6, 6))
plt.scatter(x_coords, y_coords, s=marker_sizes_cond, alpha=0.5, c='red', label='k=0 at conduction band minima')
for i in range(len(lst) - 1):
    plt.plot([lst[i][0], lst[i + 1][0]], [lst[i][1], lst[i + 1][1]], color='black')


    
plt.scatter(x_coords1, y_coords1, s=marker_sizes_val, alpha=0.5, c='blue', label='k=0 at valence band maxima')
for i in range(len(lst1) - 1):
    plt.plot([lst1[i][0], lst1[i + 1][0]], [lst1[i][1], lst1[i + 1][1]], color='black')



plt.xlabel('X',fontsize = 18)
plt.ylabel('Y',fontsize = 18)
plt.title('Probability ampltiude for $E_{ext}$ = 0.339',fontsize = 14)
plt.ylim(-2.5,2.5)
plt.legend(fontsize = 12)
plt.show()


'''
Plotting probability amplitude for unperturbed system
'''

marker_sizes_cond = [p * 2000 for p in prob_amp_8aPNR_cond_unpert]  # Adjust the multiplier for appropriate bubble sizes

marker_sizes_val = [p * 2000 for p in prob_amp_8aPNR_val_unpert]  # Adjust the multiplier for appropriate bubble sizes


plt.figure(figsize=(6, 6))
plt.scatter(x_coords, y_coords, s=marker_sizes_cond, alpha=0.5, c='red', label='k=0 at conduction band minima')
for i in range(len(lst) - 1):
    plt.plot([lst[i][0], lst[i + 1][0]], [lst[i][1], lst[i + 1][1]], color='black')

    
plt.scatter(x_coords1, y_coords1, s=marker_sizes_val, alpha=0.5, c='blue', label='k=0 at valence band maxima')
for i in range(len(lst1) - 1):
    plt.plot([lst1[i][0], lst1[i + 1][0]], [lst1[i][1], lst1[i + 1][1]], color='black')


plt.xlabel('X',fontsize = 18)
plt.ylabel('Y',fontsize = 18)
plt.title('Probability ampltiude for $E_{ext}$ = 0',fontsize = 14)
plt.ylim(-2.5,2.5)
plt.legend(fontsize = 12)
plt.show()



'''
Plotting opening of band gap after critical Electric field
'''

E_ext_pert=[0.339,0.406,0.527]
color=['green','red','blue']
lines=['--','-','-.']
for k in range(len(E_ext_pert)):
    Energies_8a=[] # Collecting eigen values for k-path grid
    for i in kblst:
        Hmltn_8aPNR = Hmltn(1.22,-3.66,0.205,0.105,0.025,i,8)   
        Eext_8aPNR = Eext_pert_mat(E_ext_pert[k],8) 
        n=len(Hmltn_8aPNR)
        overall_mat=[[Hmltn_8aPNR[i][j]+Eext_8aPNR[i][j] for j in range(n)]for i in range(n)]    
        eigen_lst_8a,eigen_vec8a = find_eigen_vec(overall_mat)
        Energies_8a.append(eigen_lst_8a)




    E_plot_8a=[] # Modifying Eigenvalue list (Energies to make it plotted)
    for i in range(len(Energies_8a[0])):
        E=[]
        for j in range(len(Energies_8a)):
            E.append(Energies_8a[j][i])
        E_plot_8a.append(E)

    k_0=int((len(E_plot_8a[0])-1)/2)
    N=int(len(E_plot_8a)/2)
    N_1=int((len(E_plot_8a)/2)-1)
    BG_1 = E_plot_8a[N][k_0] - E_plot_8a[N_1][k_0]
    frac_1 = BG_1/2
    fac = E_plot_8a[N][k_0] - frac_1 + 0.12


    for i in range(len(E_plot_8a)):
        for j in range(len(E_plot_8a[0])):
            E_plot_8a[i][j] = E_plot_8a[i][j] - fac


    plt.plot(kblst, E_plot_8a[7],c=color[k],linestyle=lines[k],label=E_ext_pert[k])
    plt.plot(kblst, E_plot_8a[8],c=color[k],linestyle=lines[k])

    # Adding labels and legend
    plt.xlabel('wave vector',fontsize = 18)
    plt.ylabel('Energy (eV)',fontsize = 18)
    plt.xlim(-math.pi,math.pi)
    plt.legend(fontsize = 14)
    # Display the plot
plt.show()




