import numpy as np
import time
from joblib import Parallel, delayed
start_time = time.time()


'''
Now we will calculate Transmission by Landeur Buttiker method for square lattice as both device and leads.
'''

# Defining function for calculation of surface green function
def surface_green(alpha,beta,threshold,pos_imag,E): # alpha = hamiltonian for unit cell of leads, beta = coupling matrix of leads, threshold is for convergance of surface green function,
                                                    # pos_imag = positive imaginary number which need to be added for stability, E = Energy of incident electrons
    N = len(alpha)
    E_matrix = np.diag([E + (10**-3)*1j] * N) # Energy matrix
    alpha_2 = np.array(alpha) 
    diff_matrix = E_matrix - alpha_2 

    try: # initial guess of surface green function
        g_s_initial = np.linalg.inv(diff_matrix)
    except np.linalg.LinAlgError:
        print("\nSum matrix is singular, inverse does not exist.")

    g_s = g_s_initial
    beta_L_dagger = np.conj(beta).T
    a=1
    k=1
    while a>(threshold):     
        E_matrix = np.diag([E + (pos_imag)*1j] * N)       
        mult_1 = np.dot(g_s,beta)
        mult_2 = np.dot(beta_L_dagger,mult_1)
        diff_matrix_2 = E_matrix - alpha_2 - mult_2

        try: 
            g_s_new = np.linalg.inv(diff_matrix_2)
        except np.linalg.LinAlgError:
            print("\nSum matrix is singular, inverse does not exist.") 
        
        a = np.sum(np.abs(g_s_new - g_s))
        g_s = g_s_new
        k=k+1
        print(k,E,a)
    return g_s


# Function for calculation of Transmission by Landeur Buttiker formalism

def landauer_buttiker(E,t,N,l): # N = width of ribbon, l = length of device (l = n means 2n layers)

    '''
    left lead
    '''

    alpha_L = [[0 for i in range(N)]for j in range(N)]
    for i in range(len(alpha_L)-1):
        alpha_L[i][i+1] = t
        alpha_L[i+1][i] = t

    beta_L = [[t if i==j else 0 for i in range(N)]for j in range(N)]
    
    g_s_L = surface_green(alpha_L,beta_L,10**-3,10**-3,E)


    '''
    Right lead
    '''       
    g_s_R = g_s_L


    '''
    coupling
    '''
    tau_L = [[0 for i in range(N)]for j in range(l*N)]
    row_L=[0]
    for i in range(N-1):
        row_L.append(row_L[i]+l)
    
    for i in range(len(row_L)):
        tau_L[row_L[i]][i] = t

    tau_R = [[0 for i in range(N)]for j in range(l*N)]
    row_R = [l-1]
    for i in range(N-1):
        row_R.append(row_R[i]+l)
    
    for i in range(len(row_R)):
        tau_R[row_R[i]][i] = t  

    '''
    Self energy for left lead
    '''
    tau_L_trans = np.conj(tau_L).T
    mult_1=np.dot(g_s_L,tau_L_trans)
    sigma_L=np.dot(tau_L,mult_1)
    sigma_L_trans = np.conj(sigma_L).T

    '''
    Self energy for right lead
    '''
    tau_R_trans = np.conj(tau_R).T
    mult_1=np.dot(g_s_R,tau_R_trans)
    sigma_R=np.dot(tau_R,mult_1)
    sigma_R_trans = np.conj(sigma_R).T
    

    '''
    device hamiltonian
    '''
    H_device = [[0 for i in range(l*N)]for j in range(l*N)]
    
    H_device_1 = [[0 for i in range(l*N)]for j in range(l*N)]
    H_device_2 = [[0 for i in range(l*N)]for j in range(l*N)]

    k=l-1
    for i in range(l*N):
        if i!=k:
            H_device_1[i][i+1] = t
        elif i==k:
            k=k+l
    
    for i in range(l*N-l):
        H_device_1[i][i+l] = t

    for i in range(len(H_device_1)):
        for j in range(len(H_device_1)):
            H_device_2[i][j] = np.conj(H_device_1[j][i])



    for i in range(len(H_device)):
        for j in range(len(H_device)):
            H_device[i][j] = H_device_1[i][j] + H_device_2[i][j]


    E_matrix = np.diag([E + (10**(-4))*1j] * l*N) # Energy matrix 
    diff_matrix = E_matrix - np.array(H_device) - np.array(sigma_L) - np.array(sigma_R)
    G_retarted =  np.linalg.inv(diff_matrix)
    G_retarted_trans = np.conj(G_retarted).T

    '''
    Broadening matrices
    '''

    subt_L = np.array(sigma_L) - np.array(sigma_L_trans) 
    gamma_L = np.array(subt_L)*1j

    subt_R = np.array(sigma_R) - np.array(sigma_R_trans) 
    gamma_R = np.array(subt_R)*1j



    '''
    Writing matrix for transmission coefficient
    '''

    mult1 = np.dot(gamma_R,G_retarted_trans)
    mult2 = np.dot(G_retarted,mult1)
    mult3 = np.dot(gamma_L,mult2)

    transmission_coeff=0
    for i in range(len(mult3)):
        transmission_coeff = transmission_coeff + mult3[i][i]
    


    return np.real(transmission_coeff)

Energy = np.linspace(-6,6,300) 


'''
Parallelising the landauer_buttiker function on a list Energy having grid of energy of incident electrons. We are using 8 processors here. (num_processors = 8, in 1201 line) 
'''

num_processors = 8  # Number of processors or cores to use
T_lst = Parallel(n_jobs=num_processors)(delayed(landauer_buttiker)(E,-1.22,10,6) for E in Energy)


import matplotlib.pyplot as plt
plt.plot(Energy, T_lst, color = 'red',linewidth = 2)

# Adding labels and title
plt.xlabel('Energy',fontsize=18)
plt.ylabel('Transmission',fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.show()

# End time
end_time = time.time()

# Calculate and print the elapsed time
elapsed_time = end_time - start_time
print(f"Elapsed time: {elapsed_time} seconds")