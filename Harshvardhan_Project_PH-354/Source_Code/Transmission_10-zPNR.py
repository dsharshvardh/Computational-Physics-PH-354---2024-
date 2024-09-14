import numpy as np
import time
from joblib import Parallel, delayed
start_time = time.time()
import matplotlib.pyplot as plt
import math



'''
For effect of in plane external electric field in z-PNR
'''

def Eext_pert_mat(theta,a,b,Eext,N,l): # a = t1 hopping bond distance, b = t2 hopping bond distance, theta is angle between two t1 bonds, N = width
    mat=[[0 for i in range(2*l*N)]for j in range(2*l*N)] 
    dist=0.797284
    e = 1 # Finding energy in eV
    f=0
    for i in range(2*l,2*l*N): # Taking 0th atom as reference

        if f<2*l:
           f=f+1
        
        else:
           f=0
           dist = dist + 0.797284
        
        mat[i][i] = -e*Eext*dist
    
    return mat  # mat is the perturbation matrix that we have to add with over tight binding Hamiltonian




'''
Now we will calculate Transmission by Landeur Buttiker method for 10-zPNR as both device and leads.
'''


def landauer_buttiker(E,N,l,t1,t2,t3,t4,t5):  # Here E is energy, N is width of the nanoribbon and l is the length of the device
    
    '''
    Left lead unit cell hamiltonian matrix, alpha_L
    '''
    alpha_L = [[0 for i in range(2*N)]for j in range(2*N)]


    # t1 and t2
    H_t1_t2 = [] # Hamiltonian for hopping of atoms within the unit cell
    row_0=[]
    for i in range(2*N):
        if i!=1:
            row_0.append(0)
        else:
            row_0.append(-t1)
    
    H_t1_t2.append(row_0)
    
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
        
        H_t1_t2.append(rows)
    
    row_2N=[]
    for i in range(2*N):
        if i!=2*N-2:
            row_2N.append(0)
        else:
            row_2N.append(-t1)
    
    H_t1_t2.append(row_2N)


    # t3 hopping

    H_t3 = [[0 for i in range(2*N)]for j in range(2*N)]
    H_t3_1 = [[0 for i in range(2*N)]for j in range(2*N)]
    H_t3_2 = [[0 for i in range(2*N)]for j in range(2*N)]
    for i in range(2*N-4): # I want edit rows till 2*N-5
        if i%2==1:
            H_t3_1[i][i+3] = -t3
        else:
            continue

    
    for i in range(len(H_t3_1)):
        for j in range(len(H_t3_1)):
            H_t3_2[i][j]=np.conj(H_t3_1[j][i])

    for i in range(len(H_t3)):
        for j in range(len(H_t3)):
            H_t3[i][j] = H_t3_1[i][j] + H_t3_2[i][j]

    

    # t4 hopping
            
    H_t4 = [[0 for i in range(2*N)]for j in range(2*N)]
    H_t4_1 = [[0 for i in range(2*N)]for j in range(2*N)]
    H_t4_2 = [[0 for i in range(2*N)]for j in range(2*N)]

    for i in range(2*N-2):
        H_t4_1[i][i+2] = -t4

    for i in range(len(H_t4_1)):
        for j in range(len(H_t4_1)):
            H_t4_2[i][j]=np.conj(H_t4_1[j][i])

    for i in range(len(H_t4)):
        for j in range(len(H_t4)):
            H_t4[i][j] = H_t4_1[i][j] + H_t4_2[i][j]

    
    # t5 hopping
            
    H_t5 = [[0 for i in range(2*N)]for j in range(2*N)]
    H_t5_1 = [[0 for i in range(2*N)]for j in range(2*N)]
    H_t5_2 = [[0 for i in range(2*N)]for j in range(2*N)]

    for i in range(2*N-3): # I want edit rows till 2*N-4
        if i%2==0:
            H_t5_1[i][i+3] = -t5
        else:
            continue

    
    for i in range(len(H_t5_1)):
        for j in range(len(H_t5_1)):
            H_t5_2[i][j]=np.conj(H_t5_1[j][i])

    for i in range(len(H_t5)):
        for j in range(len(H_t5)):
            H_t5[i][j] = H_t5_1[i][j] + H_t5_2[i][j]


    
   
    for i in range(len(alpha_L)):
        for j in range(len(alpha_L)):
            alpha_L[i][j] = H_t1_t2[i][j] + H_t3[i][j] + H_t4[i][j] + H_t5[i][j]


    '''
    Left lead unit cell coupling matrix, beta_L
    '''
    beta_L = [[0 for i in range(2*N)]for j in range(2*N)]
    

    # t1

    H_t1 = [[0 for i in range(2*N)]for j in range(2*N)]
    row=[]
    k=0
    if N%2==0:
        for i in range(2*N-1):
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k
    
    else:
        for i in range(2*N):
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k

    for i in range(len(row)):
        if i%2==0:
            H_t1[row[i]][row[i]-1] = -t1
        else:
            H_t1[row[i]][row[i]+1] = -t1

    
    # t2 hopping is not present out of the unit cell
            
    
    #t3
    
    H_t3 = [[0 for i in range(2*N)]for j in range(2*N)]

    for i in range(4,2*N,4):
        if i<(2*N-1):
            H_t3[i-3][i] = -t3

    for i in range(6,2*N,4):
        if i<(2*N-1):
            H_t3[i][i-3] = -t3
    


    # t4

    H_t4 = [[0 for i in range(2*N)]for j in range(2*N)]
    k=0
    row_t4=[]
    if N%2==0:
        for i in range(2*N-1):
            if k==1 or k==2:
                row_t4.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k
    
    else:
        for i in range(2*N):
            if k==1 or k==2:
                row_t4.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k
    
    for i in range(len(row_t4)):
        if i==0:
            H_t4[row_t4[i]][row_t4[i]+2] = -t4
        elif i==(len(row_t4)-1):
            H_t4[row_t4[i]][row_t4[i]-2] = -t4
        else:
            H_t4[row_t4[i]][row_t4[i]+2] = -t4
            H_t4[row_t4[i]][row_t4[i]-2] = -t4

    

    # t5
    H_t5 = [[0 for i in range(2*N)]for j in range(2*N)]
    H_t5_1 = [[0 for i in range(2*N)]for j in range(2*N)]
    H_t5_2 = [[0 for i in range(2*N)]for j in range(2*N)]

    for i in range(1,2*N-2,2):
        H_t5_1[i][i+1] = -t5
    
    for i in range(len(H_t5_1)):
        for j in range(len(H_t5_1)):
            H_t5_2[i][j]=np.conj(H_t5_1[j][i])

    for i in range(len(H_t5)):
        for j in range(len(H_t5)):
            H_t5[i][j] = H_t5_1[i][j] + H_t5_2[i][j]


    for i in range(len(beta_L)):
        for j in range(len(beta_L)):
            beta_L[i][j] = H_t1[i][j] + H_t3[i][j] + H_t4[i][j] + H_t5[i][j]


    '''
    Now we will find surface green function of left lead
    '''
    
    E_matrix = np.diag([E + (10**-4)*1j] * 2*N) # Energy matrix
    alpha_L_2 = np.array(alpha_L) 
    diff_matrix = E_matrix - alpha_L_2 

    try: # initial guess of surface green function
        g_s_L_initial = np.linalg.inv(diff_matrix)
    except np.linalg.LinAlgError:
        print("\nSum matrix is singular, inverse does not exist.")

    g_s_L = g_s_L_initial
    beta_L_dagger = np.conj(beta_L).T
    k=1
    a=1
    while a > (10**-2):      
        E_matrix = np.diag([E + ((10**(-4)))*(1j)] * 2*N) # Energy matrix            
        diff_matrix_2 = E_matrix - alpha_L_2 - g_s_L
        try: 
            g_s_L_new = np.linalg.inv(diff_matrix_2)
        except np.linalg.LinAlgError:
            print("\nSum matrix is singular, inverse does not exist.") 
        mult_1 = np.dot(g_s_L_new,beta_L_dagger)
        g_s_L_new = np.dot(beta_L,mult_1)
        a = np.sum(np.abs(g_s_L_new - g_s_L))
        g_s_L = g_s_L_new        
        k=k+1
        print(k,"zeroth_L",E,a)
    
    g_s_L = np.linalg.inv(E_matrix - alpha_L_2 - g_s_L)


    '''
    Right lead unit cell hamiltonian matrix, alpha_R
    '''

    alpha_R = alpha_L  # Onsite energies and hoppings inside the unit cell will be same as in left lead


    '''
    Right lead unit cell coupling matrix, beta_R
    '''

    beta_R = [[0 for i in range(2*N)]for j in range(2*N)]  

    # t1

    H_t1 = [[0 for i in range(2*N)]for j in range(2*N)]
    H_t1[0][1] = -t1
    row=[]
    k=0
    if N%2==0:
        for i in range(2,2*N):
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k
    
    else:
        for i in range(2,2*N-1):
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k

    for i in range(len(row)):
        if i%2==0:
            H_t1[row[i]][row[i]-1] = -t1
        else:
            H_t1[row[i]][row[i]+1] = -t1

    
    # t2 hopping is not present out of the unit cell
             
    #t3
    
    H_t3 = [[0 for i in range(2*N)]for j in range(2*N)]

    for i in range(6,2*N,4):
        if i<(2*N-1):
            H_t3[i-3][i] = -t3

    for i in range(4,2*N,4):
        if i<(2*N-1):
            H_t3[i][i-3] = -t3


    # t4

    H_t4 = [[0 for i in range(2*N)]for j in range(2*N)]

    row_t4=[0]
    k=0
    if N%2==0:
        for i in range(2,2*N):
            if k==1 or k==2:
                row_t4.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k
    
    else:
        for i in range(2,2*N-1):
            if k==1 or k==2:
                row_t4.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k
    
    for i in range(len(row_t4)):
        if i==0:
            H_t4[row_t4[i]][row_t4[i]+2] = -t4
        elif i==(len(row_t4)-1):
            H_t4[row_t4[i]][row_t4[i]-2] = -t4
        else:
            H_t4[row_t4[i]][row_t4[i]+2] = -t4
            H_t4[row_t4[i]][row_t4[i]-2] = -t4
   

    # t5
    H_t5 = [[0 for i in range(2*N)]for j in range(2*N)]
    H_t5_1 = [[0 for i in range(2*N)]for j in range(2*N)]
    H_t5_2 = [[0 for i in range(2*N)]for j in range(2*N)]

    for i in range(1,2*N-2,2):
        H_t5_1[i][i+1] = -t5
    
    for i in range(len(H_t5_1)):
        for j in range(len(H_t5_1)):
            H_t5_2[i][j]=np.conj(H_t5_1[j][i])

    for i in range(len(H_t5)):
        for j in range(len(H_t5)):
            H_t5[i][j] = H_t5_1[i][j] + H_t5_2[i][j]


    for i in range(len(beta_R)):
        for j in range(len(beta_R)):
            beta_R[i][j] = H_t1[i][j] + H_t3[i][j] + H_t4[i][j] + H_t5[i][j]



    '''
    Now we will find surface green function of right lead
    '''
    E_matrix = np.diag([E + (10**-4)*1j] * 2*N) # Energy matrix
    alpha_R_2 = np.array(alpha_R) 
    diff_matrix = E_matrix - alpha_R_2 

    try: # initial guess of surface green function
        g_s_R_initial = np.linalg.inv(diff_matrix)
    except np.linalg.LinAlgError:
        print("\nSum matrix is singular, inverse does not exist.")

    g_s_R = g_s_R_initial
    beta_R_dagger = np.conj(beta_R).T
    k=1
    a=1
    while a>(10**-2):      
        E_matrix = np.diag([E + (10**(-4))*1j] * 2*N)        
        diff_matrix_2 = E_matrix - alpha_R_2 - g_s_R

        try: 
            g_s_R_new = np.linalg.inv(diff_matrix_2)
        except np.linalg.LinAlgError:
            print("\nSum matrix is singular, inverse does not exist.") 
        
        mult_1 = np.dot(g_s_R_new,beta_R_dagger)
        g_s_R_new =  np.dot(beta_R,mult_1)

        a = np.sum(np.abs(g_s_R_new - g_s_R))
        g_s_R = g_s_R_new
        k=k+1
        print(k,"zeroth_R",E,a)

    g_s_R = np.linalg.inv(E_matrix - alpha_R_2 - g_s_R)


    '''
    Coupling matrix of left lead and device
    '''
    
    tau_L = [[0 for i in range(2*l*N)]for j in range(2*N)]
    
    
    # t1
    
    tau_t1 = [[0 for i in range(2*l*N)]for j in range(2*N)]
    row=[0]
    k=0
    if N%2==0:
        for i in range(2,2*N):
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k
    
    else:
        for i in range(2,2*N-1):
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k
    

    column=[]
    a=len(row)
    i=0
    j=0
    while i!=a:
        column.append(j)
        j=j+2*l
        i=i+1
    
    for i in range(len(row)):
        tau_t1[row[i]][column[i]] = -t1
    

    # t2 is not there in coupling
    

    # t3
    
    tau_t3 = [[0 for i in range(2*l*N)]for j in range(2*N)]
    row=[]
    k=0
    if N%2==0:
        for i in range(2,2*N-1):
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k

    else:
        for i in range(2,2*N-1):
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k

    column=[]
    j=0
    while j!=2*l*N:
        column.append(j)
        j=j+2*l   

    for i in range(len(row)):
        if i%2==0:
            tau_t3[row[i]][column[i+3]] = -t3
        else:
            tau_t3[row[i]][column[i-1]] = -t3


    # t4

    tau_t4 = [[0 for i in range(2*l*N)]for j in range(2*N)]
    row=[0]
    k=0
    if N%2==0:
        for i in range(2,2*N):
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k
    
    else:
        for i in range(2,2*N-1):
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k
    

    column=[]
    a=len(row)
    i=0
    j=0
    while i!=a:
        column.append(j)
        j=j+2*l
        i=i+1   

    for i in range(len(row)):
        if i==0:
            tau_t4[row[i]][column[i+1]] = -t4
        elif i==(len(row)-1):
            tau_t4[row[i]][column[i-1]] = -t4
        else:
            tau_t4[row[i]][column[i+1]] = -t4
            tau_t4[row[i]][column[i-1]] = -t4
    
    # t5

    tau_t5 = [[0 for i in range(2*l*N)]for j in range(2*N)]
    row=[]
    k=0
    if N%2==0:
        for i in range(0,2*N-1):
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k

    else:
        for i in range(0,2*N):
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k


    column=[]
    a=len(row)
    i=0
    j=0
    while i!=a:
        column.append(j)
        j=j+2*l
        i=i+1

    for i in range(len(row)):
        if i%2==0:
            tau_t5[row[i]][column[i+1]] = -t5
        else:
            tau_t5[row[i]][column[i-1]] = -t5

    row=[]
    k=0
    if N%2==0:
        for i in range(2,2*N-1):
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k

    else:
        for i in range(2,2*N-1):
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k

    column=[]
    a=len(row)
    i=0
    j=2*l+1
    while i!=a:
        column.append(j)
        j=j+2*l
        i=i+1

    for i in range(len(row)):
        if i%2==0:
            tau_t5[row[i]][column[i+1]] = -t5
        else:
            tau_t5[row[i]][column[i-1]] = -t5

    

    for i in range(len(tau_L)):
        for j in range(len(tau_L[0])):
            tau_L[i][j] = tau_t1[i][j] + tau_t3[i][j] + tau_t4[i][j] + tau_t5[i][j]



    '''
    Coupling matrix of Right lead and device
    '''
    
    tau_R = [[0 for i in range(2*l*N)]for j in range(2*N)]

    # t1

    tau_t1 = [[0 for i in range(2*l*N)]for j in range(2*N)]
    row=[]
    k=0
    if N%2==0:
        for i in range(0,2*N-1):
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k
    
    else:
        for i in range(0,2*N):
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k
    

    column=[]
    a=len(row)
    i=0
    j=2*l-1
    while i!=a:
        column.append(j)
        j=j+2*l
        i=i+1
    
    for i in range(len(row)):
        tau_t1[row[i]][column[i]] = -t1

    # t2 is not possible

    # t3

    tau_t3 = [[0 for i in range(2*l*N)]for j in range(2*N)]
    row=[]
    k=0
    if N%2==0:
        for i in range(0,2*N-1):
            
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k

    else:
        for i in range(0,2*N):
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k


    column=[]
    a=len(row)
    i=0
    j=2*l-1
    while i!=a:
        column.append(j)
        j=j+2*l
        i=i+1

    for i in range(len(row)):
        if i%2==0:
            if i+2 < len(column):
                tau_t3[row[i]][column[i+2]] = -t3
            else:
                continue
        
        else:
            if i-2 > 0:
                tau_t3[row[i]][column[i-2]] = -t3   
            else:
                continue

        
    # t4
    tau_t4 = [[0 for i in range(2*l*N)]for j in range(2*N)]
    row=[]
    k=0
    if N%2==0:
        for i in range(0,2*N-1):
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k
    
    else:
        for i in range(0,2*N):
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k
    

    column=[]
    a=len(row)
    i=0
    j=2*l-1
    while i!=a:
        column.append(j)
        j=j+2*l
        i=i+1       
    

    for i in range(len(row)):
        if i==0:
            tau_t4[row[i]][column[i+1]] = -t4
        elif i==(len(row)-1):
            tau_t4[row[i]][column[i-1]] = -t4
        else:
            tau_t4[row[i]][column[i+1]] = -t4
            tau_t4[row[i]][column[i-1]] = -t4

    # t5
    
    tau_t5 = [[0 for i in range(2*l*N)]for j in range(2*N)]
    row=[]
    k=0
    if N%2==0:
        for i in range(0,2*N-1):
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k

    else:
        for i in range(0,2*N):
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k


    column=[]
    a=len(row)
    i=0
    j=2*l-2
    while i!=a:
        column.append(j)
        j=j+2*l
        i=i+1

    for i in range(len(row)):
        if i%2==0:
            tau_t5[row[i]][column[i+1]] = -t5

        else:
            tau_t5[row[i]][column[i-1]] = -t5

    row=[]
    k=0
    if N%2==0:
        for i in range(2,2*N-1):
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k

    else:
        for i in range(2,2*N-1):
            if k==1 or k==2:
                row.append(i)
            k=k+1
            
            if k==4:
                k=0
            else:
                k=k


    column=[]
    a=len(row)
    i=0
    j=4*l-1
    while i!=a:
        column.append(j)
        j=j+2*l
        i=i+1
    for i in range(len(row)):
        if i%2==0:
            tau_t5[row[i]][column[i+1]] = -t5
        else:
            tau_t5[row[i]][column[i-1]] = -t5


    for i in range(len(tau_R)):
        for j in range(len(tau_R[0])):
            tau_R[i][j] = tau_t1[i][j] + tau_t3[i][j] + tau_t4[i][j] + tau_t5[i][j]

    
    '''
    Self energy for left lead
    '''
    tau_L_trans = np.conj(tau_L).T
    mult_1=np.dot(g_s_L,tau_L)
    sigma_L=np.dot(tau_L_trans,mult_1)
    sigma_L_trans = np.conj(sigma_L).T

    '''
    Self energy for right lead
    '''
    tau_R_trans = np.conj(tau_R).T
    mult_1=np.dot(g_s_R,tau_R)
    sigma_R=np.dot(tau_R_trans,mult_1)
    sigma_R_trans = np.conj(sigma_R).T


    '''
    Hamiltonian for the device
    '''

    # t1
    H_device_t1_1=[[0 for i in range(2*l*N)]for j in range(2*l*N)]
    H_device_t1_2=[[0 for i in range(2*l*N)]for j in range(2*l*N)]
    H_device_t1=[[0 for i in range(2*l*N)]for j in range(2*l*N)]
    k=0
    for i in range(2*l*N):
        if k!=2*l-1:
            H_device_t1_1[i][i+1] = -t1
            k=k+1
        else:
            k=0
        

    for i in range(len(H_device_t1_1)):
        for j in range(len(H_device_t1_1[0])):
            H_device_t1_2[i][j] = H_device_t1_1[j][i]
    
    for i in range(len(H_device_t1)):
        for j in range(len(H_device_t1)):
            H_device_t1[i][j] = H_device_t1_1[i][j] + H_device_t1_2[i][j]


    # t2

    H_device_t2_1=[[0 for i in range(2*l*N)]for j in range(2*l*N)]
    H_device_t2_2=[[0 for i in range(2*l*N)]for j in range(2*l*N)]
    H_device_t2=[[0 for i in range(2*l*N)]for j in range(2*l*N)]       
    
    k=0
    j=0
    for i in range(2*l*N - 2*l):
        if j <= (2*l-1):
            k=0
        else:
            k=1
        
        if i%2==k:
            H_device_t2_1[i][i+2*l] = -t2
        
        if j == (4*l-1):
            j=0
        else:
            j=j+1


    
    for i in range(len(H_device_t2_1)):
        for j in range(len(H_device_t2_1)):
            H_device_t2_2[i][j] = H_device_t2_1[j][i]

    for i in range(len(H_device_t2)):
        for j in range(len(H_device_t2)):
            H_device_t2[i][j] = H_device_t2_1[i][j] + H_device_t2_2[i][j]



    # t3

    H_device_t3_1=[[0 for i in range(2*l*N)]for j in range(2*l*N)]
    H_device_t3_2=[[0 for i in range(2*l*N)]for j in range(2*l*N)]
    H_device_t3=[[0 for i in range(2*l*N)]for j in range(2*l*N)]       
    k=0
    j=0
    row=[]
    for i in range(2*l*N - 4*l):
        if j <= (2*l-1):
            k=0
        else:
            k=1
        
        if i%2==k:
            row.append(i)
        
        if j == (4*l-1):
            j=0
        else:
            j=j+1

    a=0
    b=2*l-1
 
    for i in row:
        if i>a:
            a=a+2*l
        if i>b:
            b=b+2*l
    
        if i != a and i!= b:
            H_device_t3_1[i][i+4*l-1] = -t3
            H_device_t3_1[i][i+4*l+1] = -t3
        elif i == a:
            H_device_t3_1[i][i+4*l+1] = -t3
        elif i == b:
            H_device_t3_1[i][i+4*l-1] = -t3

    for i in range(len(H_device_t3_1)):
        for j in range(len(H_device_t3_1)):
            H_device_t3_2[i][j] = H_device_t3_1[j][i]

    for i in range(len(H_device_t3)):
        for j in range(len(H_device_t3)):
            H_device_t3[i][j] = H_device_t3_1[i][j] + H_device_t3_2[i][j]



    
    # t4

    H_device_t4_1=[[0 for i in range(2*l*N)]for j in range(2*l*N)]
    H_device_t4_2=[[0 for i in range(2*l*N)]for j in range(2*l*N)]
    H_device_t4=[[0 for i in range(2*l*N)]for j in range(2*l*N)]          
    k=0
    for i in range(2*l*N-2*l):
        if k!=0 and k!=2*l-1:
            H_device_t4_1[i][i+2*l-1] = -t4
            H_device_t4_1[i][i+2*l+1] = -t4
            k=k+1
        elif k==0:
            H_device_t4_1[i][i+2*l+1] = -t4
            k=k+1
        else:
            H_device_t4_1[i][i+2*l-1] = -t4
            k=0

    for i in range(len(H_device_t4_1)):
        for j in range(len(H_device_t4_1)):
            H_device_t4_2[i][j] = H_device_t4_1[j][i]

    for i in range(len(H_device_t4)):
        for j in range(len(H_device_t4)):
            H_device_t4[i][j] = H_device_t4_1[i][j] + H_device_t4_2[i][j]

    
    # t5
    H_device_t5_1=[[0 for i in range(2*l*N)]for j in range(2*l*N)]
    H_device_t5_2=[[0 for i in range(2*l*N)]for j in range(2*l*N)]
    H_device_t5=[[0 for i in range(2*l*N)]for j in range(2*l*N)]       

    k=0
    j=0
    for i in range(2*l*N - 2*l):
        if j <= (2*l-1):
            k=1
        else:
            k=0
        
        if i%2==k:
            H_device_t5_1[i][i+2*l] = -t5
        
        if j == (4*l-1):
            j=0
        else:
            j=j+1

    k=0
    j=0
    row=[]
    for i in range(2*l*N - 2*l):
        if j <= (2*l-1):
            k=0
        else:
            k=1
        
        if i%2==k:
            row.append(i)       
        
        if j == (4*l-1):
            j=0
        else:
            j=j+1

    f=0
    
    for i in range(len(row)):
        if f==0:
            H_device_t5_1[row[i]][row[i]+2*l+2] = -t5
            f=f+1
        elif f!=0 and f!=l-1:
            H_device_t5_1[row[i]][row[i]+2*l+2] = -t5
            H_device_t5_1[row[i]][row[i]+2*l-2] = -t5
            f=f+1
        elif f==l-1:
            H_device_t5_1[row[i]][row[i]+2*l-2] = -t5
            f=0
    
    for i in range(len(H_device_t5_1)):
        for j in range(len(H_device_t5_1)):
            H_device_t5_2[i][j] = H_device_t5_1[j][i]

    for i in range(len(H_device_t5)):
        for j in range(len(H_device_t5)):
            H_device_t5[i][j] = H_device_t5_1[i][j] + H_device_t5_2[i][j]



    H_device=[[0 for i in range(2*l*N)]for j in range(2*l*N)]  
    for i in range(len(H_device)):
        for j in range(len(H_device)):
            H_device[i][j] = H_device_t1[i][j] + H_device_t2[i][j] + H_device_t4[i][j] + H_device_t3[i][j] + H_device_t5[i][j]  

    


    E_matrix = np.diag([E +  (10**(-4))*1j] * 2*l*N) # Energy matrix 
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

    # Transmission = Trace(mult3)
    transmission_coeff=0
    for i in range(len(mult3)):
        transmission_coeff = transmission_coeff + mult3[i][i]    
    

    return np.real(transmission_coeff)            


Energy = np.linspace(-6,6,300)  


'''
Parallelising the landauer_buttiker function on a list Energy having grid of energy of incident electrons. We are using 8 processors here. (num_processors = 8, in 1201 line) 
'''

num_processors = 8  # Number of processors or cores to use
T_lst = Parallel(n_jobs=num_processors)(delayed(landauer_buttiker)(E,10,3,1.22,-3.66,0,0.105,0) for E in Energy)


plt.plot(Energy, T_lst,color='red',linewidth=2)

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