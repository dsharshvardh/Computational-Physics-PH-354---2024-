import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import trapz



# Question 7(a)
# I have added an image named as HW_5_Q_7(a)


# Question 7 (b) and (c)

e = 1.602*10**-19  # Electron charge in eV
V0 = 50*e  
a = 10**-11  
x0 = -10**-10
xf = 10**-10
psi_0 = 0.0
hbar = 1.05457*10**-34  # Planck constant  
m = 9.10938*10**-31  # Mass of electron
N = 1000  
h = (xf - x0) / N
xlst = [x0 + i*h for i in range(N+1)]


def V_Harmonic(x): # Harmonic Potential
    return V0*(x**2)/(a**2)

def V_Anharmonic(x): # Anharmonic Potential
    return V0*(x**4)/(a**4) 


def psi(E,potential): # Have to give which potential we want to use (Harmonic or Anharmonic)
    psi = psi_0
    phi = 1.0
    wavefunction = []
    
    for x in xlst:
        wavefunction.append(psi)
        V = potential(x)
        k1_psi = h*phi
        k1_phi = h*(2*m/hbar**2)*(V - E)*psi
        k2_psi = h*(phi + 0.5*k1_phi)
        k2_phi = h*(2*m/hbar**2)*(V - E)*(psi + 0.5*k1_psi)
        k3_psi = h*(phi + 0.5*k2_phi)
        k3_phi = h*(2*m/hbar**2)*(V - E)*(psi + 0.5*k2_psi)
        k4_psi = h*(phi + k3_phi)
        k4_phi = h*(2*m/hbar**2)*(V - E)*(psi + k3_psi)
        psi = psi + (k1_psi + 2*k2_psi + 2*k3_psi + k4_psi)/6
        phi = phi + (k1_phi + 2*k2_phi + 2*k3_phi + k4_phi)/6

    return np.array(wavefunction, float)


def root(E1, E2,potential): # Potential is function for V
    target_accuracy = e/100000
    wavefunction = psi(E1,potential)
    psi2 = wavefunction[N-1]
    while abs(E1-E2) > target_accuracy:
        wavefunction = psi(E2,potential)
        psi1, psi2 = psi2, wavefunction[N-1]
        E1, E2 = E2, E2-psi2*(E2-E1)/(psi2-psi1)

    mod_squared = wavefunction*wavefunction
    integral = h/3*(mod_squared[0] + mod_squared[N//2 - 1] + 4*sum(mod_squared[1:N//2:2]) + 2*sum(mod_squared[0:N//2 + 1:2]))

    return E2 / e, wavefunction / math.sqrt(2 * integral)



E0, psi0 = root(0, 100*e,V_Harmonic)
E1, psi1 = root(200*e, 400*e,V_Harmonic)
E2, psi2 = root(500*e, 700*e,V_Harmonic)


print('Ground_State_Energy of Harmonic oscilliator =', E0, 'eV')
print('First_Excited_State_Energy Harmonic oscilliator =', E1, 'eV')
print('Second_Excited_State_Energy Harmonic oscilliator =', E2, 'eV')



E0, psi0 = root(0, 100*e,V_Anharmonic)
E1, psi1 = root(200*e, 700*e,V_Anharmonic)
E2, psi2 = root(800*e, 1400*e,V_Anharmonic)

print('Ground_State_Energy of Anharmonic oscilliator =', E0, 'eV')
print('First_Excited_State_Energy of Anharmonic oscilliator =', E1, 'eV')
print('Second_Excited_State_Energy of Anharmonic oscilliator =', E2, 'eV')


x_range = slice(N//4, 3*N//4, 1)
plt.plot(xlst[x_range], psi0[x_range], color='b', label="ground_state")
plt.plot(xlst[x_range], psi1[x_range], color='r', label="First_Excited_state")
plt.plot(xlst[x_range], psi2[x_range], color='g', label="Second_Excited_state")
plt.xlabel('x(m)')
plt.ylabel('\u03C8',size=20)
plt.title("Anharmonic oscilliator")
plt.legend()
plt.show()
