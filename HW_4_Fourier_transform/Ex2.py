# Question 2(a)

import numpy as np
import matplotlib.pyplot as plt


piano=np.loadtxt("piano.txt",float)
plt.plot(piano,color='blue')
plt.title("Waveform of Piano",fontname='Times New Roman',fontsize=18)
plt.show()


c_k_piano=np.fft.rfft(piano)
plt.plot(abs(c_k_piano),color='red')
plt.xlim(0,10000)
plt.xlabel("k",fontname='Times New Roman',fontsize=18)
plt.ylabel("$C_{k}$", fontname='Times New Roman',fontsize=18)
plt.title("Fourier Transform of piano Signal",fontname='Times New Roman',fontsize=18)
plt.show()



trumpet=np.loadtxt("trumpet.txt",float)
plt.plot(trumpet,color='blue')
plt.title("Waveform of Trumpet",fontname='Times New Roman',fontsize=18)
plt.show()


c_k_trumpet=np.fft.rfft(trumpet)
plt.plot(abs(c_k_trumpet),color='red')
plt.xlim(0,10000)
plt.xlabel("k",fontname='Times New Roman',fontsize=18)
plt.ylabel("$C_{k}$",fontname='Times New Roman',fontsize=18)
plt.title("Fourier Transform of Trumpet Signal",fontname='Times New Roman',fontsize=18)
plt.show()


'''
Here in waveform plots the scale of the amplitude for the piano is about twice that of the trumpet, which tells the note from the piano is louder. The sharpness of the waveform
suggests that the note played by the piano is spread out over fewer harmonics than that of the trumpet and thus will be a more cleaner note.
In the plot for the Fourier coefficients, we see that there is a high peak from a single frequency and other peaks are very small in the case for the piano whereas the coefficients for
the trumpet are more spread out, as we expected.
Therefore, the harmonic overtones of a trumpet are more prominently created whereas a piano doesn't produce much overtone. The exact frequencies for the same note would
also be different in the two instruments due to a difference in the tuning.
Since the separation between the samples in the waveform is about 1/44100 seconds, the spacing between adjacent k values would be 44100/100000.
'''




# Question 2(b)

print("The frequencies for the peaks of the fourier coefficients of the piano are")

c_k_piano = [i for i in c_k_piano] # converting to list to use undexing
tmax=max(c_k_piano[1:2000])
kmax=c_k_piano.index(tmax)
print(kmax*44100/100000)
tmax=max(c_k_piano[2000:3400])
kmax=c_k_piano.index(tmax)
print(kmax*44100/100000)
tmax=max(c_k_piano[3400:4000])
kmax=c_k_piano.index(tmax)
print(kmax*44100/100000)
tmax=max(c_k_piano[4000:5500])
kmax=c_k_piano.index(tmax)
print(kmax*44100/100000)
tmax=max(c_k_piano[5500:6500])
kmax=c_k_piano.index(tmax)
print(kmax*44100/100000)
tmax=max(c_k_piano[6500:8000])
kmax=c_k_piano.index(tmax)
print(kmax*44100/100000)
print("The approximate frequency for the piano, therefore is 536.34Hz")



print("The frequencies for the peaks of the fourier coefficients of the trumpet are")

c_k_trumpet=[i for i in c_k_trumpet]
tmax=max(c_k_trumpet[1:2000])
kmax=c_k_trumpet.index(tmax)
print(kmax*44100/100000)
tmax=max(c_k_trumpet[2000:3400])
kmax=c_k_trumpet.index(tmax)
print(kmax*44100/100000)
tmax=max(c_k_trumpet[3400:4000])
kmax=c_k_trumpet.index(tmax)
print(kmax*44100/100000)
tmax=max(c_k_trumpet[4000:5500])
kmax=c_k_trumpet.index(tmax)
print(kmax*44100/100000)
tmax=max(c_k_trumpet[5500:6500])
kmax=c_k_trumpet.index(tmax)
print(kmax*44100/100000)
tmax=max(c_k_trumpet[6500:8000])
kmax=c_k_trumpet.index(tmax)
print(kmax*44100/100000) 
print("The approximate frequency for the trumpet, therefore is 522.14Hz")

# These notes are approximately double of the frequency of C(261 Hz) and hence the instruments were playing a note close to middle C(with the pitch of the piano a little
# higher than exactly middle C)