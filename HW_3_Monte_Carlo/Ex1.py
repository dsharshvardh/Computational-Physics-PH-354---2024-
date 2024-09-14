import numpy as np
# Life times in seconds

t_bi = 46*60
t_tl = 2.2*60
t_pb = 3.3*60

# Number of atoms initially 

bi = 10000
tl = 0
pb = 0
bi2 = 0

# probabilities of decay
t=1
bi_prob = 1-(2**(-t/t_bi))
tl_prob = 1-(2**(-t/t_tl))
pb_prob = 1-(2**(-t/t_pb))

# lists for plotting
t_lst = []
bi_lst = []
pb_lst = []
tl_lst = []
bi2_lst = []

for time in range(20001):
    pb = pb - pb*pb_prob
    bi2 = bi2 + pb*pb_prob
    tl = tl - tl*tl_prob
    pb = pb + tl*tl_prob

    bi = bi-bi*bi_prob # shows decay of Bi213
    rand_num = np.random.uniform(0,1)
    
    # Formation and decay of Tl209 and Pb209
    if rand_num <= 0.0209:
        tl = tl + bi*bi_prob 
    else:
        pb = pb + bi*bi_prob
    
    t_lst.append(time)
    bi_lst.append(bi)
    pb_lst.append(pb)
    tl_lst.append(tl)
    bi2_lst.append(bi2)

import matplotlib.pyplot as plt

plt.plot(t_lst, bi_lst, label='$^{213}$Bi', color='blue')
plt.plot(t_lst, pb_lst, label='$^{209}$Pb', color='green')
plt.plot(t_lst, tl_lst, label='$^{209}$Tl', color='red')
plt.plot(t_lst, bi2_lst, label='$^{209}$Bi', color='purple')

plt.xlabel('Time(second)', fontsize=18, fontname='Times New Roman')
plt.ylabel('Number of atoms', fontsize=18, fontname='Times New Roman')

plt.legend(fontsize=14)
plt.xticks(fontsize=16, fontname='Times New Roman')
plt.yticks(fontsize=16, fontname='Times New Roman')
plt.show()

