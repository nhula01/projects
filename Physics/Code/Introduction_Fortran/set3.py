import matplotlib.pyplot as plt
import numpy as np

#unload data
data = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_01/gregory_convergence.dat')

# unpack value
n = data[:,0]
pi = data[:,1]

#define y=pi function
def pi_function(n):
    return n*0 + np.pi
# record on k
k = pi_function(n)

# plot the repetition vs 
plt.plot(n,pi,label='Gregory series')
plt.xlabel('N [repetitions]')
plt.ylabel('Value [arb. unit]')

plt.plot(n,k, color='red', label='pi')

plt.title(r'Gregory Series')

plt.legend()
plt.grid(True)

plt.show()