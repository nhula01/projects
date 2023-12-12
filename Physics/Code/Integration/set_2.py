import matplotlib.pyplot as plt
import numpy as np

N_points_value = np.loadtxt("/Users/phihung/Documents/PHY/first/homework/hw_04/N_points_vs_value.dat")
n = N_points_value[:,0]
value = N_points_value[:,1]

upper_limit = np.loadtxt("/Users/phihung/Documents/PHY/first/homework/hw_04/upper_limit.dat")
n = upper_limit[:,0]
value = upper_limit[:,1]
plt.plot(n, value, color='y')
plt.xlabel('b upper limit [arb.unit]')
plt.ylabel('value [arb.unit]')
plt.title('value vs b')

"""real = [2]*len(n)

plt.plot(n, real, label='real value', color='r')
plt.plot(n, value,label='spatial resolution', color='y')
plt.xlabel('n points [arb.unit]')
plt.ylabel('value [arb.unit]')
plt.title('value vs n')
plt.ylim(0, 2.5)"""
plt.legend(loc='lower right')
plt.grid(True)

plt.show()