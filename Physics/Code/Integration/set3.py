import matplotlib.pyplot as plt
import numpy as np

exp_data = np.loadtxt("/Users/phihung/Documents/PHY/first/homework/hw_04/exp.dat")
x_data = np.loadtxt("/Users/phihung/Documents/PHY/first/homework/hw_04/1_x.dat")

#plot expo
"""
exp_space = exp_data[:,0]
exp_value = exp_data[:,1]

plt.plot(exp_space, exp_value,label='spatial resolution', color='y')
plt.xlabel('n [points]')
plt.ylabel('value')
plt.title('exp')
plt.legend(loc='upper left')
plt.grid(True)
"""
#plot 1_x
x_space = x_data[:,0]
x_value = x_data[:,1]

plt.plot(x_space, x_value,label='spatial resolution', color='y')
plt.xlabel('n [points]')
plt.ylabel('value')
plt.title('1_x')
plt.legend(loc='upper left')
plt.grid(True)

plt.show()