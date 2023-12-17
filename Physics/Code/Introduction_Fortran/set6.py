import matplotlib.pyplot as plt
import numpy as np

#load data
x_data = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_01/x_t.dat')
y_data = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_01/y_t.dat')
v_data = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_01/vy_t.dat')

#plot graph
plt.plot(x_data[:,2],y_data[:,2], label='exact')
plt.plot(x_data[:,1],y_data[:,1],label='approximation')
plt.title(r'Projectile motion at dt=.05')

plt.xlabel('Distance x [m]')
plt.ylabel('Height y [m]')

plt.grid(True)
plt.legend()
plt.show()