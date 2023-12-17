import matplotlib.pyplot as plt
import numpy as np

#load data
data1 = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_01/convergence_1.dat')
data10 = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_01/convergence_10.dat')
data100 = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_01/convergence_100.dat')
data1000 = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_01/convergence_1000.dat')
data10000 = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_01/convergence_10000.dat')

#unpack data
x_1 = data1[:,0]
y_1 = data1[:,1]
x_10 = data10[:,0]
y_10 = data10[:,1]
x_100 = data100[:,0]
y_100 = data100[:,1]
x_1000 = data1000[:,0]
y_1000 = data1000[:,1]
x_10000 = data10000[:,0]
y_10000 = data10000[:,1]

plt.plot(x_1,y_1, label='N=1')
plt.plot(x_10,y_10, label='N=10')
plt.plot(x_100,y_100, label='N=100')
plt.plot(x_1000,y_1000, label='N=1000')
plt.plot(x_10000,y_10000, label='N=10000')
plt.xlabel('x [arb. unit]')
plt.ylabel('y [arb. unit]')

legend = plt.legend()
frame = legend.get_frame()
frame.set_alpha(0.8) #set tranparency

plt.title(r'Function')

#plt.legend()
plt.grid(True)

plt.show()