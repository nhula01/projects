#import necessary libs
import matplotlib.pyplot as plt
import numpy as np

#unpack data
data1 = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_01/cartesian_coordinates_1.dat')
data2 = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_01/cartesian_coordinates_2.dat')
data3 = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_01/cartesian_coordinates_3.dat')
#unload (x1,y1), (x2,y2), (x3,y3)
x_1 = data1[:,0]
y_1 = data1[:,1]
x_2 = data2[:,0]
y_2 = data2[:,1]
x_3 = data3[:,0]
y_3 = data3[:,1]

#plot 3 graphs
plt.plot(x_1,y_1,label='a=2')
plt.plot(x_2,y_2,label='a=3')
plt.plot(x_3,y_3,label='a=4')

#label x and y axis
plt.xlabel('x [arb. unit]')
plt.ylabel('y [arb. unit]')

#title
plt.title(r'Function: $x^{3} + y^{3} - 3xy = 0$')


plt.legend(loc='lower left')
plt.grid(True)

plt.show()