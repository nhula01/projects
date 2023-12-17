#import necessary libs
import matplotlib.pyplot as plt
import numpy as np

#unpack data for two files
data1 = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_02/lambert_function_W1_NR.dat')
data2 = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_02/lambert_function_W0.dat')

#unload (x1,y1), (x2,y2)
x_1 = data1[:,0]
y_1 = data1[:,1]
x_2 = data2[:,0]
y_2 = data2[:,1]

#plot 2 graphs
plt.plot(x_1,y_1,label='W_neg1')
plt.plot(x_2,y_2,label='W_0')

#label x and y axis
plt.xlabel('z [arb. unit]')
plt.ylabel('W [arb. unit]')

#title
plt.title(r'Lambert Function: $we^{w}=z$')


plt.legend(loc='lower right')
plt.grid(True)

plt.show()