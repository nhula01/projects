#import necessary libs
import matplotlib.pyplot as plt
import numpy as np

#unpack data
data1 = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_02/bisectional_-0_1000.dat')
#data2 = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_02/lambert_function_W0.dat')

#unload (x1,y1), (x2,y2), (x3,y3)
x_1 = data1[:,0]
y_1 = data1[:,1]
z_1 = data1[:,2]
#x_2 = data2[:,0]
#y_2 = data2[:,1]
#plot 3 graphs
plt.subplot(1, 2, 1)  # 1 row, 2 columns, select the first subplot
plt.plot(x_1, y_1, 'o',label='iterations vs x_values', color='blue')
plt.title('x convergence')
plt.xlabel('iterations [repetitions]')
plt.ylabel('x_values [arb. unit]')
plt.legend()
plt.grid(True)
plt.legend(loc='upper right')

# Plot the second subplot on the right
plt.subplot(1, 2, 2)  # 1 row, 2 columns, select the second subplot
plt.plot(x_1, z_1, 'o',label='iterations vs y_values', color='red')
plt.title('y convergence')
plt.xlabel('iterations [repetitions]')
plt.ylabel('y-axis [arb. unit]')
plt.legend()
plt.grid(True)
plt.legend(loc='lower right')

plt.show()