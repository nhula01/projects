#import necessary libs
import matplotlib.pyplot as plt
import numpy as np

#unpack data
data1 = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_02/Newton_Raphson.dat')

#unload data
x_1 = data1[:,0]
y_1 = data1[:,1]
z_1 = data1[:,2]

#plot 1
plt.subplot(1, 2, 1)  # 1 row, 2 columns, select the first subplot
plt.plot(x_1, y_1, 'o',label='iterations vs x_values', color='blue')
plt.title('x convergence')
plt.xlabel('iterations [repetitions]')
plt.ylabel('x_values [arb. unit]')
plt.legend()
plt.grid(True)
plt.legend(loc='lower right')

# Plot the second subplot on the right
plt.subplot(1, 2, 2)  # 1 row, 2 columns, select the second subplot
plt.plot(x_1, z_1, 'o',label='iterations vs y_values', color='red')
plt.title('y convergence')
plt.xlabel('iterations [repetitions]')
plt.ylabel('y-axis [arb. unit]')
plt.legend()
plt.grid(True)
plt.legend(loc='upper right')

plt.show()