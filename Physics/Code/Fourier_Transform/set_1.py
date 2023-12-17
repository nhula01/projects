import matplotlib.pyplot as plt
import numpy as np

# plot negative first
data_neg4 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/dft_abs_neg_window4.dat')
data_neg6 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/dft_abs_neg_window6.dat')
data_neg8 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/dft_abs_neg_window8.dat')
data_neg10 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/dft_abs_neg_window10.dat')


plt.figure(1)
plt.plot(data_neg4[:,0], data_neg4[:,1], label="a=4")
plt.plot(data_neg6[:,0], data_neg6[:,1], label="a=6")
plt.plot(data_neg8[:,0], data_neg8[:,1], label="a=8")
plt.plot(data_neg10[:,0], data_neg10[:,1], label="a=10")

plt.xlabel(r'$\omega$ [Hz]')
plt.ylabel(r'$\left|f(\omega)\right|$')
plt.title(r'$\omega vs \left|f(\omega)\right|$ in negative frequencies')
plt.legend(loc='upper left')
plt.grid(True)

data_pos4 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/dft_abs_pos_window4.dat')
data_pos6 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/dft_abs_pos_window6.dat')
data_pos8 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/dft_abs_pos_window8.dat')
data_pos10 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/dft_abs_pos_window10.dat')
plt.figure(2)
plt.plot(data_pos4[:,0], data_pos4[:,1], label="a=4")
plt.plot(data_pos6[:,0], data_pos6[:,1], label="a=6")
plt.plot(data_pos8[:,0], data_pos8[:,1], label="a=8")
plt.plot(data_pos10[:,0], data_pos10[:,1], label="a=10")

plt.xlabel(r'$\omega$ [Hz]')
plt.ylabel(r'$\left|f(\omega)\right|$')
plt.title(r'$\omega vs \left|f(\omega)\right|$ in positive frequencies')
plt.legend(loc='upper right')
plt.grid(True)

plt.figure(3)
plt.plot(data_pos4[:,0], data_pos4[:,1], label="a=4")
plt.plot(data_pos6[:,0], data_pos6[:,1], label="a=6")
plt.plot(data_pos8[:,0], data_pos8[:,1], label="a=8")
plt.plot(data_pos10[:,0], data_pos10[:,1], label="a=10")
plt.plot(data_neg4[:,0], data_neg4[:,1], label="a=4")
plt.plot(data_neg6[:,0], data_neg6[:,1], label="a=6")
plt.plot(data_neg8[:,0], data_neg8[:,1], label="a=8")
plt.plot(data_neg10[:,0], data_neg10[:,1], label="a=10")
plt.xlabel(r'$\omega$ [Hz]')
plt.ylabel(r'$\left|f(\omega)\right|$')
plt.title(r'$\omega vs \left|f(\omega)\right|$ in both frequencies')
plt.legend(loc='upper left')
plt.grid(True)

plt.show()

