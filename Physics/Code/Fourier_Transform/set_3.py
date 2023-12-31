import matplotlib.pyplot as plt
import numpy as np

data_duff2 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/duffing_beta_random_phase2.dat')
data_duff8 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/duffing_beta_random_phase8.dat')
data_duff14 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/duffing_beta_random_phase14.dat')
"""
#plt.figure(1)
fig1, axs1 = plt.subplots(1,3)
#plt.figure(1)
axs1[0].plot(data_duff2[:,1], data_duff2[:,2], label="2")
axs1[0].set_xlabel(r'$x$ [m]')
axs1[0].set_ylabel(r'$v$ [m/s]')
axs1[0].set_title('Phase Diagram for Duffing')
axs1[0].legend(loc='lower right')
axs1[0].grid(True)
axs1[1].plot(data_duff8[:,1], data_duff8[:,2], label="8")
axs1[1].set_xlabel(r'$x$ [m]')
axs1[1].set_ylabel(r'$v$ [m/s]')
axs1[1].set_title('Phase Diagram for Duffing')
axs1[1].legend(loc='lower right')
axs1[1].grid(True)
axs1[2].plot(data_duff14[:,1], data_duff14[:,2], label="14")
axs1[2].set_xlabel(r'$x$ [m]')
axs1[2].set_ylabel(r'$v$ [m/s]')
axs1[2].set_title('Phase Diagram for Duffing')
axs1[2].legend(loc='lower right')
axs1[2].grid(True)

data_duff2 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/duffing_random_frequency2.dat')
data_duff8 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/duffing_random_frequency8.dat')
data_duff14 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/duffing_random_frequency14.dat')

plt.figure(2)
plt.plot(data_duff2[:,0], np.log(data_duff2[:,1]), label="2")
plt.plot(data_duff8[:,0], np.log(data_duff8[:,1]), label="8")
plt.plot(data_duff14[:,0], np.log(data_duff14[:,1]), label="14")
plt.xlabel(r'$\omega$ [Hz]')
plt.ylabel(r'$\ln(\left|f(\omega)\right|)$')
plt.title('Fourier Fast Transform')
plt.legend(loc='lower right')
plt.grid(True)
"""
"""
data_duff2 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/duffing_beta_random_phase_20_2.dat')
data_duff8 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/duffing_beta_random_phase_20_8.dat')
data_duff14 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/duffing_beta_random_phase_20_14.dat')

#plt.figure(1)
fig2, axs2 = plt.subplots(1,3)
#plt.figure(1)
axs2[0].plot(data_duff2[:,1], data_duff2[:,2], label="2")
axs2[0].set_xlabel(r'$x$ [m]')
axs2[0].set_ylabel(r'$v$ [m/s]')
axs2[0].set_title('Phase Diagram for Duffing')
axs2[0].legend(loc='lower right')
axs2[0].grid(True)
axs2[1].plot(data_duff8[:,1], data_duff8[:,2], label="8")
axs2[1].set_xlabel(r'$x$ [m]')
axs2[1].set_ylabel(r'$v$ [m/s]')
axs2[1].set_title('Phase Diagram for Duffing')
axs2[1].legend(loc='lower right')
axs2[1].grid(True)
axs2[2].plot(data_duff14[:,1], data_duff14[:,2], label="14")
axs2[2].set_xlabel(r'$x$ [m]')
axs2[2].set_ylabel(r'$v$ [m/s]')
axs2[2].set_title('Phase Diagram for Duffing')
axs2[2].legend(loc='lower right')
axs2[2].grid(True)

data_duff2 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/duffing_random_frequency_20_2.dat')
data_duff8 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/duffing_random_frequency_20_8.dat')
data_duff14 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/duffing_random_frequency_20_14.dat')

plt.figure(3)
plt.plot(data_duff2[:,0], np.log(data_duff2[:,1]), label="2")
plt.plot(data_duff8[:,0], np.log(data_duff8[:,1]), label="8")
plt.plot(data_duff14[:,0], np.log(data_duff14[:,1]), label="14")
plt.xlabel(r'$\omega$ [Hz]')
plt.ylabel(r'$\ln(\left|f(\omega)\right|)$')
plt.title('Fourier Fast Transform')
plt.legend(loc='lower right')
plt.grid(True)
"""
data_duff2 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/duffing_beta_random_phase_1_2.dat')
data_duff8 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/duffing_beta_random_phase_1_8.dat')
data_duff14 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/duffing_beta_random_phase_1_14.dat')

#plt.figure(1)
fig2, axs2 = plt.subplots(1,3)
#plt.figure(1)
axs2[0].plot(data_duff2[:,1], data_duff2[:,2], label="2")
axs2[0].set_xlabel(r'$x$ [m]')
axs2[0].set_ylabel(r'$v$ [m/s]')
axs2[0].set_title('Phase Diagram for Duffing')
axs2[0].legend(loc='lower right')
axs2[0].grid(True)
axs2[1].plot(data_duff8[:,1], data_duff8[:,2], label="8")
axs2[1].set_xlabel(r'$x$ [m]')
axs2[1].set_ylabel(r'$v$ [m/s]')
axs2[1].set_title('Phase Diagram for Duffing')
axs2[1].legend(loc='lower right')
axs2[1].grid(True)
axs2[2].plot(data_duff14[:,1], data_duff14[:,2], label="14")
axs2[2].set_xlabel(r'$x$ [m]')
axs2[2].set_ylabel(r'$v$ [m/s]')
axs2[2].set_title('Phase Diagram for Duffing')
axs2[2].legend(loc='lower right')
axs2[2].grid(True)

data_duff2 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/duffing_random_frequency_1_2.dat')
data_duff8 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/duffing_random_frequency_1_8.dat')
data_duff14 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_06/duffing_random_frequency_1_14.dat')

plt.figure(3)
plt.plot(data_duff2[:,0], np.log(data_duff2[:,1]), label="2")
plt.plot(data_duff8[:,0], np.log(data_duff8[:,1]), label="8")
plt.plot(data_duff14[:,0], np.log(data_duff14[:,1]), label="14")
plt.xlabel(r'$\omega$ [Hz]')
plt.ylabel(r'$\ln(\left|f(\omega)\right|)$')
plt.title('Fourier Fast Transform')
plt.legend(loc='lower right')
plt.grid(True)
plt.show()

