import matplotlib.pyplot as plt
import numpy as np

data_set7 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/tau_e_pop_7")

data_set5 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/tau_e_pop_5")


plt.figure(1)
plt.plot(data_set7[:,0], data_set7[:,1], label="V=.7V/nm")
plt.plot(data_set5[:,0], data_set5[:,1], label="V=.5V/nm")
plt.xlabel('tau [s]')
plt.ylabel('probability of excited state')
plt.title('tau vs probability')
plt.legend(loc='lower right')
plt.grid(True)

data_set5=np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_05/2level_system_5.dat')
data_set7=np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_05/2level_system_7.dat')

plt.figure(2)
plt.plot(data_set7[:,0], data_set7[:,1], label="V=.7V/nm")
plt.plot(data_set5[:,0], data_set5[:,1], label="V=.5V/nm")
plt.xlabel('time [s]')
plt.ylabel('probability of excited state')
plt.title('time vs probability')
plt.legend(loc='lower right')
plt.grid(True)
plt.show()

