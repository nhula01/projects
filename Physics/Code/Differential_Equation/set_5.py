import matplotlib.pyplot as plt
import numpy as np

data_set0 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/interacting_emitters-E0_scan_0.dat")
data_set1 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/interacting_emitters-E0_scan_1.dat")
data_set2 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/interacting_emitters-E0_scan_2.dat")
data_set3 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/interacting_emitters-E0_scan_3.dat")
data_set4 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/interacting_emitters-E0_scan_4.dat")
data_set5 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/interacting_emitters-E0_scan_5.dat")
data_set6 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/interacting_emitters-E0_scan_.1.dat")



plt.figure(1)
plt.plot(data_set0[:,0], data_set0[:,3], label="delta=0")
plt.plot(data_set1[:,0], data_set1[:,3], label="delta=1")
plt.plot(data_set2[:,0], data_set2[:,3], label="delta=2")
plt.plot(data_set3[:,0], data_set3[:,3], label="delta=3")
plt.plot(data_set4[:,0], data_set4[:,3], label="delta=4")
plt.plot(data_set5[:,0], data_set5[:,3], label="delta=5")

plt.plot(data_set6[:,0], data_set6[:,3], label="delta=.1")
plt.xlabel('electric field [V/m]')
plt.ylabel('probability of excited state')
plt.title('Electric field vs probability of excited state')
plt.legend(loc='lower right')
plt.grid(True)


data_set0 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/delta_E01.dat")
plt.figure(2)
plt.plot(data_set0[:,0], data_set0[:,1],color='y', label="0")
plt.xlabel('Coupling constant delta [rad/sec]')
plt.ylabel('Electric Field [V/m]')
plt.title('Coupling Constant vs Electric field')
plt.legend(loc='lower right')
plt.grid(True)
plt.show()

