import matplotlib.pyplot as plt
import numpy as np

data_set1 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/p2_dx1.dat")
data_set2 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/p2_dx2.dat")
data_set3 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/p2_dx3.dat")
data_set4 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/p2_dx4.dat")
"""
data_euler_set1 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/p2_euler_dx1.dat")
data_euler_set2 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/p2_euler_dx2.dat")
data_euler_set3 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/p2_euler_dx3.dat")
data_euler_set4 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/p2_euler_dx4.dat")
"""
plt.figure(1)
plt.plot(data_set1[:,0], data_set1[:,1],'o', color='y', label="1d-2")
plt.plot(data_set2[:,0], data_set2[:,1],'o', color='r', label="1d-1")
plt.plot(data_set3[:,0], data_set3[:,1],'o', color='b', label=".15")
plt.plot(data_set4[:,0], data_set4[:,1],'o', color='g', label=".3")
plt.xlabel('x')
plt.ylabel('y')
plt.title('x vs y')
plt.legend(loc='lower right')
plt.grid(True)

# specific value at x=.45
data_1 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/00001")
data_2 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/00002")
data_3 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/00003")
data_4 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/00004")
data_5 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/00005")
data_6 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/00006")
data_7 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/00007")
data_8 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/00008")

plt.figure(2)
#plt.scatter(data_1[0], data_1[1], label='1')
import math
y = (1/2)*np.sqrt(np.pi)*np.exp(.45**2)*math.erf(.45)
#plt.scatter(data_2[0], data_2[1], label='2')
#plt.scatter(data_3[0], data_3[1], label='3')
plt.scatter(.45,y,s=500, alpha=.1, label = "exact")
plt.scatter(data_1[0], data_1[1],s=100, alpha=.5, label='euler 1d-2')
plt.scatter(data_2[0], data_2[1],s=200, alpha=.4, label='euler 1d-1')
plt.scatter(data_3[0], data_3[1],s=300, alpha=.3, label='euler .15')
plt.scatter(data_4[0], data_4[1],s=400, alpha=.2, label='euler .3')
plt.scatter(data_5[0], data_5[1],s=100, alpha=.4, label='rk4 1d-2')
plt.scatter(data_6[0], data_6[1],s=200, alpha=.3, label='rk4 1d-1')
plt.scatter(data_7[0], data_7[1],s=300, alpha=.2, label='rk4 .15')
plt.scatter(data_8[0], data_8[1],s=200, alpha=.2, label='rk4 .3')
plt.xlabel('x')
plt.ylabel('y [arb.unit]')
plt.title('x vs y')
plt.legend(loc='upper left')
plt.grid(True)

plt.show()

