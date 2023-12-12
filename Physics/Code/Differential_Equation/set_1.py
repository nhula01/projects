import matplotlib.pyplot as plt
import numpy as np

def exact(x):
    return np.log(x+1)
x_axis = np.linspace(0,1,100)
y_axis = exact(x_axis)
# Euler
euler1 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/euler_1.dat")
euler2 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/euler_2.dat")
euler3 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/euler_3.dat")

rk21 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/runge_kutta_2nd_1.dat")
rk22 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/runge_kutta_2nd_2.dat")
rk23 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/runge_kutta_2nd_3.dat")

rk41 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/runge_kutta_4th_1.dat")
rk42 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/runge_kutta_4th_2.dat")
rk43 = np.loadtxt("/Users/phihung/NumMethod/first/homework/hw_05/runge_kutta_4th_3.dat")

fig1, axs1 = plt.subplots(1,3)
#plt.figure(1)
axs1[0].plot(euler1[:,0], euler1[:,1],'o', color='y', label="euler")
axs1[0].plot(rk21[:,0], rk21[:,1],'o', color='r', label="rk2")
axs1[0].plot(rk41[:,0], rk41[:,1],'o', color='b',label="rk4")
axs1[0].plot(x_axis, y_axis, label="exact")
axs1[0].set_xlabel('x')
axs1[0].set_ylabel('y')
axs1[0].set_title('1d-2')
axs1[0].legend(loc='lower right')
axs1[0].grid(True)

axs1[1].plot(euler2[:,0], euler2[:,1],'o', color='y', label="euler")
axs1[1].plot(rk22[:,0], rk22[:,1],'o', color='r', label="rk2")
axs1[1].plot(rk42[:,0], rk42[:,1],'o', color='b',label="rk4")
axs1[1].plot(x_axis, y_axis, label="exact")
axs1[1].set_xlabel('x')
axs1[1].set_ylabel('y')
axs1[1].set_title('5d-2')
axs1[1].legend(loc='lower right')
axs1[1].grid(True)

axs1[2].plot(euler3[:,0], euler3[:,1],'o', color='y', label="euler")
axs1[2].plot(rk23[:,0], rk23[:,1],'o', color='r', label="rk2")
axs1[2].plot(rk43[:,0], rk43[:,1],'o', color='b',label="rk4")
axs1[2].plot(x_axis, y_axis, label="exact")
axs1[2].set_xlabel('x')
axs1[2].set_ylabel('y')
axs1[2].set_title('.1')
axs1[2].legend(loc='lower right')
axs1[2].grid(True)

plt.show()