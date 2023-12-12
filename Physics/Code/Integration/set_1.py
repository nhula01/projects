import matplotlib.pyplot as plt
import numpy as np

pi_3 = np.loadtxt("/Users/phihung/Documents/PHY/first/homework/hw_04/integral_record_set1a.dat")
pi_4 = np.loadtxt("/Users/phihung/Documents/PHY/first/homework/hw_04/integral_record_set1b.dat")
pi_40 = np.loadtxt("/Users/phihung/Documents/PHY/first/homework/hw_04/integral_record_set1c.dat")
specific_coor = np.loadtxt("/Users/phihung/Documents/PHY/first/homework/hw_04/specific_comparison_1.dat")
theta_T = np.loadtxt("/Users/phihung/Documents/PHY/first/homework/hw_04/theta_vs_T.dat")

a_convergence = np.loadtxt("/Users/phihung/Documents/PHY/first/homework/hw_04/a_convergence.dat")
b_convergence = np.loadtxt("/Users/phihung/Documents/PHY/first/homework/hw_04/b_convergence.dat")
c_convergence = np.loadtxt("/Users/phihung/Documents/PHY/first/homework/hw_04/c_convergence.dat")

a_x = a_convergence[:,0]
a_y= a_convergence[:,1]

b_x = b_convergence[:,0]
b_y= b_convergence[:,1]

c_x = c_convergence[:,0]
c_y= c_convergence[:,1]
plt.plot(a_x,a_y)
plt.plot(b_x,b_y)
plt.plot(c_x,c_y)

plt.xlabel('integration step [arb. unit]')
plt.ylabel('K [arb. unit]]') 

"""
theta_specific = specific_coor[:,0]
T_specific = specific_coor[:,1]
plt.scatter(theta_specific,T_specific,s=500,label='specific coor', alpha=1,color='b')

plt.xlabel('z [unitless]')
plt.ylabel('period [s]') 
"""

"""
# plot pi_3
iter = pi_3[:,0]
K_value = pi_3[:,1]
K_actual = 1.68575
plt.plot(iter, K_value,label='pi/3')
plt.scatter(iter[-1], K_actual, label='actual at pi/3', color='m')

# plot pi_4
iter = pi_4[:,0]
K_value = pi_4[:,1]
K_actual = 1.63359
plt.plot(iter, K_value,label='pi/4')
plt.scatter(iter[-1], K_actual, label='actual at pi/4', color='k')

# plot pi_40
iter = pi_40[:,0]
K_value = pi_40[:,1]
K_actual = 1.5714
plt.plot(iter, K_value,label='pi/40')
plt.scatter(iter[-1], K_actual, label='actual at pi/40', color='b')

plt.xlabel('x [arb. unit]')
plt.ylabel('K [arb. unit]]') 
"""

"""

theta_specific = specific_coor[:,0]
T_specific = specific_coor[:,1]
plt.scatter(theta_specific,T_specific,s=500,label='specific coor', alpha=1,color='b')
theta = theta_T[:,0]
T = theta_T[:,1]
plt.plot(theta, T,'o',label='general', color='r')

T = 2*np.pi*np.sqrt(.25/9.8)
x = np.linspace(0,3.145,1000)
K = []
n=0
for i in x:
    K.append(T)
plt.plot(x,K,label='Ground state')

plt.xlabel('theta [radians]')
plt.ylabel('T [period]]')
plt.title('T vs theta')"""

plt.legend(loc='upper left')
plt.grid(True)

plt.show()