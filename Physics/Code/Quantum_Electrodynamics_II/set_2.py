import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.animation import FuncAnimation

n0 = [(i+1)*1.0e24 for i in range(10)]
omega_minus = [1.3226365114271206,1.3122872429496000,1.3060776818630879,1.2998681207765754,1.2936585596900634,1.2895188522990551,1.2853791449080469,1.2791695838215345,1.2770997301260303,1.2729600227350224]
omega_plus = [1.3723130001192190,1.3826622685967396,1.3909416833787560,1.3971512444652683,1.4033608055517806,1.4095703666382928,1.4137100740293009,1.4178497814203093,1.4219894888113174,1.4261291962023259]
omega_diff = [omega_plus[i] - omega_minus[i] for i in range(10)]
slope = (np.log(omega_diff[1])-np.log(omega_diff[0]))/(np.log(n0[1])-np.log(n0[0]))
slopes = [(np.log(omega_diff[i+1])-np.log(omega_diff[i]))/(np.log(n0[i+1])-np.log(n0[i])) for i in range(9)]
print(slopes)
#data_1 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_10/example3a_spectrum1_0d25.dat')
#plt.scatter(n0,omega_minus,label=r'$\Omega$-')
#plt.scatter(n0,omega_plus,label=r'$\Omega$+')

#plt.scatter(n0,omega_diff,label=r'($\Omega+$) - ($\Omega-$)')
plt.xlabel(r'n0 ($m^{-3}$)')
plt.ylabel('frequency (eV)')
plt.legend()
plt.show()