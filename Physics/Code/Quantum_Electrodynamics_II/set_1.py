import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.animation import FuncAnimation
"""
data_1 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_10/e0_rho.dat')
plt.plot(data_1[:,0],data_1[:,1],label='gamma_dp=0')
plt.plot(data_1[:,0],data_1[:,2],label='gamma_dp=2.0*pi/20.0D-15')
plt.plot(data_1[:,0],data_1[:,3],label='gamma_dp=2.0*pi/200.0D-15')
plt.plot(data_1[:,0],data_1[:,4],label='gamma_dp=2.0*pi/400.0D-15')
plt.xlabel('E0 (V/nm)')
plt.ylabel('abs(rho)')
plt.legend()"""
data_1 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_10/tau_rho.dat')
plt.plot(data_1[:,0],data_1[:,1],label='gamma_dp=0')
plt.plot(data_1[:,0],data_1[:,2],label='gamma_dp=2.0*pi/20.0D-15')
plt.plot(data_1[:,0],data_1[:,3],label='gamma_dp=2.0*pi/200.0D-15')
plt.plot(data_1[:,0],data_1[:,4],label='gamma_dp=2.0*pi/400.0D-15')
plt.xlabel('tau (fs)')
plt.ylabel('abs(rho)')
plt.legend()
plt.show()