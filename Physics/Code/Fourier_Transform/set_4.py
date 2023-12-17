import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks


data_name = ["000"+str(i*10+10) for i in range(7)]
data = []
for name in data_name:
    name = '/Users/phihung/NumMethod/first/homework/hw_06/' + name
    data.append(np.loadtxt(name))

def frequency(kb1):
    return (1/(2*np.pi))*np.sqrt(2/kb1)
i=0
v=0

plt.figure(1)
omega_for_plot = [] 
for data_file in data:
    lab = 10 + 20*i
    plt.plot(data_file[:,0], data_file[:,1], label=str(lab))
    peaks, _ = find_peaks(data_file[:,1],prominence=5000)
    plt.plot(data_file[peaks[0],0],data_file[peaks[0],1],"ro")
    plt.plot(data_file[peaks[1],0],data_file[peaks[1],1],"ro")
    #find the frequency
    v = frequency(lab)
    #append the frequencies at peaks
    omega_for_plot.append([data_file[peaks[0],0],v])
    omega_for_plot.append([data_file[peaks[1],0],v])
    i+=1
#plt.plot(data[1][:,0], np.log(data[1][:,1]), label="30")

plt.xlabel(r'$\omega$ [Hz]')
plt.ylabel(r'$\left|x1(\omega)\right|^2$')
plt.title('Fourier Fast Transform')
plt.legend(loc='lower right')
plt.grid(True)

plt.figure(2)
for omega in omega_for_plot:
    plt.plot(omega[1], omega[0],'ro')
plt.ylabel(r'$\omega$ [Hz]')
plt.xlabel(r'$v1$ [Hz]')
plt.title('Frequency of Oscillator 1')
plt.legend(loc='upper right')
plt.grid(True)

plt.show()
