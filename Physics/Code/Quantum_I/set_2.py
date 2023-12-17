import matplotlib.pyplot as plt
import numpy as np
"""
plt.figure(1)
for i in range(1,11):
    file_name = '/Users/phihung/NumMethod/first/homework/hw_07/converged_energy' + str(i) + '.dat'
    label_name = 'd=' + str(i)
    data = np.loadtxt(file_name)
    plt.scatter(i, data[1], label=label_name)
plt.xlabel(r'distance [a.u.]')
plt.ylabel(r'Ground Energy [a.u.]')
plt.title('Energy vs distance between wells')
plt.legend(loc='lower right')
plt.grid(True)

plt.figure(2)
for i in range(1,11):
    file_name = '/Users/phihung/NumMethod/first/homework/hw_07/2_wells' + str(i) + '.dat'
    label_name = 'd=' + str(i)
    opacity = .1 * i
    data = np.loadtxt(file_name)
    plt.plot(data[:,0], data[:,1],alpha=opacity, label=label_name)
plt.xlabel(r'distance [a.u.]')
plt.ylabel(r'Potential [a.u.]')
plt.title('Well')
plt.legend(loc='upper right')
plt.grid(True)


plt.figure(3)
for i in range(1,11):
    file_name = '/Users/phihung/NumMethod/first/homework/hw_07/ground_state_p2' + str(i) + '.dat'
    label_name = 'd=' + str(i)
    opacity = .1 * i
    data = np.loadtxt(file_name)
    plt.plot(data[:,0], data[:,1], alpha=opacity,label=label_name)
plt.xlabel(r'distance [a.u.]')
plt.ylabel(r'$\psi(x)$')
plt.title('Energies State')
plt.legend(loc='lower right')
plt.grid(True)
"""
"""
plt.figure(2)
for i in range(1,11):
    file_name = '/Users/phihung/NumMethod/first/homework/hw_07/energies_p2' + str(i) + '.dat'
    data = np.loadtxt(file_name)
    plt.scatter([i]*20, data[:20,1])
plt.xlabel(r'distance [a.u.]')
plt.ylabel(r'eigenenergies [a.u.]')
plt.title('Four Eigenenergies')
plt.legend(loc='upper right')
plt.grid(True)
"""
"""
depth = [.6,.4,.2,0,-.2,-.4,-.6,-.8,-1,-1.2]
plt.figure(2)
for i in range(1,11):
    file_name = '/Users/phihung/NumMethod/first/homework/hw_07/converged_energy_depth_p2' + str(i) + '.dat'
    data = np.loadtxt(file_name)
    plt.scatter(depth[i-1], data[1])
plt.xlabel(r'distance [a.u.]')
plt.ylabel(r'eigenenergies [a.u.]')
plt.title('Ground Eigenenergies')
plt.legend(loc='upper right')
plt.grid(True)
plt.figure(3)
for i in range(1,11):
    file_name = '/Users/phihung/NumMethod/first/homework/hw_07/ground_state_depth_p2' + str(i) + '.dat'
    data = np.loadtxt(file_name)
    opacity = .1 * i
    plt.plot(data[:,0], data[:,1],alpha=opacity, label='ground'+str(round(.8-i*.2,2)))
plt.xlabel(r'distance [a.u.]')
plt.ylabel(r'$\psi(x)$')
plt.title('Ground State')
plt.legend(loc='upper right')
plt.grid(True)
""""""
depth = [.6,.4,.2,0,-.2,-.4,-.6,-.8,-1,-1.2]
plt.figure(2)
for i in range(1,11):
    file_name = '/Users/phihung/NumMethod/first/homework/hw_07/energies__depth_p2' + str(i) + '.dat'
    data = np.loadtxt(file_name)
    plt.scatter([depth[i-1]]*10, data[:10,1])
plt.xlabel(r'distance [a.u.]')
plt.ylabel(r'eigenenergies [a.u.]')
plt.title('More Eigenenergies')
plt.legend(loc='upper right')
plt.grid(True)"""

for i in range(1,11):
    file_name = '/Users/phihung/NumMethod/first/homework/hw_07/states1,234_depth_p2' + str(i) + '.dat'
    data = np.loadtxt(file_name)
    opacity = .1 * i
    plt.figure(i)
    plt.plot(data[:,0], data[:,1],alpha=opacity, label='ground '+str(round(.8-i*.2,2)))
    plt.plot(data[:,0], data[:,2],alpha=opacity, label='first '+str(round(.8-i*.2,2)))
    plt.plot(data[:,0], data[:,3],alpha=opacity, label='second '+str(round(.8-i*.2,2)))
    plt.plot(data[:,0], data[:,4],alpha=opacity, label='third '+str(round(.8-i*.2,2)))
    plt.xlabel(r'distance [a.u.]')
    plt.ylabel(r'$\psi(x)$')
    plt.title('Three excited state and Ground state')
    plt.legend(loc='upper right')
    plt.grid(True)
"""
plt.figure(3)
for i in range(1,11):
    file_name = '/Users/phihung/NumMethod/first/homework/hw_07/states123_p2' + str(i) + '.dat'
    data = np.loadtxt(file_name)
    opacity = .1 * i
    plt.plot(data[:,0], data[:,1],alpha=opacity, label='ground'+str(i))
plt.xlabel(r'distance [a.u.]')
plt.ylabel(r'First Eneries [a.u.]')
plt.title('Well')
plt.legend(loc='upper right')
plt.grid(True)
"""
plt.show()
