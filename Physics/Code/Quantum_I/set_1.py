import matplotlib.pyplot as plt
import numpy as np


energies_50 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/energies50.dat')
energies_10 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/energies10.dat')
energies_5 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/energies5.dat')
energies_2 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/energies2.dat')

plt.figure(1)
plt.plot(energies_50[:,0], energies_50[:,1], 'o', label="num_50") #plot the numerical solution
plt.plot(energies_50[:,0], energies_50[:,2], 'o', label="exact")
plt.plot(energies_10[:,0], energies_10[:,1], 'o', label="num_10") #plot the numerical solution
plt.plot(energies_5[:,0], energies_5[:,1], 'o', label="num_5") #plot the numerical solution
plt.plot(energies_2[:,0], energies_2[:,1], 'o', label="num_2") #plot the numerical solution
plt.xlabel(r'quantum number')
plt.ylabel(r'oscillator energy (a.u.)')
plt.title('Energy')
plt.legend(loc='upper left')
plt.grid(True)

relative_error50 = []
relative_error10 = []
relative_error5 = []
relative_error2 = []
distance = [100, 20, 10, 4]
for i in range(5):
    exact = (.5+i)*.1
    relative_error50.append(np.abs(energies_50[i,3])/exact*100)
    relative_error10.append(np.abs(energies_10[i,3])/exact*100)
    relative_error5.append(np.abs(energies_5[i,3])/exact*100)
    relative_error2.append(np.abs(energies_2[i,3])/exact*100)
plt.figure(3)
plt.plot([distance[0]]*5, relative_error50, 'o', label="100") #plot the numerical solution
plt.plot([distance[1]]*5, relative_error10, 'o', label="20") #plot the numerical solution
plt.plot([distance[2]]*5, relative_error5, 'o', label="10") #plot the numerical solution
plt.plot([distance[3]]*5, relative_error2, 'o', label="4") #plot the numerical solution
plt.xlabel(r'distance b-a (a.u.)')
plt.ylabel(r'relative error')
plt.title('Relative error vs distance')
plt.legend(loc='upper right')
plt.grid(True)
"""
energies_50 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/converged_energy50.dat')
energies_20 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/converged_energy20.dat')
energies_10 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/converged_energy10.dat')
energies_5 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/converged_energy5.dat')
distance = [100, 40, 20, 10]
plt.figure(1)
plt.scatter(distance[0], energies_50[0], label=".1") #plot the numerical solution
plt.scatter(distance[1], energies_20[0], label=".1") #plot the numerical solution
plt.scatter(distance[2], energies_10[0], label=".1") #plot the numerical solution
plt.scatter(distance[3], energies_5[0], label=".1") #plot the numerical solution
plt.xlabel(r'distance [a.u.]')
plt.ylabel(r'$iterations$')
plt.title('Iterations vs Distance')
plt.legend(loc='upper left')
plt.grid(True)

energies_50 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/converged_energy505.dat')
energies_20 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/converged_energy205.dat')
energies_10 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/converged_energy105.dat')
energies_5 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/converged_energy55.dat')
distance = [100, 40, 20, 10]
plt.figure(1)
plt.scatter(distance[0], energies_50[0], label=".5") #plot the numerical solution
plt.scatter(distance[1], energies_20[0], label=".5") #plot the numerical solution
plt.scatter(distance[2], energies_10[0], label=".5") #plot the numerical solution
plt.scatter(distance[3], energies_5[0], label=".5") #plot the numerical solution
plt.legend(loc='upper left')

energies_50 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/ground_state50.dat')
energies_20 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/ground_state20.dat')
energies_10 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/ground_state10.dat')
energies_5 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/ground_state5.dat')
plt.figure(2)
plt.plot(energies_50[:,0], energies_50[:,1], label="num_50_.1") #plot the numerical solution
plt.plot(energies_20[:,0], energies_20[:,1], label="num_20_.1") #plot the numerical solution
plt.plot(energies_10[:,0], energies_10[:,1], label="num_10_.1") #plot the numerical solution
plt.plot(energies_5[:,0], energies_5[:,1], label="num_5_.1") #plot the numerical solution
plt.xlabel(r'distance [a.u.]')
plt.ylabel(r'$\psi(x)$')
plt.title('Eigenstates')
plt.legend(loc='upper left')
plt.grid(True)

energies_50 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/ground_state505.dat')
energies_20 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/ground_state205.dat')
energies_10 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/ground_state105.dat')
energies_5 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/ground_state55.dat')
plt.figure(2)
plt.plot(energies_5[:,0], energies_5[:,1], label="num_5_.5") #plot the numerical solution
plt.legend(loc='upper left')
plt.figure(3)
ground_state = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/ground_state.dat')
plt.plot(ground_state[:,0], ground_state[:,1], label="num_2") #plot the numerical solution
plt.xlabel(r'quantum number')
plt.ylabel(r'oscillator energy (a.u.)$')
plt.title('Energies State')
plt.legend(loc='upper left')
plt.grid(True)


states123_50 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/states123_50.dat')
states123_10 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/states123_10.dat')
states123_5 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/states123_5.dat')
states123_2 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/states123_2.dat')

plt.plot(states123_50[:,0], states123_50[:,1], label="ground state") #plot the numerical solution

plt.plot(states123_50[:,0], states123_50[:,2], label="first state") #plot the numerical solution

plt.plot(states123_50[:,0], states123_50[:,3], label="second state") #plot the numerical solution
plt.xlabel(r'x [a.u.]')
plt.ylabel(r'$\psi(x)$')
plt.title('Eigenfunctions')
plt.legend(loc='lower right')
plt.grid(True)


fig1, axs1 = plt.subplots(1,4)
#plt.figure(1)
axs1[0].plot(states123_50[:,0], states123_50[:,1], label="ground")
axs1[0].plot(states123_50[:,0], states123_50[:,2], label="first")
axs1[0].plot(states123_50[:,0], states123_50[:,3], label="second")
axs1[0].set_xlabel(r'x')
axs1[0].set_ylabel(r'$\psi(x)$')
axs1[0].set_title('Wave Function')
axs1[0].legend(loc='upper right')
axs1[0].grid(True)
axs1[1].plot(states123_10[:,0], states123_10[:,1], label="ground")
axs1[1].plot(states123_10[:,0], states123_10[:,2], label="first")
axs1[1].plot(states123_10[:,0], states123_10[:,3], label="second")
#axs1[1].set_xlabel(r'x')
#axs1[1].set_ylabel(r'$\psi(x)')
#axs1[1].set_title('Wave Function')
axs1[1].legend(loc='upper right')
axs1[1].grid(True)
axs1[2].plot(states123_5[:,0], states123_5[:,1], label="ground")
axs1[2].plot(states123_5[:,0], states123_5[:,2], label="first")
axs1[2].plot(states123_5[:,0], states123_5[:,3], label="second")
#axs1[2].set_xlabel(r'x')
#axs1[2].set_ylabel(r'$\psi(x)')
#axs1[2].set_title('Wave Function')
axs1[2].legend(loc='upper right')
axs1[2].grid(True)
axs1[3].plot(states123_2[:,0], states123_2[:,1], label="ground")
axs1[3].plot(states123_2[:,0], states123_2[:,2], label="first")
axs1[3].plot(states123_2[:,0], states123_2[:,3], label="second")
#axs1[3].set_xlabel(r'x')
#axs1[3].set_ylabel(r'$\psi(x)')
#axs1[3].set_title('Wave Function')
axs1[3].legend(loc='upper right')
axs1[3].grid(True)
"""

plt.show()
