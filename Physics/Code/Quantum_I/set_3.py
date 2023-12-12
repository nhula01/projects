import matplotlib.pyplot as plt
import numpy as np

def plot(name):
    plt.figure(1)
    file_name = '/Users/phihung/NumMethod/first/homework/hw_07/wells1_' + str(name) + '.dat'
    file_name1 = np.loadtxt(file_name)
    #file_name2 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/10_wells2.dat')
    plt.plot(file_name1[:,0], file_name1[:,1], label='well=.'+str(name))
    #plt.plot(file_name2[:,0], file_name2[:,1], label='weell2')
    plt.xlabel(r'distance [a.u.]')
    plt.ylabel(r'Ground Energy [a.u.]')
    plt.title('Energy vs distance between wells')
    plt.legend(loc='lower right')
    plt.grid(True)

    plt.figure(1)
    file_name = '/Users/phihung/NumMethod/first/homework/hw_07/psi0_' + str(name) + '.dat'
    file_name1 = np.loadtxt(file_name)
    file_name = '/Users/phihung/NumMethod/first/homework/hw_07/psi1_' + str(name) + '.dat'
    file_name2 = np.loadtxt(file_name)
    plt.plot(file_name1[:,0], file_name1[:,1], label='ground.' + str(name))
    plt.plot(file_name2[:,0], file_name2[:,1], label='excited.' + str(name))
    plt.xlabel(r'distance [a.u.]')
    plt.ylabel(r'$\psi(x)$')
    plt.title('Energy vs distance between wells')
    plt.legend(loc='lower right')
    plt.grid(True)
    # Create a secondary y-axis on the right
    ax2 = plt.twinx()
    # Set labels and title for the secondary y-axis
    ax2.set_ylabel('Potential Energy [a.u.]', color='blue')

    # Add a legend for the secondary plot
    ax2.legend(loc='upper right')
    plt.show()
    return None

def plotenergy(name,opacity,size):
    plt.figure(2)
    file_name_0 = '/Users/phihung/NumMethod/first/homework/hw_07/converged_energy_E0_' + str(name) + '.dat'
    file_name_1 = '/Users/phihung/NumMethod/first/homework/hw_07/converged_energy_E1_' + str(name) + '.dat'
    data1 = np.loadtxt(file_name_0)
    data2 = np.loadtxt(file_name_1)
    #file_name2 = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_07/10_wells2.dat')
    plt.scatter(data1[0], data1[1],alpha=opacity, s=size,label='ground/'+str(name))
    plt.scatter(data2[0], data2[1],alpha=opacity,s=size, label='excited/'+str(name))
    #plt.plot(file_name2[:,0], file_name2[:,1], label='weell2')
    plt.xlabel(r'repetitions [a.u.]')
    plt.ylabel(r'Energy [a.u.]')
    plt.title('Energy vs repetitions')
    plt.legend(loc='lower right')
    plt.grid(True)
    return None
#plot(2)
plotenergy('normal',1,10)
plotenergy('1imp',.8,80)
plotenergy('1imp2',.8,80)
plotenergy('2end',.6,160)
plotenergy('2mid',.6,160)
plotenergy('3imp',.4,320)
plotenergy('3imp2',.4,320)
plt.show()