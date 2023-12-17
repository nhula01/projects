import numpy as np
import matplotlib.pyplot as plt

# Example Z values (replace with your own data)
def plot1(): #positive
    file_name = '/Users/phihung/NumMethod/first/homework/hw_08/Nguyen_HW8/Evspop3.dat'
    Z = np.loadtxt(file_name)
    plt.plot(Z[:,0],Z[:,1], label='ground')
    plt.plot(Z[:,0],Z[:,2],label='excited')
    plt.xlabel('10*E [V/nm]')
    plt.ylabel('Electronic Population')
    plt.grid()
    plt.legend()
    plt.title('Electric Field vs Probability')
    return None
 
plot1()
plt.show()