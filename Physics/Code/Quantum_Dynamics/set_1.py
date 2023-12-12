import numpy as np
import matplotlib.pyplot as plt

# Example Z values (replace with your own data)
def plot1(name): #positive
    file_name = '/Users/phihung/NumMethod/first/homework/hw_08/psi_xt-state' + name + '.dat'
    Z = np.loadtxt(file_name)
    xshape, tshape= Z.shape
    a = -50
    b = 50
    tp = 50
    dx = (b-a)/xshape
    x = [(a+dx*i) for i in range(xshape)]
    dt = tp/tshape
    t = [(0+dt*i) for i in range(tshape)]
    X, T = np.meshgrid(x, t) 
    Z = np.transpose(Z)
    # Create the contour plot
    contours = plt.contour(X, T, Z, levels=80)  # You can adjust the number of levels as needed
    plt.xlabel('x [a.u.]')
    plt.ylabel('t [s]')
    plt.title('Propagation of Wave Function')
    plt.colorbar(label='Probability')  # Add a color bar to indicate Z-values
    return None

def plot3(name): #plot the well
    file_name = '/Users/phihung/NumMethod/first/homework/hw_08/well' + name + '.dat'
    Z = np.loadtxt(file_name)
    plt.plot(Z[:,0],Z[:,1])
    return None
 
def plot4(name):
    file_name = '/Users/phihung/NumMethod/first/homework/hw_08/x_vs_t-state' + str(name) +'.dat'
    Z = np.loadtxt(file_name) #'/Users/phihung/NumMethod/first/homework/hw_08/x_vs_t-state1'
    plt.plot(Z[:,0],Z[:,1])
    plt.xlabel('t [a.u.]')
    plt.ylabel('<x> [a.u.]')
    plt.title('Propagation of Wave Function')
    return None
#plot1('1oscillate2')
#plot1('2oscillate2')
plot4('1oscillate2')
#plot4('2oscillate5')
plt.show()