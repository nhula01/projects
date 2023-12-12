import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.animation import FuncAnimation

dz=1.0e-9
Nz=3000 
z0=-dz*Nz*0.5
zM=dz*Nz*0.5-dz 
tp=35

detection = [-1+i*.2 for i in range(10)]
dposi= [zM-1500*dz]*10
sposi= [z0+100*dz]*10
def coordinate():
    E_part = '/Users/phihung/NumMethod/first/homework/hw_09/electric_t.dat'
    H_part = '/Users/phihung/NumMethod/first/homework/hw_09/magnetic_t.dat'
    detection_evol = '/Users/phihung/NumMethod/first/L14-Electrodynamics/detection.dat'
    E = np.loadtxt(E_part)
    H = np.loadtxt(H_part)
    detect = np.loadtxt(detection_evol)
    zshape, tshape= E.shape
    tp = 35
    dz = (zM-z0)/zshape
    z = [(z0+dz*i) for i in range(zshape)]
    zH = [(z0+dz*i+dz/2) for i in range(zshape-1)]
    dt = tp/tshape
    t = [(0+dt*i) for i in range(tshape)]
    num_points = 500
    #create the metal
    #define metal
    metal=[z[-1]]*10
    # Initialize the line object
    fig, (ax1, ax2) = plt.subplots(2, 1)
    line1, = ax1.plot(z, E[:,0],label=r'E')
    line2, = ax1.plot(zH, H[:,0],label=r'H')
    line3, = ax1.plot(dposi, detection,label=r'Detector')
    line5, = ax1.plot(sposi, detection,label=r'Source')
    line6, = ax1.plot(metal, detection,label=r'Metal')
    #line4, = ax2.plot(detect[:,0], detect[:,1],label=r'Detection')
    #ax1.legend(fontsize='small')
    ax1.set_ylim(-3, 3)
    ax2.legend()
    ax2.set_ylim(-1, 1)
    title1 = ax1.set_title('')
    title2 = ax2.set_title('')
    # Function to update the plot for each frame of the animation
    def update(frame):
        frame = frame + 1
        t = frame*dt
        line1.set_ydata(E[:,frame])
        line2.set_ydata(H[:,frame])
        #line4.set_data(detect[0:frame,0], detect[0:frame,1]) 
        title1.set_text('Time = {:.2f}'.format(frame*dt))
        #title2.set_text('Time = {:.2f}'.format(frame*dt))
        return line1, line2,line3,line5,line6, title1#, title2,line4
    # Create the animation
    num_frames = tshape-1
    ani = FuncAnimation(fig, update, frames=num_frames, interval=80)
    ax1.set_xlabel('z')
    ax1.set_ylabel('E/H')
    ax2.set_xlabel('t')
    ax2.set_ylabel('E')
    #ani.save('/Users/phihung/NumMethod/first/L14-Electrodynamics/animated_plot_9_wavepacket.gif', writer='pillow')
    #ax1.set_ylim(-.2, .2)
    #ani.save('/Users/phihung/NumMethod/first/L14-Electrodynamics/animated_plot_10_wavepacket.gif', writer='pillow')
    # Show the plot
    plt.show()
    return None

data = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_09/speed_distance.dat')
plt.scatter(data[:,0],data[:,1],label='estimate')

def actual(x):
    return 299792458
x = np.linspace(-50,50,20)
y = [actual(x) for i in range(20)]
plt.plot(x,y,label='real')
plt.xlabel('distance from 500nm (nm)')
plt.ylabel('speed of light (m/s)')
plt.legend()
plt.grid()
plt.show()