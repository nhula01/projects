import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.animation import FuncAnimation
from scipy.signal import find_peaks

dz=1.0e-9
Nz=3000 
z0=-dz*Nz*0.5
zM=dz*Nz*0.5-dz 
tp=35

detection = [-1+i*.2 for i in range(10)]
dposi= [zM-1500*dz]*10
sposi= [z0+100*dz]*10
def coordinate():
    E_part = '/Users/phihung/NumMethod/first/L14-Electrodynamics/electric_t_5.dat'
    H_part = '/Users/phihung/NumMethod/first/L14-Electrodynamics/magnetic_t_5.dat'
    #detection_evol = '/Users/phihung/NumMethod/first/L14-Electrodynamics/detection4a.dat'
    E = np.loadtxt(E_part)
    H = np.loadtxt(H_part)
    #detect = np.loadtxt(detection_evol)
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

dz=1.0e-9
Nz=3000 
z0=-dz*Nz*0.5
zM=dz*Nz*0.5-dz 
tp=1000
detection = [-1+i*.2 for i in range(10)]
dposi= [zM-69*dz]*10 #detector position
sposi= [z0+100*dz]*10
metal1= [-300*dz]*10
metal2= [300*dz]*10
metal3= [-300*dz + 600*dz*i/10 for i in range(10)]
def coordinate1():
    E_part = '/Users/phihung/NumMethod/first/homework/hw_09/electric_t_2_glass.dat'
    H_part = '/Users/phihung/NumMethod/first/homework/hw_09/magnetic_t_2_glass.dat'
    detection_evol = '/Users/phihung/NumMethod/first/homework/hw_09/detection2_glass.dat'
    E = np.loadtxt(E_part)
    H = np.loadtxt(H_part)
    detect = np.loadtxt(detection_evol)
    zshape, tshape= E.shape
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
    fig, (ax1, ax2) = plt.subplots(2,1)
    line1, = ax1.plot(z, E[:,0],label=r'E')
    line2, = ax1.plot(zH, H[:,0],label=r'H')
    line3, = ax1.plot(dposi, detection,label=r'Detector')
    line5, = ax1.plot(sposi, detection,label=r'Source')
    line6, = ax1.plot(metal1, detection,label=r'Metal')
    line7, = ax1.plot(metal2, detection,label=r'Metal')
    line8, = ax1.plot(metal3, [0]*10,label=r'Dielectric')
    line4, = ax2.plot(detect[:,0], detect[:,1],label=r'Detection')
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
        line4.set_data(detect[0:frame,0], detect[0:frame,1]) 
        title1.set_text('Time = {:.2f}'.format(frame*dt))
        #title2.set_text('Time = {:.2f}'.format(frame*dt))
        return line1, line2,line3,line5,line6,line7, title1,line4,line8
    # Create the animation
    num_frames = tshape-1
    ani = FuncAnimation(fig, update, frames=200, interval=80)
    ax1.set_xlabel('z')
    ax1.set_ylabel('E/H')
    ax2.set_xlabel('t')
    ax2.set_ylabel('E')
    ax2.set_ylim(-.2, .2)
    ax2.set_xlim(0, 35)
    ani.save('/Users/phihung/NumMethod/first/homework/hw_09/problem2_1_5.gif', writer='pillow')
    ax1.set_ylim(-.2, .2)
    ani.save('/Users/phihung/NumMethod/first/homework/hw_09/problem2_zoomedin_1_5.gif', writer='pillow')
    plt.show()
    return None
#coordinate1()

#data = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_09/010002_0')
#plt.plot(data[:,0],np.log(data[:,1]))
emptypeaks=[]
def findpeeak(bisect1,bisect2,name,number):
    #emptypeaks=[]
    name1=name
    data_name = ["00"+str(i*100+300)+name for i in range(7)] +['01000' + name]
    data = []
    #add all data files
    for name in data_name:
        name = '/Users/phihung/NumMethod/first/homework/hw_09/' + name
        data.append(np.loadtxt(name))
    i=0
    #plt.figure(number)
    for data_file in data:
        lab = 300 + 100*i
        #newdata = data[300:,]
        #plt.plot(data_file[bisect1:bisect2,0], np.log(data_file[bisect1:bisect2,1]), label='distance '+ str(lab) +' '+name1)
        peaks, _ = find_peaks(np.log(data_file[bisect1:bisect2,1]),prominence=.5) #find peaks indices
        emptypeaks.append(data_file[peaks[0]+bisect1,0])
        #plt.plot(data_file[peaks[0]+bisect1,0],np.log(data_file[peaks[0]+bisect1,1]),"ro") #plot peak1
        try:
         #   plt.plot(data_file[peaks[1]+bisect1,0],np.log(data_file[peaks[1]+bisect1,1]),"ro") # plot peak2
            emptypeaks.append(data_file[peaks[1]+bisect1,0])
        except IndexError:
            print('opps')
        #find the frequency
        i+=1
    #plt.plot(data[1][:,0], np.log(data[1][:,1]), label="30")

    #plt.xlabel(r'frequency [eV]')
    #plt.ylabel(r'transmission')
    #plt.title('')
    #plt.legend(loc='lower right')
    #plt.grid(True)
    return None

#find analytical solution
refractive_index = [1.5,2,2.5,3,3.5]
distance = [i*100+300 for i in range(8)] 
analytical =[]
for refractive in refractive_index: #for each refractive index
    for d in distance: #for each distance
        analytical.append(1239.8/(d*2*refractive))

findpeeak(0,-1,'1_5',1)
dialec_constant_1_5 = emptypeaks
#plt.plot(distance, analytical[:8], label='analytical at 1.5')
#plt.plot(distance, [dialec_constant_1_5[i*2] for i in range(8)], label='numercal at 1.5')


emptypeaks=[]
findpeeak(0,-1,'2_0',2)
dialec_constant_2_0 = emptypeaks 
#plt.plot(distance, analytical[8:16], label='analytical at 2.0')
#plt.plot(distance, [dialec_constant_2_0[i*2] for i in range(8)], label='numercal at 2.0')

emptypeaks=[]
findpeeak(0,-1,'2_5',3)
dialec_constant_2_5 = emptypeaks 
#plt.plot(distance, analytical[16:24], label='analytical at 2.5')
#plt.plot(distance, [dialec_constant_2_5[i*2] for i in range(8)], label='numercal at 2.5')

emptypeaks=[]
findpeeak(0,-1,'3_0',4)
dialec_constant_3_0 = emptypeaks 
#plt.plot(distance, analytical[24:32], label='analytical at 3.0')
#plt.plot(distance, [dialec_constant_3_0[i*2] for i in range(8)], label='numercal at 3.0')

emptypeaks=[]
findpeeak(0,-1,'3_5',5)
dialec_constant_3_5 = emptypeaks 
#plt.plot(distance, analytical[32:40], label='analytical at 3.5')
#plt.plot(distance, [dialec_constant_3_5[i*2] for i in range(8)], label='numercal at 3.5')

#plot frequency vs distance for first peak fixed ri=1_5
"""
plt.figure(6)
plt.plot(distance,[dialec_constant_1_5[i*2] for i in range(8)],label='natural at 1.5') #
plt.plot(distance,[dialec_constant_1_5[i*2+1] for i in range(8)],label='first mode at 1.5') 

plt.plot(distance,[dialec_constant_2_0[i*2] for i in range(8)],label='natural at 2.0') #
plt.plot(distance,[dialec_constant_2_0[i*2+1] for i in range(8)],label='first mode at 2.0') 

plt.plot(distance,[dialec_constant_2_5[i*2] for i in range(8)],label='natural at 2.5') #
plt.plot(distance,[dialec_constant_2_5[i*2+1] for i in range(8)],label='first mode at 2.5') 

plt.plot(distance,[dialec_constant_3_0[i*2] for i in range(8)],label='natural at 3.0') #
plt.plot(distance,[dialec_constant_3_0[i*2+1] for i in range(8)],label='first mode at 3.0') 

plt.plot(distance,[dialec_constant_3_5[i*2] for i in range(8)],label='natural at 3.5') #
plt.plot(distance,[dialec_constant_3_5[i*2+1] for i in range(8)],label='first mode at 3.5') """
#similarly for other ri

#8 data sets, each gives 2 modes ->16, there are 5 reflective indices
#plt.figure(7)
dialec_constant = dialec_constant_1_5 + dialec_constant_2_0 + dialec_constant_2_5 + dialec_constant_3_0 +dialec_constant_3_5
#plot frequency vs ri with fixeed distance=300
plt.plot(refractive_index,[dialec_constant[i*16] for i in range(5)], label='natural frequency at 300nm') #first peak #do it forr other oarts
plt.plot(refractive_index,[dialec_constant[i*16+1] for i in range(5)],label='first mode at 300nm') #second 


plt.plot(refractive_index,[dialec_constant[i*16+4] for i in range(5)], label='natural frequency at 500nm') #first peak #do it forr other oarts
plt.plot(refractive_index,[dialec_constant[i*16+5] for i in range(5)],label='first mode at 500nm') #second 


plt.plot(refractive_index,[dialec_constant[i*16+8] for i in range(5)], label='natural frequency at 700nm') #first peak #do it forr other oarts
plt.plot(refractive_index,[dialec_constant[i*16+9] for i in range(5)],label='first mode at 700nm') #second 

plt.plot(refractive_index,[dialec_constant[i*16+12] for i in range(5)], label='natural frequency at 900nm') #first peak #do it forr other oarts
plt.plot(refractive_index,[dialec_constant[i*16+13] for i in range(5)],label='first mode at 900nm') #second 

plt.xlabel(r'refractive index')
plt.ylabel(r'frequency [eV]')
plt.title('')
plt.legend(loc='upper right')

plt.grid(True)
plt.show()
