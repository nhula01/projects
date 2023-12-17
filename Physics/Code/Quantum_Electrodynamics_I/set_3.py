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
metal1= [-.4e-6]*10
metal2= [.4e-6]*10
metal3= [0]*10
def coordinate1():
    E_part = '/Users/phihung/NumMethod/first/homework/hw_09/electric_t_3.dat'
    H_part = '/Users/phihung/NumMethod/first/homework/hw_09/magnetic_t_3.dat'
    detection_evol = '/Users/phihung/NumMethod/first/homework/hw_09/detection3.dat'
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
    line8, = ax1.plot(metal3, detection,label=r'Metal')
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
    ani.save('/Users/phihung/NumMethod/first/homework/hw_09/problem3.gif', writer='pillow')
    ax1.set_ylim(-.2, .2)
    ani.save('/Users/phihung/NumMethod/first/homework/hw_09/problem3_zoomedin.gif', writer='pillow')
    # Show the plot
    plt.show()
    return None
#coordinate1()

#data = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_09/example_5-880nm-20nm_middle.dat')
#plt.plot(data[:300,0],np.log(data[:300,1]))

emptypeaks=[]
def findpeeak(bisect1,bisect2):

    data_name = ["000"+str(i*10+10) for i in range(7)] 
    data = []
    #add all data files
    for name in data_name:
        name = '/Users/phihung/NumMethod/first/homework/hw_09/' + name
        data.append(np.loadtxt(name))

    i=0
    for data_file in data:
        lab = 10 + 10*i
        #newdata = data[300:,]
        plt.plot(data_file[bisect1:bisect2,0], np.log(data_file[bisect1:bisect2,1]), label='thickness '+ str(lab))
        peaks, _ = find_peaks(np.log(data_file[bisect1:bisect2,1])) #find peaks indices
        print(peaks)
        emptypeaks.append(data_file[peaks[0]+bisect1,0])
        plt.plot(data_file[peaks[0]+bisect1,0],np.log(data_file[peaks[0]+bisect1,1]),"ro") #plot peak1
        try:
            plt.plot(data_file[peaks[1]+bisect1,0],np.log(data_file[peaks[1]+bisect1,1]),"ro") # plot peak2
            emptypeaks.append(data_file[peaks[1]+bisect1,0])
        except IndexError:
            print('opps')
        #find the frequency
        i+=1
    #plt.plot(data[1][:,0], np.log(data[1][:,1]), label="30")

    plt.xlabel(r'frequency [eV]')
    plt.ylabel(r'transmission')
    plt.title('')
    plt.legend(loc='lower right')
    plt.grid(True)
    return None
#findpeeak(0,300)
#print(emptypeaks)
peaks_1 = [1.1384195325272557, 1.4406181720708544,
            1.2543313394754854, 1.4406181720708544, 
            1.3247063651226247, 1.4447578794618625, 
           1.366103439032707, 1.4488975868528708,
           1.405, 1.4571770016348873, 
           1.4321,1.4654564164169035,
           1.45, 1.47373583119892] #last 3 peaks can not find first peak
#findpeeak(300,-1)
print(emptypeaks)
peaks_2 = [2.4051699941757656, 2.881236344141709, 
           2.541780338079036, 2.881236344141709,
             2.645273022854241, 2.889515758923725, 
           2.719787755892389, 2.8977951737057417, 
           2.7943024889305366, 2.9143540032697746,
             2.848118685013643, 2.930912832833807, 
           2.8975, 2.951611369788848] #add the last first peeak

distance = [i*10+10 for i in range(7)]
peak1diff = [peaks_1[i*2+1]-peaks_1[i*2] for i in range(int(len(peaks_1)/2))]
peak2diff = [peaks_2[i*2+1]-peaks_2[i*2] for i in range(int(len(peaks_2)/2))]
plt.plot(distance, peak1diff, label='natural mode')

plt.plot(distance, peak2diff,label='first mode')
plt.xlabel(r'thickness [nm]')
plt.ylabel(r'splitting [eV]')
plt.title('Energy Splitting')
plt.legend(loc='upper right')
plt.grid(True)
plt.show()
