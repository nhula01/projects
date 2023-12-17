import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.animation import FuncAnimation
from scipy.signal import find_peaks

emptypeaks=[]
data_name = ["00"+str(i*16+380) for i in range(11)] 
def findpeeak(bisect1,bisect2):
    data_name = ["00"+str(i*16+380) for i in range(11)] 
    data = []
    #add all data files
    for name in data_name:
        name = '/Users/phihung/NumMethod/first/homework/hw_10/' + name
        data.append(np.loadtxt(name))
    i=0
    for data_file in data:
        lab = 300 + i*16
        #newdata = data[300:,]
        #plt.plot(data_file[bisect1:bisect2,0], (data_file[bisect1:bisect2,1]), label='distance '+ str(lab))
        peaks, _ = find_peaks((data_file[bisect1:bisect2,1])) #find peaks indices
        print(peaks)
        emptypeaks.append(data_file[peaks[0],0])
        #plt.plot(data_file[peaks[0],0],(data_file[peaks[0],1]),"ro") #plot peak1
        #find the frequency
        i+=1
    #plt.plot(data[1][:,0], np.log(data[1][:,1]), label="30")

    plt.xlabel(r'frequency [eV]')
    plt.ylabel(r'transmission')
    plt.title('')
    plt.legend(loc='lower right')
    plt.grid(True)
    return None
findpeeak(0,-1)
plt.figure(1)
distances = [i*16+380 for i in range(11)] 
omega_minus=[1.3329857799046410,1.3205666577316164,1.2917,1.2586,1.2238,1.1885,1.1536,1.1232,1.0881,1.0578,1.0374]
omega_plus=[1.4281990498978299,1.3930115370742600,1.3686,1.36,1.3584,1.3565,1.3546,1.3545,1.3529,1.3525,1.3524]
line = [1.34699295]*11
plt.plot(distances, emptypeaks,label='cavity (off)')
#distance = np.loadtxt('/Users/phihung/NumMethod/first/homework/hw_10/e00540')
plt.plot(distances, omega_minus,label=r'$\Omega-$ (on)')
plt.plot(distances, omega_plus,label=r'$\Omega+$ (on)')
plt.plot(distances, line,label='avoided cross')
plt.xlabel(r'distance [nm]')
plt.ylabel(r'frequency [eV]')
plt.title('')
plt.legend(loc='lower left')
plt.grid(True)
plt.show()