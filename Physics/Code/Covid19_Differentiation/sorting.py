import matplotlib.pyplot as plt
import numpy as np

file_path = '/Users/phihung/Documents/PHY/first/homework/hw_03/df_dx_central.dat'
file_path2 = '/Users/phihung/Documents/PHY/first/homework/hw_03/df_dx_central_ordered.dat'
data_central = np.loadtxt(file_path)
x2 = data_central[:,0]
df_dx2 = data_central[:,1]
#sort central
x2,df_dx2 = zip(*sorted(zip(x2,df_dx2)))
#reweite the file
with open(file_path2,'w') as file:
    for i in range(len(x2)):
        file.write("{}\t{}\n".format(x2[i], df_dx2[i]))