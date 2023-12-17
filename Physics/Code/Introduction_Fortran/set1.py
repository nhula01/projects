#import necessary libs
import matplotlib.pyplot as plt
import numpy as np

#load data
data = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_01/temp.dat')

#unload data
celsius = data[:,0] # load first column
fahrenheit = data[:,1] # load second colum

#plot the line
plt.plot(celsius,fahrenheit)

#plot the points
plt.scatter(fahrenheit[0],celsius[0],color = 'red',label='avg_temp of j_f_m' )
plt.scatter(fahrenheit[1],celsius[1],color = 'orange',label='avg_temp of a_m_j' )
plt.scatter(fahrenheit[2],celsius[2],color = 'blue',label='avg_temp of j_a_s' )
plt.scatter(fahrenheit[3],celsius[3],color = 'green',label='avg_temp of o_n_d' )

#label
plt.xlabel(r'Fahrenheit [$\degree$F]')
plt.ylabel(r'Celsius [$\degree$C]')

#plot frrom -20C->-4F to 45C->113
y = [-20, 45]
x = [-4, 113]
plt.plot(x,y)

#title of graph
plt.title('C=(F-32)/1.8')

#add legend
plt.legend(loc='upper left')
plt.grid(True)

plt.show()