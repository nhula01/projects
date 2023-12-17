import matplotlib.pyplot as plt
import numpy as np


# df/dx data
df_dx_7 = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_03/CovidRaw_akimadf_dx_7')
df_dx_10 = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_03/CovidRaw_akimadf_dx_10')

day7 = df_dx_7[:,0]
case_change7 = df_dx_7[:,1]

day10 = df_dx_10[:,0]
case_change10 = df_dx_10[:,1]

plt.plot(day7,case_change7,label='covid change rate/7')
#plt.plot(day10,case_change10,label='covid change rate/10')


# infection rate
# df/dx data
d2f_dx2_7 = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_03/CovidRaw_akimad2f_dx2_7')
d2f_dx2_10 = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_03/CovidRaw_akimad2f_dx2_10')

x1 = d2f_dx2_7[:,0]
y1 = d2f_dx2_7[:,1]

x2 = d2f_dx2_10[:,0]
y2 = d2f_dx2_10[:,1]

plt.plot(x1,y1,label='covid infection rate/7')
#plt.plot(x2,y2,label='covid infeection rate/10')

# averrage
avg7 = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_03/CovidRaw_akimaavg7')
avg10 = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_03/CovidRaw_akimaavg10')

x1 = avg7[:,0]
y1 = avg7[:,1]

x2 = avg10[:,0]
y2 = avg10[:,1]

plt.plot(x1,y1,label='covid case running 7')
#plt.plot(x2,y2,label='covid case running 10')
plt.gca().xaxis.set_ticks_position('top')
plt.xticks([1, 34, 41, 116], ["April 3rd", "May 8th", "May 15th", "July 29th"],rotation=90)


#labelS
plt.xlabel('days')
plt.ylabel('cases')

#title
plt.title(r'COVID 19')

plt.legend()
plt.grid(True)

plt.show()