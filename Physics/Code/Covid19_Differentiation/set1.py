import matplotlib.pyplot as plt
import numpy as np

data_forward = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_03/df_dx_forward.dat')
data_backward = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_03/df_dx_backward.dat')
data_central = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_03/df_dx_central.dat')
data_thirdorder = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_03/df_dx_thirdorder.dat')
data_fourthorder= np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_03/df_dx_fourthorder.dat')
data_exact= np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_03/df_dx_exact.dat')
root = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_03/a_1.dat')
data_akima= np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_03/df_dx_akima.dat')
#forward

def function(x):
    return -np.sin(x) * np.sinh(x)+np.cos(x) * np.cosh(x)


x = data_forward[:,0]

df_dx = data_forward[:,1]
#backward
x1 = data_backward[:,0]
df_dx1 = data_backward[:,1]
#central
x2 = data_central[:,0]
df_dx2 = data_central[:,1]

#thirdorder
x3 = data_thirdorder[:,0]
df_dx3 = data_thirdorder[:,1]
#fourthorder
x4 = data_fourthorder[:,0]
df_dx4 = data_fourthorder[:,1]
#exacr
x5 = data_exact[:,0]
df_dx5 = data_exact[:,1]
#akima
x6 = data_akima[:,0]
df_dx6 = data_akima[:,1]

plt.plot(x,df_dx,label='forward')
plt.plot(x1,df_dx1,label='backward')
plt.plot(x2,df_dx2,label='central')
plt.plot(x3,df_dx3,label='third order')
plt.plot(x4,df_dx4,label='fourth order')
plt.plot(x5,df_dx5,label='exact')
plt.scatter(root, function(root), s=100, alpha=1, label='root')
plt.plot(x6,df_dx6,label='akima')
#labelS
plt.xlabel('x [arb. unit]')
plt.ylabel('dy/dx [arb. unit]')

#title
plt.title(r'Derivative')

plt.legend(loc='lower left')
plt.grid(True)

plt.show()
