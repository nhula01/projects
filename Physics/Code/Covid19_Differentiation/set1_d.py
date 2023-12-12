import matplotlib.pyplot as plt
import numpy as np

root = np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_03/afx.dat')

a = root[:,0]
value = root[:,1]

plt.plot(a,value,label='root(a)')
#labelS
plt.xlabel('a [arb. unit]')
plt.ylabel('roots [arb. unit]')

#title
plt.title(r'Root locations with respect to a values')

plt.legend(loc='lower left')
plt.grid(True)

plt.show()


