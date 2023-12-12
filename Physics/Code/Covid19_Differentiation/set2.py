import matplotlib.pyplot as plt
import numpy as np

data_fourthorder=np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_03/df_dx_fourthorder_problem2.dat')
data_exact=np.loadtxt('/Users/phihung/Documents/PHY/first/homework/hw_03/df_dx_exact_problem2.dat')

#fourth order
x = data_fourthorder[:,0]
df_dx = data_fourthorder[:,1]

#exact
x1 = data_exact[:,0]
df_dx1 = data_exact[:,1]

plt.plot(x,df_dx,label='fourth_order')
plt.plot(x1,df_dx1,label='exact')

#labelS
plt.xlabel('x [arb. unit]')
plt.ylabel('dy/dx [arb. unit]')

#title
plt.title(r'Derivative')

plt.legend()
plt.grid(True)

plt.show()