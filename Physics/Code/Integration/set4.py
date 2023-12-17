#import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
"""N_points_value = np.loadtxt("/Users/phihung/Documents/PHY/first/homework/hw_04/double_integration.dat")
y = N_points_value[:,0]
value = N_points_value[:,1]
plt.plot(y, value, color='y')
plt.xlabel('y [arb.unit]')
plt.ylabel('value [arb.unit]')
plt.title('value vs y')
plt.legend(loc='lower right')
plt.grid(True)"""
y_values = [1, 2, 3, 4, 5]
x_values = [1, 2, 3, 4, 5]
z_values = [10, 8, 6, 4, 2]

fig = go.Figure(data=[go.Scatter3d(x=x_values, y=y_values, z=z_values, mode='markers')])
fig.update_layout(scene=dict(xaxis_title='X Axis', yaxis_title='Y Axis', zaxis_title='Z Axis'))
fig.update_layout(title='3D Scatter Plot with Y Values')

fig.show()
#plt.show()