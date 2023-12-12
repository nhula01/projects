import plotly.express as px

import numpy as np
import pandas as pd

N_points_value = np.loadtxt("/Users/phihung/Documents/PHY/first/homework/hw_04/double_integration_3d.dat")
x = N_points_value[:,0]
y = N_points_value[:,1]
z = N_points_value[:,2]

data = {'x_column_name': x,
        'y_column_name': y,
        'z_column_name': z}

df = pd.DataFrame(data)

fig = px.scatter_3d(df, x='x_column_name', y='y_column_name', z='z_column_name')
fig.show()
