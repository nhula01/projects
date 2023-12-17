import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

min_val, max_val = 1, 3

intersection_matrix = np.array([[7, 6, 5],
                                [-2, 0, 2],
                                [-5, -6, -7]])

ax.matshow(intersection_matrix, cmap=plt.cm.Blues)

for i in range(3):
    for j in range(3):
        c = intersection_matrix[j,i]
        ax.text(i, j, str(c), va='center', ha='center')
plt.show()
