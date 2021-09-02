import numpy as np
import matplotlib.pyplot as plt
import math

num_of_basis_func = 4

x1 = 0.
x2 = 1.

h = x2 - x1

delta = 101

x_vec = np.zeros((delta))
N = np.zeros((delta, num_of_basis_func))

for i in range(delta):
    x = (h / (delta - 1.)) * i

    x_vec[i] = x

    N[i][0] = ((-((x - x2) ** 2.))*(-h + 2.*(x1 - x))) / (h ** 3.)

    N[i][1] = ((x - x1)*((x - x2) ** 2.)) / (h ** 2.)

    N[i][2] = (((x - x1) ** 2.)*(h + 2.*(x2 - x))) / (h ** 3.)

    N[i][3] = (((x - x1) ** 2.)*(x - x2)) / (h ** 2.)

    # if i == int(delta/2):
    #     N_total = N[i][0] + N[i][1] + N[i][2] + N[i][3]
    #     print(N_total)

    N_total = N[i][0] + N[i][1] + N[i][2] + N[i][3]
    N_total_2 = 2 * (x**3) - 3 * (x**2) + x + 1.

    print(N_total - N_total_2)

# Define color vector
color = np.array(["r", "g", "b", "c", "m", "y", "k"])

# 描写
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)

for i in range(num_of_basis_func):
    ax1.plot(x_vec[:], N[:, i], c=color[i %
             (color.shape[0])], marker="", linewidth=0.7)

# 描写
ax1.set_aspect('equal', adjustable='box')
ax1.set_axisbelow(True)
fig.set_figheight(9)
fig.set_figwidth(12)
ax1.grid()
ax1.set_xlim(-.5, 1.5)
ax1.set_ylim(-.5, 1.5)
plt.show()
