import FEM
import numpy as np
import matplotlib.pyplot as plt

coords = np.array([[0, 0, 0.0], [1, 1, 0.5], [2, 2, 2]])
gdl = np.array([[1, 2, 3, 4, 5, 6], [7, 8, 9, 10, 11, 12],
               [13, 14, 15, 16, 17, 18]]).T
e = FEM.QuadraticElement(coords, gdl, 8, boundary=True)

# Plot 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
coords = e.T(e.domain.T)[0].T
jacs = e.J(e.domain.T)[0]

ax.plot(*coords, 'o-')
for base, jac in zip(coords.T, jacs):

    ax.plot(*coords, 'o-')
    ax.quiver(*base, *jac[0]/3, color="red")
plt.show()
