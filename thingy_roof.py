import matplotlib.pyplot as plt
import numpy as np


a = 1
b = 1
L = (a**2 + b**2)**0.5
young = 210e3  # Young's modulus in MPa
h = 0.5
t = 0.1
A = h*t
EA = young*A

sintheta = b/L


def material_model(E):
    Sx = young*E
    return Sx


def get_force(delta):
    LP = (L**2 + delta**2 - 2*L*sintheta*delta)**0.5
    sinthetap = (L*sintheta-delta)/LP
    eng_strain = (L-LP)/L
    E = eng_strain + 1/2*eng_strain**2
    stress = material_model(E)
    N = stress*A
    P = 2*N*sinthetap
    return P


def get_force2(delta):
    LP = (L**2 + delta**2 - 2*L*sintheta*delta)**0.5
    sinthetap = (L*sintheta-delta)/LP
    eng_strain = (L-LP)/L
    E = eng_strain
    stress = material_model(E)
    N = stress*A
    P = 2*N*sinthetap
    return P


def get_force3(delta):
    LP = (L**2 + delta**2 - 2*L*sintheta*delta)**0.5
    sinthetap = (L*sintheta-delta)/LP
    eng_strain = (L-LP)/L
    E = np.log(eng_strain + 1)
    stress = material_model(E)
    N = stress*A
    P = 2*N*sinthetap
    return P


def get_force4(delta):
    LP = (L**2 + delta**2 - 2*L*sintheta*delta)**0.5
    sinthetap = (L*sintheta-delta)/LP
    E = np.log(LP/L)
    stress = material_model(E)
    N = stress*A
    P = 2*N*sinthetap
    return -P


U = np.linspace(0, 8, 100)
F = np.array([get_force(u) for u in U])
F2 = np.array([get_force2(u) for u in U])
F3 = np.array([get_force3(u) for u in U])
F4 = np.array([get_force4(u) for u in U])
# np.savetxt('./Examples/examples_results/roof_analytical.csv', np.array([U, F]).T,
#    header='Displacement\tLoad', delimiter=',')
plt.plot(U, F, label='Green Lagrange strain')
plt.plot(U, F2, label='Engineering strain')
plt.plot(U, F3, label='"True" strain')
plt.plot(U, F4, label='Logarithmic strain')

plt.xlabel('Displacement')
plt.ylabel('Load')
plt.title('Analytical solution')
plt.grid()
plt.legend()
plt.show()
