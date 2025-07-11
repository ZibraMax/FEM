import matplotlib.pyplot as plt
import numpy as np


a = 400
b = 20
L = (a**2 + b**2)**0.5
young = 20500  # Young's modulus in MPa
h = 1
t = 653
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


fact = EA/1000
U = np.linspace(0, 3*b, 100)
F = np.array([get_force(u) for u in U])/fact
F2 = np.array([get_force2(u) for u in U])/fact
F3 = np.array([get_force3(u) for u in U])/fact
F4 = np.array([get_force4(u) for u in U])/fact
# np.savetxt('./Examples/examples_results/roof_analytical.csv', np.array([U, F]).T,
#    header='Displacement\tLoad', delimiter=',')
U = U/b
fig = plt.figure(figsize=(6, 6))
plt.plot(U, F, label='Green Lagrange strain')
plt.plot(U, F2, label='Engineering strain')
plt.plot(U, F3, label='"True" strain')
plt.plot(U, F4, label='Logarithmic strain')

plt.xlabel('Displacement')
plt.ylabel('Load')
plt.xlim(0, 3)
plt.ylim(-0.1, 0.1)
plt.savefig('roof_analytical.png', dpi=300, transparent=True)
plt.grid()
plt.legend()
plt.show()
