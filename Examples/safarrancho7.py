if __name__ == '__main__':
    from FEM import MGDCM, ContinumTotalLagrangian, QuadMembraneLinear, QuadShellLinear, Geometry3D
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib.animation as animation

    # Example usage. Quadrilateral shell element diagonal
    coords = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 1], [0, 1, 1]])

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # plt.plot(*coords.T, 'ro-')
    # plt.title('Quadrilateral Shell Element')
    # plt.axis('equal')
    # plt.show()

    gdl = np.array([[0, 1, 2, 3], [4, 5, 6, 7], [
                   4, 5, 6, 7], [4, 5, 6, 7], [4, 5, 6, 7]])
    element = QuadShellLinear(coords, gdl)

    element.set_thickness(0.1)
    element.get_jacobians(0.5)
    element.calculate_deformation_gradients(0.5)
    # element.set_constitutive_model(cm)
    a = 0
