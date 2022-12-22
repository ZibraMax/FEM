
if __name__ == '__main__':
    import time
    import matplotlib.pyplot as plt
    from FEM.Geometry.Geometry import Geometry3D
    from FEM.Utils import plot_list_elements
    FILENAME = "Examples\Mesh_tests\SPHERE_FERNANDO.json"
    E = 21000000.0
    v = 0.2
    h = 0.6
    b = 0.3
    L = 3.5
    gamma = 23.54

    geometria = Geometry3D.importJSON(FILENAME, fast=True)
    e = geometria.elements[1]
    #fig = plt.figure()
    #ax = fig.add_subplot(projection="3d")
    # geometria.OctTree.graph_query_range(
    #    e.T(e.center.T)[0].flatten(), 2*d)
    st = time.time()
    be = geometria.detectBorderElementsLegacy()
    et1 = time.time()
    be2 = geometria.detectBorderElements()
    et2 = time.time()
    print(len(be), et1-st, len(be2), et2-et1)
    geometria.exportJSON("ESFERA_CON_BORDES.json")
    be = [geometria.elements[i] for i in be]
    be2 = [geometria.elements[i] for i in be2]
    plot_list_elements(geometria.elements)
    plot_list_elements(be, c="red", acum=True)
    plt.show()

    plot_list_elements(geometria.elements)
    plot_list_elements(be2, c="green", acum=True)
    plt.show()
