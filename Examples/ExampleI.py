
if __name__ == '__main__':
    import json
    import numpy as np
    import matplotlib.pyplot as plt
    from FEM.Elasticity3D import Elasticity
    from FEM.Geometry.Geometry import Geometry3D
    from FEM.Utils import testNeighborg, plot_list_elements
    FILENAME = "big_esfera.json"
    geometria = Geometry3D.importJSON(FILENAME, fast=True, max_depth=10)
    p = geometria.elements[0]
    be = geometria.detectBorderElements()
    geometria.exportJSON("BIG_ESFERA_CON_BORDES.json")
    be = [geometria.elements[i] for i in be]
    print(len(be))
    plot_list_elements(geometria.elements)
    plot_list_elements(be, c="red", acum=True)
    plt.show()
