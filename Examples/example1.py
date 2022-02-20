if __name__ == '__main__':
    # Creación de elementos
    import numpy as np

    from FEM.Elements.E2D.Serendipity import Serendipity
    from FEM.Elements.E2D.Quadrilateral import Quadrilateral
    from FEM.Elements.E2D.LTriangular import LTriangular
    from FEM.Elements.E2D.QTriangular import QTriangular

    from FEM.Elements.E1D.LinealElement import LinealElement
    from FEM.Elements.E1D.QuadraticElement import QuadraticElement

    ELEMENTOS = []

    coords = [[1, 1], [3, 2], [3.5, 3], [0.5, 4], [2, 1.5],
              [3.25, 2.5], [(3.5+.5)/2, 3.5], [(0.5+1)/2, 5/2]]
    gdl = np.array([[1, 2, 3, 4, 5, 6, 7, 8]])
    eRS = Serendipity(coords, gdl)
    eRS2 = Serendipity(coords, gdl)
    ELEMENTOS.append(eRS)

    coords = [[1, 1], [3, 2], [3.5, 3], [0.5, 4]]
    gdl = np.array([[1, 2, 3, 4]])
    eRC = Quadrilateral(coords, gdl)
    ELEMENTOS.append(eRC)

    # [[-1.5, -1], [1, 0], [0, 5], [-0.25, -0.5], [0.5, 2.5], [-0.75, 2.0]]
    coords = [[-1.5, -1], [1, 0], [0, 5]]
    for i in range(len(coords)-1):
        coords += [[(coords[i][0]+coords[i+1][0])/2,
                    (coords[i][1]+coords[i+1][1])/2]]
    coords += [[(coords[i+1][0]+coords[0][0])/2,
                (coords[i+1][1]+coords[0][1])/2]]
    gdl = np.array([[1, 2, 3, 4, 5, 6]])
    eTC = QTriangular(coords, gdl)
    ELEMENTOS.append(eTC)

    coords = [[1.4, 1.5], [3.1, 2.4], [2, 3]]
    gdl = np.array([[1, 2, 3]])
    eTL = LTriangular(coords, gdl)
    ELEMENTOS.append(eTL)

    coords = [[3], [4], [5]]
    gdl = np.array([[3, 4, 5]])
    e = QuadraticElement(coords, gdl, 3)
    ELEMENTOS.append(e)

    coords = [[3], [5]]
    gdl = np.array([[3, 5]])
    e = LinealElement(coords, gdl, 3)
    ELEMENTOS.append(e)

    # Dibujar elementos
    # for e in ELEMENTOS:
    # 	e.draw()
    # 	plt.show()

    # Esta Adentro
    # print(e.isInside(np.array([1,2,3,4,5])))

    # Inverse Mapping
    z = eTL.inverseMapping(np.array([[1.3, 2.5, 3.5], [1.5, 2.6, 8.5]]))
    # z = eTL.inverseMapping(np.array([[1.3],[1.5]]))
    print(z)
    # print(eTL.isInside(np.array([3.5,2.5])))

    # Jacobian Graphs
    # for e in ELEMENTOS:
    # 	e.jacobianGraph()
    # 	plt.show()
