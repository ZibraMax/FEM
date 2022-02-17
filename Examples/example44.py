# %% Functions
from FEM.Geometry.Delaunay import Delaunay
from FEM.Torsion2D import Torsion2D
from FEM.Utils.polygonal import generatePolygon
import numpy as np
import triangle as tr
import matplotlib.pyplot as plt


def generatePair(cx=10.0, cy=10.0, r=5.0, irre=0.5, spik=0.3, vert=8, areaLDUM=0.5, areaHDUM=0.0425, bc=None, valuebc=0.0):
    """Generates a random geometry pair. Low Density Uniform Mesh (LDUM) and High Density Uniform Mesh (HDUM)

    Args:
        cx (float, optional): X centroid of genometry. Defaults to 10.0.
        cy (float, optional): Y centroid of geometry. Defaults to 10.0.
        r (float, optional): Polygon equivalent radius. Defaults to 5.0.
        irre (float, optional): Irregularity of polygon [0,1]. Defaults to 0.5.
        spik (float, optional): Spikeness of polygon [0,1]. Defaults to 0.3.
        vert (int, optional): Number of vertices of polygon. Defaults to 8.
        areaLDUM (float, optional): Maximun trinagle area for the LDUM. Defaults to 0.5.
        areaHDUM (float, optional): Maximun trinagle area for the HDUM. Defaults to 0.0425.
        bc (list or None, optional): List of borders in which border conditions will be applied. If none, the border conditions is applied in all borders. Defaults to None.
        valuebc (float, optional): Border condition value. Defaults to 0.0.

    Returns:
        Tuple[Torsion2D]: LDUM and HDUM as Torsion2D objects.
    """
    polygon = generatePolygon(cx, cy, r, irre, spik, vert)
    paramsLD = Delaunay._strdelaunay(a=areaLDUM, o=1)
    paramsHD = Delaunay._strdelaunay(a=areaHDUM, o=1)
    LDUMgeometry = Delaunay(polygon, paramsLD)
    LDUMgeometry.maskFromSegments()
    HDUMgeometry = Delaunay(polygon, paramsHD)
    HDUMgeometry.maskFromSegments()
    G = 0.5  # -∇/G²Ψ = 2ϕ  #If you want to set the poisson equation to a specific value you can change the ϕ value
    phi = 1
    LDUM = Torsion2D(LDUMgeometry, G, phi)
    HDUM = Torsion2D(HDUMgeometry, G, phi)
    LDUM.geometry.cbe = []
    HDUM.geometry.cbe = []
    if bc:
        for s in bc:
            LDUM.geometry.cbe += LDUM.geometry.cbFromSegment(s, valuebc)
            HDUM.geometry.cbe += HDUM.geometry.cbFromSegment(s, valuebc)
    else:
        for s in range(len(LDUM.geometry.regions)):
            LDUM.geometry.cbe += LDUM.geometry.cbFromSegment(s, valuebc)
            HDUM.geometry.cbe += HDUM.geometry.cbFromSegment(s, valuebc)
    return LDUM, HDUM


def error(LDUM: Torsion2D, HDUM: Torsion2D):
    """Calculates the error between two Torsion2D objetcs

    Args:
        LDUM (Torsion2D): The LDUM
        HDUM (Torsion2D): The HDUM

    Returns:
        np.ndarray: An array with the calculated error in the LDUM nodes
    """
    print('Solving LDUM')
    LDUM.solve(plot=False)
    print('Solving HDUM')
    HDUM.solve(plot=False)
    errors = []
    for i in range(HDUM.ngdl):
        coord = HDUM.geometry.gdls[i]
        uHDUM = HDUM.U[i, 0]
        for e in LDUM.elements:
            if e.isInside(coord):
                uLDUM = inverseMappingSolution(e, *coord)
                uLDUM = uLDUM[0]
                if abs(uHDUM) == 0:
                    errors += [0]
                else:
                    errors += [np.abs((uLDUM-uHDUM)/uHDUM)]
                break
    return np.array(errors)


def inverseMappingSolution(e, x, y):
    """Calculates the interpolated solution in a given point. This functions only works with lineal triangular elements

    Args:
        e (Element): Elemento for interpolate the solution
        x (float): X coordinate
        y (float): Y coordinate

    Returns:
        float: The solution ibnterpolate value
    """
    corners = np.array(e.coords)
    e.alpha = [corners[1][0]*corners[2][1]-corners[2][0]*corners[1][1], corners[2][0]*corners[0]
               [1]-corners[0][0]*corners[2][1], corners[0][0]*corners[1][1]-corners[1][0]*corners[0][1]]
    e.beta = [corners[1][1]-corners[2][1], corners[2]
              [1]-corners[0][1], corners[0][1]-corners[1][1]]
    e.gamma = [-(corners[1][0]-corners[2][0]), -(corners[2][0] -
                                                 corners[0][0]), -(corners[0][0]-corners[1][0])]
    e.area2 = np.sum(e.alpha)
    e.psisLocales = lambda x, y: np.array([1/e.area2*(e.alpha[0]+e.beta[0]*x+e.gamma[0]*y),
                                          1/e.area2 *
                                           (e.alpha[1]+e.beta[1]
                                            * x+e.gamma[1]*y),
                                          1/e.area2*(e.alpha[2]+e.beta[2]*x+e.gamma[2]*y)]).T
    return e.Ue@e.psisLocales(x, y).T


def areaCorrection(errors, K=0.5, alpha=1, area0=1.67):
    """Calculates the new area based on lineal model of the paper

    Args:
        errors (np.ndarray): Errors in the LDUM nodes
        K (float, optional): K parameter. Defaults to 0.5.
        alpha (int, optional): Alpha value. Defaults to 1.
        area0 (float, optional): Maximun area value. Defaults to 1.67.

    Returns:
        np.ndarray: Area value in the LDUM nodes
    """
    Areas = K*1/errors**alpha * \
        (K/errors**alpha <= area0)+(K/errors**alpha > area0)*area0
    Areas = np.nan_to_num(Areas, nan=area0)
    return Areas


def areaRefiner(LDUM):
    """Extract the min of area for each element.


    Args:
        LDUM (Torsion2D): The LDUM

    Returns:
        list: The list of the new areas of the triangulation
    """
    a = []
    for e in LDUM.elements:
        a.append(np.min(e.Ue))
    return a


def extractTrainDataset(model, areasR):
    """Generates the dataset over a single model. One row for element

    Args:
        model (Torsion2D): The LDUM
        areasR (list): The new areas generated by the areaRefiner function

    Returns:
        list: Dataset matrix
    """
    m = []
    for i, e in enumerate(model.elements):
        coordsModel = model.geometry.original['vertices'].flatten().tolist()
        cx = [np.average(e.coords[:, 0])]
        cy = [np.average(e.coords[:, 1])]
        coords = e.coords.flatten().tolist()
        cb = [np.any(np.isin(e.gdl, np.array(model.cbe)[:, 0]))*1]
        y = [areasR[i]]
        fila = coordsModel+coords+cx+cy+cb+y
        m.append(fila)
    return m


# %% FEM structure generation
LDUM, HDUM = generatePair()
# %% Area correction calculation
errors = error(LDUM, HDUM)
errors.reshape([len(errors), 1])
LDUM.geometry.mask = None
LDUM.solveFromArray(errors, derivatives=False)
areas = areaCorrection(errors, K=0.02, alpha=1, area0=0.5)
LDUM.solveFromArray(areas, derivatives=False, plot=False)
plt.show()  # Delete thsi line if you gon on for.
# %% Generate new triangulation

areasR = areaRefiner(LDUM)
LDUM.geometry.triangulation['triangle_max_area'] = np.array(areasR).reshape([
    len(areasR), 1])
tnueva = tr.triangulate(LDUM.geometry.triangulation, 'ra')
tr.compare(plt, LDUM.geometry.triangulation, tnueva)

# %% Save the dataset results in csv file and png file.
plt.savefig('PROBLEMA.png')
np.savetxt("PROBLEMA.csv",
           np.array(extractTrainDataset(LDUM, areasR)), delimiter=",", fmt='%s')
plt.close('all')
