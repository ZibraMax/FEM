<!-- LTeX: language=en -->

[![Build status](https://github.com/ZibraMax/FEM/actions/workflows/python-publish.yml/badge.svg)](https://github.com/ZibraMax/FEM/actions/workflows/python-publish.yml)
[![Docs](https://github.com/ZibraMax/FEM/actions/workflows/docs.yml/badge.svg)](https://github.com/ZibraMax/FEM/actions/workflows/docs.yml)
[![PyPI version](https://badge.fury.io/py/AFEM.svg)](https://badge.fury.io/py/AFEM)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/854107ce95794d28beac5ea5c44e1dd2)](https://www.codacy.com/gh/ZibraMax/FEM/dashboard?utm_source=github.com&utm_medium=referral&utm_content=ZibraMax/FEM&utm_campaign=Badge_Grade)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/ZibraMax/FEM/blob/master/LICENSE)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![GitHub release](https://img.shields.io/github/v/release/ZibraMax/FEM.svg)](https://github.com/ZibraMax/FEM/releases/)

A [Python](https://www.python.org/) FEM implementation.

N dimensional FEM implementation for M variables per node problems.

## Installation

Use the package manager [pip](https://pypi.org/project/AFEM/) to install AFEM.

```bash
pip install AFEM
```

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## [Full Docs](https://zibramax.github.io/FEM/)

## Tutorial

### Using pre implemented equations

Avaliable equations:

-   1D 1 Variable ordinary diferential equation
-   1D 1 Variable 1D Heat with convective border
-   1D 2 Variable Euler Bernoulli Beams
-   1D 3 Variable Non-linear Euler Bernoulli Beams
-   2D 1 Variable Torsion
-   2D 1 Variable Poisson equation
-   2D 1 Variable second order PDE
-   2D 1 Variable 2D Heat with convective borders
-   2D 2 Variable Plane Strees
-   2D 2 Variable Plane Strees Orthotropic
-   2D 2 Variable Plane Strain
-   3D 3 variables per node isotropic elasticity

Numerical Validation:

-   [x] 1D 1 Variable ordinary diferential equation
-   [ ] 1D 1 Variable 1D Heat with convective border
-   [ ] 1D 2 Variable Euler Bernoulli Beams
-   [ ] 1D 3 Variable Non-linear Euler Bernoulli Beams
-   [x] 2D 1 Variable Torsion
-   [ ] 2D 1 Variable 2D Heat with convective borders
-   [x] 2D 2 Variable Plane Strees
-   [x] 2D 2 Variable Plane Strain

#### Steps:

-   Create geometry (From coordinates or GiD)
-   Create Border Conditions (Point and regions supported)
-   Solve!
-   For example: Example 2, Example 5, Example 11-14

#### Example without geometry file (Test 2):

```python
import matplotlib.pyplot as plt #Import libraries
from FEM.Torsion2D import Torsion2D #import AFEM Torsion class
from FEM.Mesh.Delaunay import Delaunay #Import Meshing tools

#Define some variables with geometric properties
a = 0.3
b = 0.3
tw = 0.05
tf = 0.05

#Define material constants
E = 200000
v = 0.27
G = E / (2 * (1 + v))
phi = 1 #Rotation angle

#Define domain coordinates
vertices = [
        [0, 0],
        [a, 0],
        [a, tf],
        [a / 2 + tw / 2, tf],
        [a / 2 + tw / 2, tf + b],
        [a, tf + b],
        [a, 2 * tf + b],
        [0, 2 * tf + b],
        [0, tf + b],
        [a / 2 - tw / 2, tf + b],
        [a / 2 - tw / 2, tf],
        [0, tf],
]

#Define triangulation parameters with `_strdelaunay` method.
params = Delaunay._strdelaunay(constrained=True, delaunay=True,
                                                                        a='0.00003', o=2)
#**Create** geometry using triangulation parameters. Geometry can be imported from .msh files.
geometry = Delaunay(vertices, params)

#Save geometry to .msh file
geometry.saveMesh('I_test')

#Create torsional 2D analysis.
O = Torsion2D(geometry, G, phi)
#Solve the equation in domain.
#Post process and show results
O.solve()
plt.show()

```

#### Example with geometry file (Example 2):

```python
import matplotlib.pyplot as plt #Import libraries
from FEM.Torsion2D import Torsion2D #import AFEM
from FEM.Mesh.Geometry import Geometry #Import Geometry tools

#Define material constants.
E = 200000
v = 0.27
G = E / (2 * (1 + v))
phi = 1 #Rotation angle

#Load geometry with file.
geometry = Geometry.loadmsh('I_test.msh')

#Create torsional 2D analysis.
O = Torsion2D(geometry, G, phi)
#Solve the equation in domain.
#Post process and show results
O.solve()
plt.show()


```

### Creating equation classes

Note: Don't forget the docstring!

#### Steps

1. Create a Python flie and import the libraries:

    ```python
    from .Core import *
    from tqdm import tqdm
    import numpy as np
    import matplotlib.pyplot as plt
    ```

    - Core: Solver
    - Numpy: Numpy data
    - Matplotlib: Matplotlib graphs
    - Tqdm: Progressbars

2. Create a Python class with Core inheritance
    ```python
    class PlaneStress(Core):
        def __init__(self,geometry,*args,**kargs):
        #Do stuff
        Core.__init__(self,geometry)
    ```
    It is important to manage the number of variables per node in the input geometry.
3. Define the matrix calculation methods and post porcessing methods.

    ```python
    def elementMatrices(self):
    def postProcess(self):
    ```

4. The `elementMatrices` method uses gauss integration points, so you must use the following structure:

    ```python

    for e in tqdm(self.elements,unit='Element'):
        _x,_p = e.T(e.Z.T) #Gauss points in global coordinates and Shape functions evaluated in gauss points
        jac,dpz = e.J(e.Z.T) #Jacobian evaluated in gauss points and shape functions derivatives in natural coordinates
        detjac = np.linalg.det(jac)
        _j = np.linalg.inv(jac) #Jacobian inverse
        dpx = _j @ dpz #Shape function derivatives in global coordinates
        for k in range(len(e.Z)): #Iterate over gauss points on domain
            #Calculate matrices with any finite element model
        #Assign matrices to element
    ```

A good example is the `PlaneStress` class in the `Elasticity2D.py` file.

## Roadmap

2. 2D elastic plate theory
3. Geometry class modification for hierarchy with 1D, 2D and 3D geometry child classes
4. Transient analysis (Core modification)
5. Non-Lineal for 2D equation (All cases)
6. Testing and numerical validation (WIP)

## Example index:

-   Example 1: Preliminar geometry test

-   Example 2: 2D Torsion 1 variable per node. H section-Triangular Quadratic.

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example2.png">

-   Example 3: 2D Torsion 1 variable per node. Square section-Triangular Quadratic.

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example3.png">

-   Example 4: 2D Torsion 1 variable per node. Mesh from internet-Square Lineal.

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example4.png">

-   Example 5: 2D Torsion 1 variable per node. Creating and saving mesh-Triangular Quadratic.

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example5.png">

-   Example 6: 1D random differential equation 1 variable per node. Linear Quadratic.

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example6.png">

-   Example 7: GiD Mesh import test — Serendipity elements

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example7.png">

-   Example 8: Plane Stress 2 variable per node. Plate in tension — Serendipity.

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example8.png">

-   Example 9: Plane Stress 2 variable per node. Simple Supported Beam — Serendipity.

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example9.png">

-   Example 10: Plane Stress 2 variable per node. Cantilever Beam — Triangular Quadratic.

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example10.png">

-   Example 11: Plane Stress 2 variable per node. Fixed-Fixed Beam — Serendipity.

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example11.png">

-   Example 12: Plane Strain 2 variable per node. Embankment from GiD — Serendipity.

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example12.png">

-   Example 13: Plane Strain 2 variable per node. Embankment — Triangular Quadratic.

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example13_geometry.png">
      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example13.png">

-   Example 14: Plane Stress 2 variable per node. Cantilever Beam — Serendipity.

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example14.png">

-   Example 15: Profile creation tool. Same as Example 14

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example15.png">

-   Example 16: Non-Local Plane Stress. [WIP]
-   Example 17: 1D Heat transfer.

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example17.png">

-   Example 18: 2D border elements creation.

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example18.png">

-   Example 19: Apply loads on regions. `loadOnRegion` method on Test 11

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example19.png">

-   Example 20: Reddy's Example 11.7.1 Ed 3
-   Example 21: Example 20 with serendipity elements.
-   Example 22: Example 20 with refined mesh.

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example22.png">

-   Example 23: Reddy's Problem 11.1 Ed 3 Plain Strain

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example23.png">

-   Example 24: Example 23 with refined mesh

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example24.png">

-   Example 25: Holes concept. With Example 24

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example25_geometry.png">
      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example25.png">

-   Example 26: Fillets concept.

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example26_geometry.png">
      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example26.png">

-   Example 27: Combination of Holes Fillets, Plane Stress

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example27_geometry.png">

-   Example 28: Fillets and Holes mesh files of Example 27

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example28.png">

-   Example 29: Fillets and Holes in Example 13

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example29_geometry.png">
      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example29.png">

-   Example 30: Border conditions and loads in holes

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example30_geometry.png">
      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example30.png">

-   Example 31: 2D Heat with convective borders

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example31_geometry.png">
      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example31.png">

-   Example 32: Border conditions and loads in holes

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example32_geometry.png">
      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example32.png">

-   Example 33: Example 30 with Heat

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example33_geometry.png">
      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example33.png">

-   Example 34: Custom plots, Beam-Girder steel plate connection

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example34_geometry.png">
      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example34.png">

-   Example 35: Torsion with fillets

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example35_geometry.png">
      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example35.png">

-   Example 36: Convective Heat Transfer from [Samson-Mano's software](https://github.com/Samson-Mano/2D_Heat_transfer)

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example36.png">

-   Example 37: Convective Heat Transfer from [Samson-Mano's software](https://github.com/Samson-Mano/2D_Heat_transfer)

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example37_geometry.png">
      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example37.png">

-   Example 38: Elements with different properties: Torsion with holes

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example38_geometry.png">
      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example38.png">

-   Example 37: Elements with different properties: Torsion with holes Symetrical

      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example39_geometry.png">
      <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example39.png">

-   Example 38 & 39: Polar moment of inertia for hollow sections
-   Example 40 & 41: Euler Bernoulli beams, linear and non-linear
-   Example 42: Non-linear equation solver test
-   Example 43: Orthotripic plane stress
-   Example 44: MeshingNet data creation code

## References

J. N. Reddy. Introduction to the Finite Element Method, Third Edition (McGraw-Hill Education: New York, Chicago, San Francisco, Athens, London, Madrid, Mexico City, Milan, New Delhi, Singapore, Sydney, Toronto, 2006). https://www.accessengineeringlibrary.com/content/book/9780072466850

Jonathan Richard Shewchuk, (1996) Triangle: Engineering a 2D Quality Mesh Generator and Delaunay Triangulator

Ramirez, F. (2020). ICYA 4414 Modelación con Elementos Finitos [Class handout]. Universidad de Los Andes.

## License

[MIT](https://choosealicense.com/licenses/mit/)
