<!-- LTeX: language=en -->
[![Build status](https://github.com/ZibraMax/FEM/actions/workflows/python-publish.yml/badge.svg)](https://github.com/ZibraMax/FEM/actions/workflows/python-publish.yml)
[![Docs](https://github.com/ZibraMax/FEM/actions/workflows/docs.yml/badge.svg)](https://github.com/ZibraMax/FEM/actions/workflows/docs.yml)
[![PyPI version](https://badge.fury.io/py/AFEM.svg)](https://badge.fury.io/py/AFEM)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/854107ce95794d28beac5ea5c44e1dd2)](https://www.codacy.com/gh/ZibraMax/FEM/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ZibraMax/FEM&amp;utm_campaign=Badge_Grade)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/ZibraMax/FEM/blob/master/LICENSE)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![GitHub release](https://img.shields.io/github/v/release/ZibraMax/FEM.svg)](https://github.com/ZibraMax/FEM/releases/)

A [Python](https://www.python.org/) FEM implementation.

N dimensional FEM implementation for M variables per node problems.

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install AFEM.

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
- 1D 1 Variable ordinary diferential equation
- 1D 1 Variable 1D Heat with convective border
- 1D 2 Variable Euler Bernoulli Beams [TODO]
- 1D 2 Variable Timoshenko Beams [TODO]
- 2D 1 Variable Torsion
- 2D 2 Variable Plane Strees
- 2D 2 Variable Plane Strain

#### Steps:
- Create geometry (From coordinates or GiD)
- Create Border Conditions (Point and segment supported)
- Solve!
- For example: Test 2, Test 5, Test 11-14

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

#### Example with geometry file (Test 2):


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
    - Core: Numpy data
    - Core: Matplotlib graphs
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

A good example is the `PlaneStress` class

## Roadmap

1. Beam bending by Euler Bernoulli and Timoshenko equations
2. 2D elastic plate theory
3. 2D heat transfer
4. Geometry class modification for hierarchy with 1D, 2D and 3D geometry child classes
5. Transient analysis (Core modification)
6. Elasticity in 3D (3D meshing and post process)
7. Non-Lineal analysis for 1D equation (All cases)
8. Non-Lineal for 2D equation (All cases)
9. UNIT TESTING
10. NUMERICAL VALIDATION
11. Non-Local 2D?

## Test index:

- Test 1: Preliminar geometry test

- Test 2: 2D Torsion 1 variable per node. H section-Triangular Quadratic.

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test2.png">
- Test 3: 2D Torsion 1 variable per node. Square section-Triangular Quadratic.

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test3.png">
- Test 4: 2D Torsion 1 variable per node. Mesh from internet-Square Lineal.

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test4.png">
- Test 5: 2D Torsion 1 variable per node. Creating and saving mesh-Triangular Quadratic.

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test5.png">
- Test 6: 1D random differential equation 1 variable per node. Linear Quadratic.

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test6.png">
- Test 7: GiD Mesh import test — Serendipity elements

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test7.png">
- Test 8: Plane Stress 2 variable per node. Plate in tension — Serendipity.

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test8.png">
- Test 9: Plane Stress 2 variable per node. Simple Supported Beam — Serendipity.

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test9.png">
- Test 10: Plane Stress 2 variable per node. Cantilever Beam — Triangular Quadratic.

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test10.png">
- Test 11: Plane Stress 2 variable per node. Fixed-Fixed Beam — Serendipity.

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test11.png">
- Test 12: Plane Strain 2 variable per node. Embankment from GiD — Serendipity.

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test12.png">
- Test 13: Plane Strain 2 variable per node. Embankment — Triangular Quadratic.

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test13_geometry.png">
    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test13.png">
- Test 14: Plane Stress 2 variable per node. Cantilever Beam — Serendipity.

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test14.png">
- Test 15: Profile creation tool. Same as Test 14

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test15.png">
- Test 16: Non-Local Plane Stress. [WIP]
- Test 17: 1D Heat transfer.

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test17.png">

- Test 18: 2D border elements creation.

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test18.png">
- Test 19: Apply loads on segments. `loadOnSegment` method on Test 11

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test19.png">
- Test 20: Reddy's Example 11.7.1 Ed 3
- Test 21: Test 20 with serendipity elements.
- Test 22: Test 20 with refined mesh.

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test22.png">
- Test 23: Reddy's Problem 11.1 Ed 3 Plain Strain

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test23.png">
- Test 24: Test 23 with refined mesh

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test24.png">
- Test 25: Holes concept. With Test 24

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test25_geometry.png">
    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test25.png">
- Test 26: Fillets concept.

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test26_geometry.png">
    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test26.png">
- Test 27: Combination of Holes Fillets, Plane Stress

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test27_geometry.png">
- Test 28: Fillets and Holes mesh files of Test 27

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test28.png">
- Test 29: Fillets and Holes in Test 13

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test29_geometry.png">
    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test29.png">


- Test 30: Border conditions and loads in holes

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test30_geometry.png">
    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test30.png">

- Test 31: 2D Heat with convective borders

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test31_geometry.png">
    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test31.png">

- Test 32: Border conditions and loads in holes

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test32_geometry.png">
    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test32.png">

- Test 33: Test 30 with Heat

    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test33_geometry.png">
    <img src="https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/test_results/Test33.png">

## References

J. N. Reddy. Introduction to the Finite Element Method, Third Edition (McGraw-Hill Education: New York, Chicago, San Francisco, Athens, London, Madrid, Mexico City, Milan, New Delhi, Singapore, Sydney, Toronto, 2006). https://www.accessengineeringlibrary.com/content/book/9780072466850

Jonathan Richard Shewchuk, (1996) Triangle: Engineering a 2D Quality Mesh Generator and Delaunay Triangulator

## License
[MIT](https://choosealicense.com/licenses/mit/)