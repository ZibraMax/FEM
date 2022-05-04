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

From source:

```bash
git clone https://github.com/ZibraMax/FEM
cd FEM
python -m venv .venv
python -m pip install build
python -m build
python -m pip install -e .[docs] # Basic instalation with docs
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

-   Create geometry
-   Create Border Conditions (Point and regions supported)
-   Solve!
-   For example: Example 2, Example 5, Example 11-14

#### Example without geometry file (Test 2):

```python
import matplotlib.pyplot as plt #Import libraries
from FEM.Torsion2D import Torsion2D #import AFEM Torsion class
from FEM.Geometry import Delaunay #Import Meshing tools

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

#Save geometry to .json file
geometry.exportJSON('I_test.json')

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
from FEM.Geometry import Geometry #Import Geometry tools

#Define material constants.
E = 200000
v = 0.27
G = E / (2 * (1 + v))
phi = 1 #Rotation angle

#Load geometry with file.
geometry = Geometry.importJSON('I_test.json')

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

1. 2D elastic plate theory
2. Transient analysis (Core modification)
3. Non-Lineal for 2D equation (All cases)
4. Testing and numerical validation (WIP)

## [Examples](https://zibramax.github.io/FEM/Examples.html)

## References

J. N. Reddy. Introduction to the Finite Element Method, Third Edition (McGraw-Hill Education: New York, Chicago, San Francisco, Athens, London, Madrid, Mexico City, Milan, New Delhi, Singapore, Sydney, Toronto, 2006). https://www.accessengineeringlibrary.com/content/book/9780072466850

Jonathan Richard Shewchuk, (1996) Triangle: Engineering a 2D Quality Mesh Generator and Delaunay Triangulator

Ramirez, F. (2020). ICYA 4414 Modelación con Elementos Finitos [Class handout]. Universidad de Los Andes.

## License

[MIT](https://choosealicense.com/licenses/mit/)
