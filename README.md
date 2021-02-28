# FEM
N dimensional FEM implementation for M variables per node problems.

## Tutorial

### Using pre implemented equations

Avaliable equations:
- 1D 1 Variable ordinary diferential equation
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

### Creating equation classes

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
3. Define the matrix calculation methods and post porcessing methods
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
3. 1D and 2D heat transfer
4. Geometry class modification for hierarchy with 1D, 2D and 3D geometry child classes
4. Transient analysis (Core modification)
5. Elasticity in 3D (3D meshing and post process)
6. Non Lineal analysis for 1D equation (All cases)
7. Non Lineal for 2D equation (All cases)
8. UNIT TESTING
9. NUMERICAL VALIDATION
10. Non Local 2D?

## Test index:

- Test 1: Preliminar geometry test
- Test 2: 2D Torsion 1 variable per node. H section - Triangular Quadratic
- Test 3: 2D Torsion 1 variable per node. Square section - Triangular Quadratic
- Test 4: 2D Torsion 1 variable per node. Mesh from internet - Square Lineal
- Test 5: 2D Torsion 1 variable per node. Creating and saving mesh - Triangular Quadratic 
- Test 6: 1D random differential equation 1 variable per node. Linear Quadratic
- Test 7: GiD Mesh import test - Serendipity elements
- Test 8: Plane Stress 2 variable per node. Plate in tension - Serendipity
- Test 9: Plane Stress 2 variable per node. Simple Supported Beam - Serendipity
- Test 10: Plane Stress 2 variable per node. Cantilever Beam - Triangular Quadratic
- Test 11: Plane Stress 2 variable per node. Fixed-Fixed Beam - Serendipity
- Test 12: Plane Strain 2 variable per node. Embankment from GiD - Serendipity
- Test 13: Plane Strain 2 variable per node. Embankment - Triangular Quadratic
- Test 14: Plane Stress 2 variable per node. Cantilever Beam - Serendipity

## References

J. N. Reddy. Introduction to the Finite Element Method, Third Edition (McGraw-Hill Education: New York, Chicago, San Francisco, Athens, London, Madrid, Mexico City, Milan, New Delhi, Singapore, Sydney, Toronto, 2006). https://www.accessengineeringlibrary.com/content/book/9780072466850

Jonathan Richard Shewchuk, (1996) Triangle: Engineering a 2D Quality Mesh Generator and Delaunay Triangulator
