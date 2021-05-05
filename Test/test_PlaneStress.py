import unittest
import os
import numpy as np
import sys
import inspect
currentdir = os.path.dirname(os.path.abspath(
    inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
from FEM.PlaneStrain import PlaneStress
from FEM.Mesh.Geometry import Geometry
import matplotlib.pyplot as plt


FILENAME = 'Test/resources/beam_ws.msh'
TOL = 1*10**(-5)


class TestPlaneStress(unittest.TestCase):
    def test_cantilever_beam_uniform_1(self):
        geometry = Geometry.loadmsh(FILENAME)
        cb = geometry.cbFromSegment(3,0.0,1)
        cb += geometry.cbFromSegment(3,0.0,2)
        geometry.cbe = cb
        geometry.maskFromSegments()
        self.assertTrue(len(geometry.cbe)>0)
        h = 0.5
        L = 2.0
        E = 21000000
        v = 0.2
        b = 0.3
        I = b*h**3/12
        A = b*h
        k = 5/6
        G = E/(2*(1+v))
        gamma = 23.54
        W = A*gamma
        O = PlaneStress(geometry,E,v,b,fy=lambda x: -gamma)
        O.solve(plot=False)
        X,U = O.profile([0,0.5*h], [L,0.5*h])
        U = np.array(U)[:,1]
        g = lambda x: W*x**2/24/E/I*(x**2+6*L**2-4*L*x)+W*(L*x-x**2/2)/k/A/G
        plt.close()
        err = U + g(np.array(X))
        errores = np.max(err**2)
        self.assertTrue(errores<TOL)


    def test_cantilever_beam_uniform_2(self):
        geometry = Geometry.loadmsh(FILENAME)
        cb = geometry.cbFromSegment(3,0.0,1)
        cb += geometry.cbFromSegment(3,0.0,2)
        geometry.cbe = cb
        geometry.maskFromSegments()
        self.assertTrue(len(geometry.cbe)>0)
        h = 0.5
        L = 2.0
        E = 21000000
        v = 0.2
        b = 0.3
        I = b*h**3/12
        A = b*h
        gamma = 23.54
        k = 5/6
        G = E/(2*(1+v))
        W = A*gamma
        O = PlaneStress(geometry,E,v,b)
        O.geometry.loadOnSegment(2,fy= lambda x: -W)
        O.solve(plot=False)
        X,U = O.profile([0,0.5*h], [L,0.5*h])
        U = np.array(U)[:,1]
        g = lambda x: W*x**2/24/E/I*(x**2+6*L**2-4*L*x)+W*(L*x-x**2/2)/k/A/G
        plt.close()
        err = U + g(np.array(X))
        errores = np.max(err**2)
        self.assertTrue(errores<TOL)

    def test_cantilever_beam_triangular_3(self):
        geometry = Geometry.loadmsh(FILENAME)
        cb = geometry.cbFromSegment(3,0.0,1)
        cb += geometry.cbFromSegment(3,0.0,2)
        geometry.cbe = cb
        geometry.maskFromSegments()
        self.assertTrue(len(geometry.cbe)>0)
        h = 0.5
        L = 2.0
        E = 21000000
        v = 0.2
        b = 0.3
        I = b*h**3/12
        A = b*h
        gamma = 23.54
        W = A*gamma
        G = E/(2*(1+v))
        k = 5/6
        O = PlaneStress(geometry,E,v,b)
        O.geometry.loadOnSegment(2,fy= lambda x: -W*(x/L))
        O.solve(plot=False)
        X,U = O.profile([0,0.5*h], [L,0.5*h])
        U = np.array(U)[:,1]
        g = lambda x: W*x**2/120/E/I/L*(10*L**3-10*L**2*x+5*L*x**2-x**3)+W*(L*x/2-x**3/6/L)/k/A/G
        plt.close()
        err = U + g(np.array(X))
        errores = np.max(err**2)
        plt.plot(X,U,'--',color='black',label='FEM')
        plt.plot(X,-g(np.array(X)),color='black',label='Exact')
        plt.grid()
        plt.legend()
        plt.show()
        self.assertTrue(errores<TOL)


    def test_cantilever_beam_point_1(self):
        h = 0.5
        L = 2.0
        E = 21000000
        v = 0.2
        b = 0.3
        I = b*h**3/12
        A = b*h
        gamma = 23.54
        W = A*gamma
        G = E/(2*(1+v))
        k = 5/6

        P = 100

        geometry = Geometry.loadmsh(FILENAME)
        cb = geometry.cbFromSegment(3,0.0,1)
        cb += geometry.cbFromSegment(3,0.0,2)
        geometry.cbe = cb
        
        geometry.maskFromSegments()
        self.assertTrue(len(geometry.cbe)>0)

        O = PlaneStress(geometry,E,v,b)
        O.geometry.loadOnSegment(1,fy= lambda x: -P/A*b)
        O.solve(plot=False)
        X,U = O.profile([0,0.5*h], [L,0.5*h])
        U = np.array(U)[:,1]
        g = lambda x: P*x**2*(3*L-x)/6/E/I+P*x/k/A/G
        plt.close()
        err = U + g(np.array(X))
        errores = np.max(err**2)
        self.assertTrue(errores<TOL)

    def test_cantilever_beam_point_2(self):
        h = 0.5
        L = 2.0
        E = 21000000
        v = 0.2
        b = 0.3
        I = b*h**3/12
        A = b*h
        gamma = 23.54
        W = A*gamma
        G = E/(2*(1+v))
        k = 5/6

        P = 100

        geometry = Geometry.loadmsh(FILENAME)
        cb = geometry.cbFromSegment(3,0.0,1)
        cb += geometry.cbFromSegment(3,0.0,2)
        cbn = geometry.generateBCFromCoords(L,h,-P,2)
        geometry.cbe = cb
        geometry.cbn = cbn
        
        geometry.maskFromSegments()
        self.assertTrue(len(geometry.cbe)>0)

        O = PlaneStress(geometry,E,v,b)
        O.solve(plot=False)
        X,U = O.profile([0,0.5*h], [L,0.5*h])
        U = np.array(U)[:,1]
        g = lambda x: P*x**2*(3*L-x)/6/E/I+P*x/k/A/G
        plt.close()
        err = U + g(np.array(X))
        errores = np.max(err**2)
        self.assertTrue(errores<TOL)

if __name__ == '__main__':
    unittest.main()