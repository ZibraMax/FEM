import numpy as np
import matplotlib.pyplot as plt
from ..Elements.E3D.Brick import Brick


class Quadrant3D(Brick):

    def __init__(self, p: tuple, dim: tuple) -> None:
        x, y, z = p
        w, h, d = dim
        self.x, self.y, self.z = x, y, z
        self.w, self.h, self.d = w, h, d
        coords = np.array([
            [x-w, y-h, z-d],
            [x+w, y-h, z-d],
            [x+w, y+h, z-d],
            [x-w, y+h, z-d],
            [x-w, y-h, z+d],
            [x+w, y-h, z+d],
            [x+w, y+h, z+d],
            [x-w, y+h, z+d]])
        Brick.__init__(self, coords, np.array(
            [[-1]*8]), border=True, fast=True)

    def contains(self, e: Brick) -> bool:
        x, _ = e.T(e.center.T)
        x = x.flatten()

        maximos_self = np.max(self.coords, axis=0)
        minimos_self = np.min(self.coords, axis=0)

        superior = (maximos_self-x) >= 0
        inferior = (x-minimos_self) >= 0

        return superior.all() and inferior.all()

    def boxes_disjoint(self, e):

        maximos_self = np.max(self.coords, axis=0)
        minimos_self = np.min(self.coords, axis=0)

        maximos_e = np.max(e.coords, axis=0)
        minimos_e = np.min(e.coords, axis=0)

        maxx1, maxy1, maxz1 = maximos_self
        minx1, miny1, minz1 = minimos_self

        maxx2, maxy2, maxz2 = maximos_e
        minx2, miny2, minz2 = minimos_e

        return (maxx2 <= minx1 or maxx1 <= minx2
                or maxy2 <= miny1 or maxy1 <= miny2
                or maxz2 <= minz1 or maxz1 <= minz2)

    def intesects_quadrant(self, e: Brick) -> bool:
        return not self.boxes_disjoint(e)

    def subdivide(self) -> list:
        divs = []
        nw = self.w/2
        nh = self.h/2
        nd = self.d/2
        x, y, z = self.x, self.y, self.z

        divs.append(Quadrant3D((x+nw, y+nw, z+nd), (nw, nh, nd)))
        divs.append(Quadrant3D((x+nw, y+nw, z-nd), (nw, nh, nd)))
        divs.append(Quadrant3D((x+nw, y-nw, z+nd), (nw, nh, nd)))
        divs.append(Quadrant3D((x+nw, y-nw, z-nd), (nw, nh, nd)))
        divs.append(Quadrant3D((x-nw, y+nw, z+nd), (nw, nh, nd)))
        divs.append(Quadrant3D((x-nw, y+nw, z-nd), (nw, nh, nd)))
        divs.append(Quadrant3D((x-nw, y-nw, z+nd), (nw, nh, nd)))
        divs.append(Quadrant3D((x-nw, y-nw, z-nd), (nw, nh, nd)))
        return divs

    def draw(self, ax):
        ax.plot([self.x, self.x], [self.y, self.y],
                [self.z-self.d, self.z+self.d], c="k", alpha=0.4)
        ax.plot([self.x-self.w, self.x+self.w],
                [self.y, self.y], [self.z, self.z], c="k", alpha=0.4)
        ax.plot([self.x, self.x], [self.y-self.h,
                self.y+self.h], [self.z, self.z], c="k", alpha=0.4)

    def draw_(self, ax):
        ax.plot(*self._coordsg.T, c="red")


class Quadrant3DSpherical(Quadrant3D):
    def __init__(self, p: tuple, r: tuple):
        dim = [r, r, r]
        self.r = r
        Quadrant3D.__init__(self, p, dim)

    def contains(self, e: Brick) -> bool:
        xself = self.T(self.center.T)[0].flatten()
        x = e.T(e.center.T)[0].flatten()

        return (sum((xself-x)**2) <= self.r**2)


class Geometree():

    def __init__(self, boundary, n: int = 1, depth: int = 1) -> None:
        self.boundary = boundary
        self.points = []
        self.n = n
        self.divided = False
        self.children = []
        self.depth = depth

    def draw_points(self, ax):

        if not self.divided:
            for e in self.points:
                center = e.T(e.center.T)[0].flatten()
                plt.plot(*center, "o", c="black", alpha=0.5)
        for c in self.children:
            c.draw_points(ax)

    def draw(self, ax):

        if self.divided:
            self.boundary.draw(ax)
        for c in self.children:
            c.draw(ax)

    def contains(self, p: tuple) -> bool:
        return self.boundary.contains(p)

    def subdivide(self) -> None:
        self.divided = True
        self.children = []

        divs = self.boundary.subdivide()
        for d in divs:
            self.children.append(Geometree(d, self.n, self.depth+1))

    def add_point(self, p: tuple) -> bool:
        if not self.contains(p):
            return False
        if len(self.points) < self.n and not self.divided:
            self.points.append(p)
            return True
        if not self.divided:
            self.subdivide()
            for p2 in self.points[::-1]:
                for sq in self.children:
                    if sq.add_point(p2):
                        self.points.pop()
                        break
        for sq in self.children:
            if sq.add_point(p):
                return True
        raise Exception("This should never happen")

    def query_range(self, quadrant) -> bool:
        result = []
        if not self.boundary.intesects_quadrant(quadrant):
            return result
        for p in self.points:
            if quadrant.contains(p):
                result.append(p)
        if not self.divided:
            return result
        for sq in self.children:
            result += sq.query_range(quadrant)
        return result

    def query_range_point_radius(self, p, r, plot=False):
        q = Quadrant3DSpherical(p, r)
        selected = self.query_range(q)
        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(projection="3d")
            self.draw_points(ax)
            q.draw_(ax)
            for e in selected:
                center = e.T(e.center.T)[0].flatten()
                plt.plot(*center, "o", c="yellow")
            plt.show()
        return selected