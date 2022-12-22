"""
DEPRECATED

This class was my introduction to OctTrees. I implemented all geometry with this class. Then, I realized SciPy has a KDTree implementation, faster and better. IDK what to tink about it
"""

import numpy as np
import matplotlib.pyplot as plt


class Quadrant3D():

    def __init__(self, p: tuple, dim: tuple) -> None:
        x, y, z = p
        w, h, d = dim
        self.x, self.y, self.z = x, y, z
        self.w, self.h, self.d = w, h, d
        self.coords = np.array([
            [x-w, y-h, z-d],
            [x+w, y-h, z-d],
            [x+w, y+h, z-d],
            [x-w, y+h, z-d],
            [x-w, y-h, z+d],
            [x+w, y-h, z+d],
            [x+w, y+h, z+d],
            [x-w, y+h, z+d]])
        self.maximos_self = np.max(self.coords, axis=0)
        self.minimos_self = np.min(self.coords, axis=0)
        self._xcenter = np.array([x, y, z])

    def contains(self, e) -> bool:
        x = e._xcenter

        superior = (self.maximos_self-x) >= 0
        inferior = (x-self.minimos_self) >= 0

        return superior.all() and inferior.all()

    def boxes_disjoint(self, e):

        maxx1, maxy1, maxz1 = self.maximos_self
        minx1, miny1, minz1 = self.minimos_self

        maxx2, maxy2, maxz2 = e.maximos_self
        minx2, miny2, minz2 = e.minimos_self

        return (maxx2 <= minx1 or maxx1 <= minx2
                or maxy2 <= miny1 or maxy1 <= miny2
                or maxz2 <= minz1 or maxz1 <= minz2)

    def intesects_quadrant(self, e) -> bool:
        return not self.boxes_disjoint(e)

    def subdivide(self) -> list:
        divs = []
        nw = self.w/2
        nh = self.h/2
        nd = self.d/2
        x, y, z = self.x, self.y, self.z

        divs.append(Quadrant3D((x+nw, y+nh, z+nd), (nw, nh, nd)))
        divs.append(Quadrant3D((x+nw, y+nh, z-nd), (nw, nh, nd)))
        divs.append(Quadrant3D((x+nw, y-nh, z+nd), (nw, nh, nd)))
        divs.append(Quadrant3D((x+nw, y-nh, z-nd), (nw, nh, nd)))
        divs.append(Quadrant3D((x-nw, y+nh, z+nd), (nw, nh, nd)))
        divs.append(Quadrant3D((x-nw, y+nh, z-nd), (nw, nh, nd)))
        divs.append(Quadrant3D((x-nw, y-nh, z+nd), (nw, nh, nd)))
        divs.append(Quadrant3D((x-nw, y-nh, z-nd), (nw, nh, nd)))
        return divs

    def draw(self, ax):
        ax.plot([self.x, self.x], [self.y, self.y],
                [self.z-self.d, self.z+self.d], c="k", alpha=0.4)
        ax.plot([self.x-self.w, self.x+self.w],
                [self.y, self.y], [self.z, self.z], c="k", alpha=0.4)
        ax.plot([self.x, self.x], [self.y-self.h,
                self.y+self.h], [self.z, self.z], c="k", alpha=0.4)

    def draw_(self, ax):
        ax.plot(*self.coords.T, c="red")


class Quadrant3DSpherical(Quadrant3D):
    def __init__(self, p: tuple, r: tuple):
        dim = [r, r, r]
        self.r = r
        Quadrant3D.__init__(self, p, dim)

    def contains(self, e) -> bool:

        return (sum((self._xcenter-e._xcenter)**2) <= self.r**2)


class Geometree():
    min_search_size = -1

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
                plt.plot(*e._xcenter, "o", c="black", alpha=0.5)
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
        dist = p.coords-p._xcenter
        min_search_size = max(np.sum(dist**2, axis=1)**0.5)
        self.min_search_size = max(min_search_size, self.min_search_size)
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

    def query_range(self, quadrant, plot=False, ax=None) -> bool:
        result = []
        if not self.boundary.intesects_quadrant(quadrant):
            return result

        for p in self.points:
            if plot:
                ax.plot(*p._xcenter, "o", c="green", alpha=1)
            if quadrant.contains(p):
                result.append(p)
        if not self.divided:
            return result
        for sq in self.children:
            result += sq.query_range(quadrant, plot=plot, ax=ax)
        return result

    def query_range_point_radius(self, p, r=None, plot=False, ax=None):
        if r == None:
            r = 2*self.min_search_size
        q = Quadrant3DSpherical(p, r)
        selected = self.query_range(q, plot, ax)
        return selected

    def graph_query_range(self, p, r):
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
        self.draw_points(ax)
        result = self.query_range_point_radius(p, r, True, ax)
        for p in result:
            ax.plot(*p._xcenter, "o", c="yellow", alpha=1)
        plt.show()

    def query_first_point_set(self):
        if self.divided:
            for ch in self.children:
                if ch.children or ch.points:
                    return ch.query_first_point_set()
        else:
            return self.points
        raise Exception("This should not happen")
