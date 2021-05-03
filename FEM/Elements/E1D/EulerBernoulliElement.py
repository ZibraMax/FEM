from .QuadraticElement import QuadraticElement, Element1D


class EulerBernoulliElement(QuadraticElement):
    """docstring for EulerBernoulliElement
    """

    def __init__(self, coords, gdl):
        QuadraticElement.__init__(self, coords, gdl, n=3)
