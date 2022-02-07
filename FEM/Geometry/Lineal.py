"""Generate a 1D Domain"""


from .Geometry import Geometry, np


class Lineal(Geometry):

    """Generate a evenly spaced elements domain

    Args:
        lenght (float): Domain lenght
        n (int): Number of elements
        o (int): Element order, can be 1 or 2
        nvn (int, optional): Number of variables per node. Defaults to 1.
    """

    def __init__(self, lenght: float, n: int, o: int, nvn: int = 1) -> None:
        """Generate a evenly spaced elements domain

        Args:
            lenght (float): Domain lenght
            n (int): Number of elements
            o (int): Element order, can be 1 or 2
            nvn (int, optional): Number of variables per node. Defaults to 1.
        """

        self.lenght = lenght
        dictionary = []
        gdls = []
        he = self.lenght / (n)
        for i in range(0, n):
            xa = i * he
            if o == 1:
                gdls += [xa]
                dictionary += [[i, i+1]]
            elif o == 2:
                gdls += [xa, xa+he/2]
                dictionary += [[i*o, i*o+1, i*o+2]]
            else:
                gdls += [xa, xa+he/3, xa+2*he/3]
                dictionary += [[i*o, i*o+1, i*o+2, i*o+3]]
        gdls += [self.lenght]
        if o == 1:
            tipo = 'L1V'
        elif o == 2:
            tipo = 'L2V'
        else:
            tipo = 'L3V'
        types = [tipo]*len(dictionary)
        gdls = np.array(gdls).reshape([len(gdls), 1])
        Geometry.__init__(self, dictionary, gdls, types, nvn=nvn, segments=[])
