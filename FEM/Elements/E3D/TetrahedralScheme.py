"""Define the tetrahedral scheme used by tetrahedral elements
"""


import numpy as np


class TetrahedralScheme():
    """Generate a tetrahedral integration scheme

    Args:
        n (int): Number of gauss points
    """

    def __init__(self, n: int, **kargs) -> None:
        """Generate a tetrahedral integration scheme

        Args:
            n (int): Number of gauss points
        """

        self.Z = np.array([[0.01583591, 0.328054697, 0.328054697],
                           [0.328054697, 0.01583591, 0.328054697],
                           [0.328054697, 0.328054697, 0.01583591],
                           [0.328054697, 0.328054697, 0.328054697],
                           [0.679143178, 0.106952274, 0.106952274],
                           [0.106952274, 0.679143178, 0.106952274],
                           [0.106952274, 0.106952274, 0.679143178],
                           [0.106952274, 0.106952274, 0.106952274]])
        self.W = np.array([0.023087995, 0.023087995, 0.023087995, 0.023087995,
                          0.018578672, 0.018578672, 0.018578672, 0.018578672])
        self.center = np.array([[1.0, 1.0, 1.0]])/3.0
        self.domain = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [
                               0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] + self.Z.tolist())
