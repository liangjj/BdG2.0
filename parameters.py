""" Auxiliary classes for sets of parameters.
"""
from collections import namedtuple

Material = namedtuple("Material",
                      ["carrier_density", "debye_energy", "lambda", "V0"])

Shape = namedtuple("Shape", ["xsize", "ysize", "zsize"])


class Parameters:
    """ Wrapper class for a whole set of parameters.

    Extra parameters are calculated from the Material and Shape properties
    (stuff like Fermi velocity, Debye temperature, chemical potential, etc...
    """
    def __init__(self, material, shape):
        self.material = material
        self.shape = shape
