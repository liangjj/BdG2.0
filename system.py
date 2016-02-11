from parameters import Parameters
from constants import *
from tcsolver import TcSolver


class System:
    """The main interface with which to interact with the code.

        A system allows for calculating the superconducting gap, Tc, etc...
    """
    def __init__(self, material, shape):
        self.material = material
        self.shape = shape
        self.parameters = Parameters(material, shape)
        self.solver = TcSolver(self.parameters)


class Hamiltonian:  # Have Hamiltonian take a Shape object?
    """ A Hamiltonian basically stores information on the spectrum and DOS.
    """
    def __init__():
        pass

    def dos(i, e):  # Correct definition of DOS?
        return theta(e - h22m*(i+1)**2)


def theta(x):
    """Heaviside step function"""
    return 0.5*(1 + np.sign(x))
