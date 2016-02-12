from parameters import Parameters, Shape
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
        self.hamiltonian = Hamiltonian(self.shape)
        self.tcsolver = TcSolver(self.parameters, self.hamiltonian)


class Hamiltonian:  # Have Hamiltonian take a Shape object?
    """ A Hamiltonian basically stores information on the spectrum and DOS.
    """
    def __init__(self, shape):
        self.shape = shape

    def dos(self, i, e):  # Correct definition of DOS?
        ''' Density of states, *ABSOLUTE*, not relative to the fermi energy.
        '''
        return theta(e - h22m*np.pi**2*(i+1)**2/self.shape.Lz**2)


def theta(x):
    """Heaviside step function"""
    return 0.5*(1 + np.sign(x))
