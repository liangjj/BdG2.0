import constants as c
import numpy as np


class TcSolver:
    ''' Class in charge of solving for Tc.
    '''
    def __init__(self, parameters, hamiltonian):
        self.parameters = parameters
        self.phi = self.overlaps()
        self.F = self.weight
        self.H = hamiltonian

    def weight(self, T, ksi):
        return np.tanh(ksi/c.kB/T)/ksi

    def overlaps(self):
        N = self.parameters.nu
        return (np.ones((N, N)) + np.identity(N))/self.parameters.shape.Lz

    def energy_integration(self, fun, emin, emax, dos):
        ne = 1000
        e = np.linspace(emin, emax, ne)
        return N

