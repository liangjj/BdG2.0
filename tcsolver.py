import constants as c
import numpy as np
from functools import partial


class TcSolver:
    ''' Class in charge of solving for Tc.
    '''
    def __init__(self, parameters, hamiltonian):
        self.parameters = parameters
        self.F = self.weight
        self.H = hamiltonian

    def weight(self, T, ksi):
        return np.tanh(ksi/c.kB/T)/ksi

    def overlaps(self):
        N = self.parameters.nu
        return (np.ones((N, N)) + np.identity(N))/self.parameters.shape.Lz

    def energy_integration(self, fun, emin, emax, dos):
        """ Takes a partially applied DOS.

        If we partially apply the range of bands to the DOS, we get a DOS that
        returns an array of subband densities.
        Use this DOS (so, the DOS no longer needs to take an i argument, it's
        already included.
        """
        ne = 1000  # Makes more sense to define a de (delta e) instead?
        e = np.linspace(emin, emax, ne)[:, None]
        F = fun(1, e)*dos(e)
        I = np.trapz(F, e, axis=0)
        return I

    def determinant(self, T):
        ''' Calculate the determinant of the matrix obtained by linearizing the
        the SCF equation.
        '''

        phi = self.overlaps()
        emin = -self.parameters.material.debye_energy
        emax = self.parameters.material.debye_energy
        i = np.arange(0, self.parameters.nu)[None, :]

        F = partial(self.weight, T)
        dosi = partial(self.ham.dos, i)
        M = self.energy_integration(F, emin, emax, dosi)
        return np.det(M*phi)

    def calc_Tc(self):
        return np.brentq(self.determinant, 0, 10)
