import constants as c
import numpy as np


""" Auxiliary classes for sets of parameters.
"""
from collections import namedtuple

Material = namedtuple("Material",
                      ["carrier_density", "debye_energy", "lam", "V0"])

Shape = namedtuple("Shape", ["Lx", "Ly", "Lz"])


class Parameters:
    """ Wrapper class for a whole set of parameters.

    Extra parameters are calculated from the Material and Shape properties
    (stuff like Fermi velocity, Debye temperature, chemical potential, etc...
    """
    def __init__(self, material, shape):
        self.material = material
        self.shape = shape

        self.mu = self.chemical_potential(self.material.carrier_density,
                                          self.shape.Lz)
        self.nu = self.get_imax(self.mu, self.shape.Lz)
        self.kmax = self.get_kmax(self.material.debye_energy, self.mu,
                                  self.shape.Lz)

    def print_parameters(self):
        s = """
        ------------------------------------------------------------------------
        -- Parameters ----------------------------------------------------------
        ------------------------------------------------------------------------
        -- Material parameters -------------------------------------------------
        ------------------------------------------------------------------------
        -- n \t \t= {material.carrier_density}, \t\t lambda = {material.lam}
        -- debye_energy = {material.debye_energy}, \t\t V0 \t= {material.V0}
        ------------------------------------------------------------------------
        -- Shape parameters ----------------------------------------------------
        -- Lx = {shape.Lx},\t\t Ly = {shape.Ly},\t\t Lz = {shape.Lz}
        ------------------------------------------------------------------------
        -- Derived parameters --------------------------------------------------
        ------------------------------------------------------------------------
        -- Mu \t= {mu}, \t\t nu = {nu},
        -- kmax = {kmax}
        ------------------------------------------------------------------------
        """.format(material=self.material,
                   shape=self.shape,
                   mu=self.mu,
                   nu=self.nu,
                   kmax=self.kmax)
        print(s)

    # --------------------------------------------------------------------------
    # Auxiliary routines
    # --------------------------------------------------------------------------
    def calculate_chempot(self, n, nu, L):
        ''' Calculate the 2D chemical potential.

        Calculate chemical potential from parameters according to [ref].
        '''

        mu = 2 * c.h22m * np.pi * L / nu * (n + np.pi / (6*L**3) * nu *
                                            (nu + 0.5) * (nu + 1))
        return mu

    def chemical_potential(self, n, L):
        ''' Find the 2D chemical potential self-consistently.

        '''

        tol = 0.0001
        mu_old = 0
        mu = 2 * c.h22m * np.pi * n
        while (np.abs(mu - mu_old) > tol):
            mu_old = mu
            nu = self.get_imax(mu_old, L)
            mu = self.calculate_chempot(n, nu, L)

        return mu

    def get_imax(self, mu, L):
        '''Figure out the number of bands to consider.

        Input:  - mu:   Chemical potential
                - L:    Film thickness

        Output: imax: Highest occupied band index.
        '''
        # TODO: Generalize this to arbitrary dimensions!
        return np.floor(L/np.pi * np.sqrt(1/c.h22m * mu))

    def get_kmax(self, hw_debye, mu, L):
        ''' Get the 2D maximal k-vector.

        '''
        # TODO: Generalize this to arbitrary dimensions.
        return np.sqrt((mu + hw_debye - np.pi**2/L**2)/c.h22m)
