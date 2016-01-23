''' Module containing auxiliary class definitions.

Several useful wrapper classes are defined that take some of the leveraging:
    - Material: Holds information such as material, mass, electron-phonon
                coupling, etc...
    - Structure: Holds dimensional information: what are the x, y, and z
                 dimensions? This determines if the system is bulk, film, wire
                 or dot.
    - Parameters: From Material and Structure, we can calculate all the relevant
                  parameters of a corresponding system: chemical potential,
                  fermi momentum, etc...
'''

import numpy as np
import string as s
from constants import *
from textwrap import dedent


class Material:
    ''' Class that holds material information.

    Several properties define a material:
    - Carrier density
    - el-phon coupling
    - Debye energy
    - Electron-electron coupling V_0
    '''

    def __init__(self, n, hw_debye, lam, V0):
        self.n = n
        self.hw_debye = hw_debye
        self.lam = lam
        self.V0 = V0


class Geometry:
    ''' Class that holds the dimensional information
    '''

    def __init__(self, Lx, Ly, Lz):
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz


class Parameters:
    ''' Class that holds all derived system parameters.

    From material parameters and dimensionality held in Material and Structure
    objects, we can calculate a whole bunch of derived properties, such as the
    chemical potential, Fermi momentum (max k-vector), maximal band index,
    theoretical Fermi-energies in different dimensionalities, DOS at Fermi
    level, ...
    '''

    def __init__(self, material, geometry):
        self.material = material
        self.geometry = geometry

        self.carrier_density = material.n / (nm**3)
        # self.hw_debye = material.hw_debye
        # self.lam = material.lam
        # self.V0 = material.V0

        # Theoretical bulk values
        self.EF_3D = h22m * (3.0 * pi**2.0 * self.carrier_density)**(2.0/3.0)
        self.N0_3D = (1 + self.material.lam) * np.sqrt(self.EF_3D / h22m) / (
            h22m * 4.0 * pi**2.0)

        # TODO: fix the dimensionality here
        self.N0 = 1/(4*pi*h22m) / self.geometry.Lz
        self.Mu = self.chemical_potential(self.carrier_density,
                                          self.geometry.Lz)
        self.nu = self.get_imax(self.Mu, self.geometry.Lz)
        self.kmax = self.get_kmax(self.material.hw_debye,
                                  self.Mu, self.geometry.Lz)

    # Auxiliary routines
    def calculate_chempot(self, n, nu, L):
        ''' Calculate the 2D chemical potential.

        Calculate chemical potential from parameters according to [ref].
        '''
        mu = 2 * h22m * np.pi * L / nu * (n + np.pi / (6*L**3) * nu *
                                          (nu + 0.5) * (nu + 1))
        return mu

    def chemical_potential(self, n, L):
        ''' Find the 2D chemical potential self-consistently.
        '''
        tol = 0.0001
        mu_old = 0
        mu = 2 * h22m * np.pi * n
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
        return np.floor(L/pi * np.sqrt(1/h22m * mu))

    def get_kmax(self, hw_debye, mu, L):
        ''' Get the 2D maximal k-vector.

        '''
        # TODO: Generalize this to arbitrary dimensions.
        return np.sqrt((mu + hw_debye - pi**2/L**2)/h22m)

    def print_parameters(self):
        # TODO: Use string formatting after all...
        str = s.Template("""
        *** Fundamental parameters ********************************************
        Debye energy    = $debye Ha,          Carrier density = $n Bohr^-3
        BCS lambda      = $lam (dim'less),  Disorder V0     = $V0 (dim'less)
        *** Derived parameters ************************************************
        gN              = $gN,              D0      = $D0 Ha^-1
        Max band index  = $nu,                  g       = $g
        EF              = $mu Ha,     N(0)    = $N0 Ha^-1
        EF_3D           = $EF_3D,    N(0)_3D = $N0_3D
        v_F = p_F       = $vF,      Max k vector = $kmax
        """)

        return dedent(str.substitute(dict(debye=self.hw_debye,
                                          n=self.carrier_density,
                                          lam=self.lam,
                                          V0=self.V0,
                                          D0=0.0,
                                          gN=self.lam,
                                          nu=self.nu,
                                          g=self.lam/self.N0,
                                          mu=self.Mu,
                                          N0=self.N0,
                                          EF_3D=self.EF_3D,
                                          N0_3D=self.N0_3D,
                                          vF=0.0,
                                          kmax=self.kmax)))
