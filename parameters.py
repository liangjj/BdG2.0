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
        self.N0_3D = (1 + self.lam) * np.sqrt(self.EF_3D / h22m) / (
            h22m * 4.0 * pi**2.0)

        #TODO: fix the dimensionality here 
        self.N0 = 1/(4*pi*h22m) / self.geometry.Lz
        self.Mu = self.chemical_potential(self.carrier_density, 
                self.geometry.Lz)
        self.nu = self.get_imax(self.Mu, self.geometry.Lz)
        self.kmax = self.get_kmax(self.material.hw_debye,
                self.Mu, self.geometry.Lz)
