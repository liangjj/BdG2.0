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


class Structure:
    ''' Class that holds the dimensional information

    Sets the dimensions of the system. Can be an array of length 0, 1, 2 or 3,
    depending on the dimensionality of the system.
    This is important to charactarize k-space and k-space integrations.

    TODO: Maybe I should make this responsible for making a "k-space measure"?
    '''

    def __init__(self, dims):
        self.dims = dims
        self.ndim = 3-dims.size
        self.cont_dims = [0,1,2][:self.ndim]
        self.discr_dims = [0,1,2][self.ndim:]

class Parameters:
    ''' Class that holds all derived system parameters.

    From material parameters and dimensionality held in Material and Structure
    objects, we can calculate a whole bunch of derived properties, such as the
    chemical potential, Fermi momentum (max k-vector), maximal band index,
    theoretical Fermi-energies in different dimensionalities, DOS at Fermi
    level, ...
    '''
    # TODO: Make parameters depend on model!
    def __init__(self, material, structure):
        self.carrier_density = material.n / (nm**3)
        self.hw_debye = material.hw_debye
        self.lam = material.lam
        self.V0 = material.V0

        # Theoretical bulk values
        self.EF_3D = h22m * (3.0 * pi**2.0 * self.carrier_density)**(2.0/3.0)
        self.N0_3D = (1 + self.lam) * np.sqrt(self.EF_3D / h22m) / (
            h22m * 4.0 * pi**2.0)

        # Dimension-dependent stuff
        if structure.ndim == 0:
            # Quantum dot
            self.Lx = structure.dims[0]
            self.Ly = structure.dims[1]
            self.Lz = structure.dims[2]
            self.Mu = 1.0
            self.N0 = 1.0
            self.kmax = 1
            self.nu = 1
        elif structure.ndim == 1:
            # Quantum wire
            self.Lx = structure.dims[0]
            self.Ly = structure.dims[1]
            self.Mu = 1.0
            self.N0 = 1.0
            self.kmax = 1
            self.nu = 1
        elif structure.ndim == 2:
            # Film
            self.Lz = structure.dims[0]
            # self.EF = h22m * 2.0 * pi * self.carrier_density
            self.N0 = 1/(4 * pi * h22m)/self.Lz
            self.Mu = self.chemical_potential(self.carrier_density, self.Lz)
            self.nu = self.get_imax(self.Mu, self.Lz)
            self.imax = self.nu # Just because I tend to use both...
            self.kmax = self.get_kmax(self.hw_debye, self.Mu, self.Lz)
        elif structure.ndim == 3:
            self.Mu = self.EF_3D
            self.N0 = self.N0_3D
            self.nu = 0.1
            self.kmax = 10.0

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
                                          kmax = self.kmax)))


class Integrator:
    """ Class that has knowledge of k-space and can do integrations etc...

    From a Structure, we can deduce the k-space structure. We can use an
    Integrator object to do summations over (regions of) k-space.  """
    N = 50 # Number of samples for continuous integration

    def __init__(self, struct):
        if struct.ndim == 0:
            self.Lx, self.Ly, self.Lz = struct.dims[:]
            self.dx, self.dy, self.dz = 1/struct.dims[:]
        elif struct.ndim == 1:
            self.Lx = 0
            self.Ly, self.Lz = struct.dims[:]
            self.dx = 1/self.N
            self.dy, self.dz = 1/struct.dims[:]
        elif struct.ndim == 2:
            self.Lx = self.Ly = 0
            self.Lz = struct.dims[0]
            self.dx = self.dy = 1/self.N
            self.dz = 1/struct.dims[0]
        elif struct.ndim == 3:
            self.Lx = self.Ly = self.Lz = 0
            self.dx = self.dy = self.dz = 1/self.N

    def create_arrays(self, lims, axes):
        print(axes)
        if 0 in axes:
            idx = axes.index(0)
            X = np.arange(0, lims[idx], lims[idx]*self.dx)[:, None, None]
        else:
            X = np.array([0])[:,None, None] 

        if 1 in axes:
            idy = axes.index(1)
            Y = np.arange(0, lims[idy], lims[idy]*self.dy)[None, :, None]
        else:
            Y = np.array([0])[None, :, None] 

        if 2 in axes:
            idz = axes.index(2)
            Z = np.arange(0, lims[idz], lims[idz]*self.dz)[None, None, :]
        else:
            Z = np.array([0])[None, None, :] 

        return X, Y, Z
    
    def integrate(self, f, lims, axes = [0,1,2]):
        """ Integrate over a subset of dimensions (hence selective).

        Integrate the passed function (assumed to be R3 -> R) over the subset 
        of axes passed in 'axes', with limits going from 0 to 'klims'.
        The trick is that skipped dimensions are represented by a 1-element
        array [[[0]]]. This is necessary, because we assume the passed function
        f to accept 3 arrays, so we can't skip any...
        """
        X, Y, Z = self.create_arrays(lims, axes)
        fvals = f(X, Y, Z)
        return self.dx*self.dy*self.dz*np.sum(fvals, keepdims=True)
