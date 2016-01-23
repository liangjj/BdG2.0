""" Module that contains core functionality.

Classes contained:
    - Hamiltonian: Initialized as a free particle Hamiltonian, which can return
    a 'Spectrum' or a 'Dos'. 'Corrections' can be applied to it, resulting in a
    modified 'Spectrum' and 'Dos'.

"""

import numpy as np
from constants import *
from parameters import Material, Structure, Parameters, Integrator


class Hamiltonian:
    ''' A Hamiltonian governs the spectrum and dos of a system.

    Fields: - Spectrum (function)
            - Dos (function)
    '''

    def __init__(self):
        self.spectrum = self.free_spectrum

    def free_spectrum(self, kx, ky, kz):
        """ Return the quadratic spectrum

        We implicitely assume any of the passed k's can be a discrete set
        of bands, though provided in the form (i+1)²pi²/L².
        """
        return h22m * (kx**2 + ky**2 + kz**2)


class System:
    ''' A system is responsible for calculating the superconducting properties.

    A system is defined by a Hamiltonian and a set of system parameters, and can
    calculate the superconducting properties from this. (Tc and Delta, and maybe
    more).
    '''

    def __init__(self, ham, par, struct):
        self.H = ham
        self.par = par
        self.struct = struct

    def overlaps(self):
        """ The interaction matrix elements.

        These are given by the wavefunction overlaps of the simple free-particle
        wavefunctions.
        At the moment, this is only implemented for 2D!
        """

        size = self.par.nu
        Lz = self.struct.dims[-1]
        return (np.ones((size, size)) + 0.5*np.identity(size))/Lz

    def thermal_weight_LO_k(self, kx, ky, kz, T):
        ''' The thermal weight appearing in the band gap equation.

        A factor proportional to tanh(beta E) weights the "propagator"
        1/E. (Where E = sqrt(ksi^2 + Delta^2). To leading order (hence LO), the
        Delta's drop out, and we get a factor tanh(beta ksi)/ksi.
        
        This is simply a function from R^3 -> R, but can be vectorized (in light
        of the composition shenanigans I intend to perform).
        '''
        ksi = self.H.spectrum(kx, ky, kz)
        return np.tanh(ksi/(2*kB*T))/ksi 

    def thermal_weight_LO(self, ksi, T):
        ''' The thermal weight appearing in the band gap equation.

        A factor proportional to tanh(beta E) weights the "propagator"
        1/E. (Where E = sqrt(ksi^2 + Delta^2). To leading order (hence LO), the
        Delta's drop out, and we get a factor tanh(beta ksi)/ksi.
        
        This is simply a function from R -> R, but can be vectorized (in light
        of the composition shenanigans I intend to perform).
        '''

        return np.tanh(ksi/(2*kB*T))/ksi 

    def thermal_weight_LO_atT(self, T):
        def thermal_weight(kx, ky, kz):
            return self.thermal_weight_LO_k(kx, ky, kz, T)
        return thermal_weight

    def getTc_det(self):
        ### Overlap matrix
        Phi = self.overlaps()

        ### Thermal factor
        integrator = Integrator(self.struct)
        T = 0.1
        tw = self.thermal_weight_LO_atT(T)

        axes = self.struct.cont_dims
        lims = self.par.kmax*np.ones(len(axes))
        A = integrator.integrate(tw, lims, axes)
        
        return A
