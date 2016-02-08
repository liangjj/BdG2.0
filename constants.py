''' Load a variety of relevant physical parameters.

All quantities are in atomic units, such that
    m_e = 1
    e = 1
    hbar = 1
    1/4\pi\epsilon = 1
'''

import numpy as np

hbar = 1.0
m_e = 1.0
h22m = hbar**2 / (2*m_e)
pi = np.pi
eV = 1/27.21138505
eV_Ha = eV
nm = 18.89726124565

kB_eV = 8.6173324e-5
kB = kB_eV * eV_Ha 
