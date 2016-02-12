from parameters import Material, Shape, Parameters
from system import System, Hamiltonian
import numpy as np
from functools import partial
import matplotlib.pyplot as plt

sn = Material(148.0, 0.000618, 0.66666, 0.1)
shape = Shape(200, 200, 10)
params = Parameters(sn, shape)

sys = System(sn, shape)

ham = Hamiltonian(shape)
a = np.array([[1, 2], [3, 4]])
i = np.arange(0, 5)[None, :]
emin = 0.0
emax = 1.0
e = np.linspace(emin, emax, 1000)[:, None]

dosi = partial(ham.dos, i)


def F(T, y):
    return 1

I = sys.tcsolver.energy_integration(F, emin, emax, dosi)
dostot = np.sum(dosi(e), axis=1)

M = sys.tcsolver.scf_matrix(ham.dos, F)
print(M)
print(M.shape)
# plt.plot(e, dostot)
# plt.show()
