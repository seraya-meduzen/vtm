import numpy as np


mtx = np.matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0.01]])

mtx = np.conjugate(mtx)


u, s, vt = np.linalg.svd(mtx, full_matrices=False, hermitian=False)

print(u)
