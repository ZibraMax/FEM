import numpy as np
from datetime import datetime
from scipy.linalg import eigh as largest_eigh
from scipy.sparse.linalg.eigen.arpack import eigsh as largest_eigsh
from scipy.sparse.linalg import eigs

np.set_printoptions(suppress=True)
np.random.seed(0)
N = 2500
k = 10
X = np.random.random((N, N)) - 0.5
X = np.dot(X, X.T)  # create a symmetric matrix

M = np.random.random((N, N)) - 0.5
M = np.dot(X, X.T)  # create a symmetric matrix


omh = np.linalg.solve(M, X)

# # Benchmark the dense routine
# start = datetime.now()
# print(start)
# evals_large, evecs_large = largest_eigh(omh, eigvals=(N-k, N-1))
# elapsed = (datetime.now() - start)
# print("eigh elapsed time: ", elapsed)
# print(evals_large)

# # Benchmark the sparse routine
# start = datetime.now()
# print(start)
# evals_large_sparse, evecs_large_sparse = largest_eigsh(omh, k, which='SM')
# elapsed = (datetime.now() - start)
# print("eigsh elapsed time: ", elapsed)
# print(evals_large)

# Benchmark the sparse routine
start = datetime.now()
print(start)
evals_large_sparse, evecs_large_sparse = eigs(X, k, M, which='SM')
elapsed = (datetime.now() - start)
print("eigs sparse KM elapsed time: ", elapsed)
print(evals_large_sparse)

# Benchmark the sparse routine
start = datetime.now()
print(start)
evals_large_sparse, evecs_large_sparse = eigs(omh, k, which='SM')
elapsed = (datetime.now() - start)
print("eigs sparse MK elapsed time: ", elapsed)
print(evals_large_sparse)
