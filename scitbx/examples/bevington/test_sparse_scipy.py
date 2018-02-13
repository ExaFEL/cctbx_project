#Small-scale test script for STRUMPACK solution, compared with scipy
from __future__ import division
import numpy as np
import scipy.sparse
import scipy.sparse.linalg

A_row_offset = np.array([0,3,5,8])
A_col_offset = np.array([0,1,2,1,2,0,1,2])
A_values = np.array([1.,2.,1.,1.,2.,3.,2.,1.])
A = scipy.sparse.csr_matrix((A_values, A_col_offset, A_row_offset))
b = np.array([2,0,1])

x = scipy.sparse.linalg.spsolve(A, b)
print x
