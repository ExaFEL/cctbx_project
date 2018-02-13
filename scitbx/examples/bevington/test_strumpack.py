from scitbx.examples.bevington import sparse_solver as ss
from cctbx.array_family import flex
from __future__ import division

n_rows = 3
n_cols = 3
nnz = 8 #Cross pattern of 1's, Identity with periodic boundary

A_row_offset = flex.int([0,3,5,8])
A_col_offset = flex.int([0,1,2,1,2,0,1,2])
A_values = flex.double([1.,2.,1.,1.,2.,3.,2.,1.])

b=flex.double([2.,0.,1.])
res = ss(n_rows, n_cols, A_row_offset, A_col_offset, A_values, b)
print list(res.x)
