from scitbx.matrix import sqr,col
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from __future__ import division
'''
/          \     /    \
| 1. 2. 1. |     | 2. |
| 2. 2. 2. | X = | 0. |
| 1. 2. 3. |     | 1. |
\         /      \   /
'''
A = sqr( [1.,2.,1.,2.,2.,2.,1.,2.,3.] ) #Dense representation of above CSR 3x3 matrix
b=col([2.,0.,1.])
s_mat_result =  ((A.inverse()*b)).as_flex_double_matrix()
print "scitbx.matrix solution:="
print list(s_mat_result)

import boost.python
ext = boost.python.import_ext("scitbx_examples_strumpack_solver_ext")

n_rows = 3
n_cols = 3
nnz = 8

A_row_offset = flex.int([0,3,6,9])
A_col_offset = flex.int([0,1,2,0,1,2,0,1,2])
A_values = flex.double([1.,2.,1.,2.,2.,2.,1.,2.,3.])

print "A_row_offset=",list(A_row_offset)
print "A_col_offset=",list(A_col_offset)
print "A_values=",list(A_values)

b=flex.double([2.,0.,1.])
res = ext.sparse_solver(n_rows, n_cols, A_row_offset, A_col_offset, A_values, b)
strum_result = res.x
print "strumpack solution:="
print list(strum_result)

assert approx_equal(s_mat_result, strum_result)
print "scitbx.matrix solution and strumpack solution agree"
