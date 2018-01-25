/*
 * STRUMPACK -- STRUctured Matrices PACKage, Copyright (c) 2014, The Regents of
 * the University of California, through Lawrence Berkeley National Laboratory
 * (subject to receipt of any required approvals from the U.S. Dept. of Energy).
 * All rights reserved.
 *
 * If you have questions about your rights to use or distribute this software,
 * please contact Berkeley Lab's Technology Transfer Department at TTD@lbl.gov.
 *
 * NOTICE. This software is owned by the U.S. Department of Energy. As such, the
 * U.S. Government has been granted for itself and others acting on its behalf a
 * paid-up, nonexclusive, irrevocable, worldwide license in the Software to
 * reproduce, prepare derivative works, and perform publicly and display publicly.
 * Beginning five (5) years after the date permission to assert copyright is
 * obtained from the U.S. Department of Energy, and subject to any subsequent five
 * (5) year renewals, the U.S. Government is granted for itself and others acting
 * on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
 * Software to reproduce, prepare derivative works, distribute copies to the
 * public, perform publicly and display publicly, and to permit others to do so.
 *
 * Developers: Pieter Ghysels, Francois-Henry Rouet, Xiaoye S. Li.
 *             (Lawrence Berkeley National Lab, Computational Research Division).
 *
 * Modified by Lee J. O'Riordan to adapt for use with cctbx via Boost.Python 
 */
#include <iostream>
#include "StrumpackSparseSolver.hpp"
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>

using namespace boost::python;
using namespace strumpack;
typedef double numType;
typedef int intType;

namespace sparse_solver {
  struct my_solver{
    /**
      Constructor for the STRUMPACK Ax=b sparse solver my_solver.
      This object will be callable from python upon importing the appropriate extension module using
      'from scitbx.examples.bevington import sparse_solver as ss'
      Result vector x_res will be initialised by constructor to the same size as n_rows to contain solution.

      @param[in]  n_rows Row count of the A matrix.
      @param[in]  n_cols Column count of the A matrix.
      @param[in]  A_row_offset Row offset indices for CSR sparse format matrix.
      @param[in]  A_col_idx Column indices for CSR sparse format matrix.
      @param[in]  A_values A matrix values in CSR sparse format matrix.
      @param[in]  b Right hand side vector of linear system.
      @return Description of returned value.
    */
    my_solver(const int n_rows, const int n_cols,
              const scitbx::af::shared<intType> A_row_offset, 
              const scitbx::af::shared<intType> A_col_idx, 
              const scitbx::af::shared<numType> A_values,
              const scitbx::af::shared<numType> b) : x_res(n_rows, 0.){//int argc, char* argv[])
      StrumpackSparseSolver<numType,intType> spss;
      spss.options().set_mc64job(0);
      //spss.options().set_reordering_method(ReorderingStrategy::GEOMETRIC);
      //spss.options().set_from_command_line(argc, argv);

      std::size_t nnz = A_values.size();
      std::size_t gridSize = n_rows * n_cols;

      //Taking raw pointer to flex data types
      const intType* r_ptr = A_row_offset.begin();
      const intType* c_ptr = A_col_idx.begin();
      const numType* val_ptr = A_values.begin();

      //Set the sparse matrix values using the CSR formatted data
      spss.set_csr_matrix(n_rows, r_ptr, c_ptr, val_ptr, false);
      spss.reorder( n_rows, n_cols );

      //Solve Ax=b
      spss.solve(b.begin(), x_res.begin());
    }
    scitbx::af::shared<numType> x_res; //Resulting x will be stored here
  }; //end of struct sparse_solver

namespace boost_python{
namespace {
  void export_strumpack_solver() {
    typedef return_value_policy<return_by_value> rbv;
    class_<sparse_solver::my_solver> ("sparse_solver", no_init)
      .def( init< int, int, 
                 scitbx::af::shared<intType>,
                 scitbx::af::shared<intType>,
                 scitbx::af::shared<numType>,
                 scitbx::af::shared<numType> 
                >
      ((  
        boost::python::arg("n_rows"), 
        boost::python::arg("n_cols"), 
        boost::python::arg("A_row_offset"), 
        boost::python::arg("A_col_idx"), 
        boost::python::arg("A_values"), 
        boost::python::arg("b")
      ))
      )
      .add_property("x",make_getter(&sparse_solver::my_solver::x_res, rbv()))
    ;
  }
}//anon
}//boost_python
}//sparse_solver
BOOST_PYTHON_MODULE(scitbx_examples_strumpack_solver_ext)
{
   sparse_solver::boost_python::export_strumpack_solver();
}

