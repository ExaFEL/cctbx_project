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
 * Modified by Lee J. O'Riordan to adapt for use with Boost.Python
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
  //template<typename numType, typename intType> 
  struct my_solver{

    /*my_solver(intType n)//int argc, char* argv[])
    {
      StrumpackSparseSolver<numType,intType> spss;
      spss.options().set_mc64job(0);
      spss.options().set_reordering_method(ReorderingStrategy::GEOMETRIC);
      //spss.options().set_from_command_line(argc, argv);

      intType N = n * n;
      intType nnz = 5 * N - 4 * n;
      CSRMatrix<numType,intType> A(N, nnz);
      intType* col_ptr = A.get_ptr();
      intType* row_ind = A.get_ind();
      numType* val = A.get_val();

      nnz = 0;
      col_ptr[0] = 0;
      for (intType row=0; row<n; row++) {
        for (intType col=0; col<n; col++) {
          intType ind = col+n*row;
          val[nnz] = 4.0;
          row_ind[nnz++] = ind;
          if (col > 0)  { val[nnz] = -1.0; row_ind[nnz++] = ind-1; } // left
          if (col < n-1){ val[nnz] = -1.0; row_ind[nnz++] = ind+1; } // right
          if (row > 0)  { val[nnz] = -1.0; row_ind[nnz++] = ind-n; } // up
          if (row < n-1){ val[nnz] = -1.0; row_ind[nnz++] = ind+n; } // down
          col_ptr[ind+1] = nnz;
        }
      }
      //A.set_symmetric_sparsity();

      std::vector<numType> b(N, numType(1.)), x(N, numType(0.));

      spss.set_csr_matrix(N, col_ptr, row_ind, val, true);
      spss.reorder(n, n);
      // spss.factor();   // not really necessary, called if needed by solve
      spss.solve(b.data(), x.data());

      // just a check, system is already solved, so solving again
      // with the solution as initial guess should stop immediately
      //spss.solve(b.data(), x.data(), true);
      b_res.swap(b);
      x_res.swap(x);

      std::cout << "# COMPONENTWISE SCALED RESIDUAL = " << A.max_scaled_residual(x_res.data(), b_res.data()) << std::endl;
    } */

    /*
    CSRMatrix<numType,intType> af_to_csr( scitbx::af::shared<intType> A_row_offset, 
                                          scitbx::af::shared<intType> A_col_offset, 
                                          scitbx::af::shared<numType> A_values)
    {
            CSRMatrix<numType,intType> A(N, nnz);

    }*/

    my_solver( const int n_rows, const int n_cols,
               const scitbx::af::shared<intType> A_row_offset, 
               const scitbx::af::shared<intType> A_col_offset, 
               const scitbx::af::shared<numType> A_values,
               const scitbx::af::shared<numType> b){//int argc, char* argv[])
      StrumpackSparseSolver<numType,intType> spss;
      spss.options().set_mc64job(0);
      //spss.options().set_reordering_method(ReorderingStrategy::GEOMETRIC);
      //spss.options().set_from_command_line(argc, argv);

      std::size_t nnz = A_values.size();
      std::size_t gridSize = n_rows * n_cols;// (A_row_offset.size()-1) * (scitbx::af::max(A_col_offset)+1); //(Length of rows -1) * (the largest index+1) gives the full grid

      //The number of row entries-1 will be the number of different equations, thus the size of the b-vector and x-vector
      std::cout << "Got here 1 \n";
      scitbx::af::shared<numType> x(3,0.0);//(A_row_offset.size()-1, numType(0.));
      std::cout << "Got here 2 \n";

      const intType* r_ptr = A_row_offset.begin();
      const intType* c_ptr = A_col_offset.begin();
      const numType* val_ptr = A_values.begin();

      spss.set_csr_matrix(3, r_ptr, c_ptr, val_ptr, false);
      std::cout << "Got here 3 \n";
      SCITBX_EXAMINE (*( r_ptr ));
      SCITBX_EXAMINE (*( c_ptr ));
      SCITBX_EXAMINE (*( val_ptr ));

      //CSRMatrix<numType,intType> AA(3, nnz);


      //spss.reorder();//( n_rows, n_cols);//(A_row_offset.size()-1), ( scitbx::af::max(A_col_offset)+1 ) );
      std::cout << "Got here 4 \n";

      x_res.push_back(0.);x_res.push_back(0.);x_res.push_back(0.);
      std::cout << x_res.size() << "\n";

      spss.solve(b.begin(), x_res.begin());
      std::cout << "Got here 5 \n";

      // just a check, system is already solved, so solving again
      // with the solution as initial guess should stop immediately
      //spss.solve(b.data(), x.data(), true);
      //x_res.swap(x);

      //std::cout << "# COMPONENTWISE SCALED RESIDUAL = " << A.max_scaled_residual(x_res.data(), b_res.data()) << std::endl;
    }

    scitbx::af::shared<numType> x_res;
    scitbx::af::shared<numType> b_res;
  }; //end of struct sparse_solver

namespace boost_python{
namespace {

  void export_strumpack_solver()
  {
    typedef return_value_policy<return_by_value> rbv;
    //templated constructor for int and size_t flex arrays of id
    class_<sparse_solver::my_solver> ("sparse_solver", no_init)
      .def( init< int, int, 
                 scitbx::af::shared<intType>,
                 scitbx::af::shared<intType>,
                 scitbx::af::shared<numType>,
                 scitbx::af::shared<numType> >
      ((  
        boost::python::arg("n_rows"), 
        boost::python::arg("n_cols"), 
        boost::python::arg("A_row_offset"), 
        boost::python::arg("A_col_offset"), 
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

