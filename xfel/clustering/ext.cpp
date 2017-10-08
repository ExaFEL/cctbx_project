// Author:Nick Sauter
#include <cctbx/boost_python/flex_fwd.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>

namespace sx_clustering {
  void foo2(){
    std::cout<<"HELLO foo2"<<std::endl;
  }

  struct Rodriguez_Laio_clustering_2014 {

    Rodriguez_Laio_clustering_2014( ):rho_set(false){}

    Rodriguez_Laio_clustering_2014(double d_c):d_c(d_c),rho_set(false){}

    Rodriguez_Laio_clustering_2014(
      scitbx::af::flex_double distance_matrix,
      double d_c):
        distance_matrix( distance_matrix ),
        d_c(d_c),
        rho_set(false){
        NN = distance_matrix.accessor().focus()[0];
        SCITBX_ASSERT(NN == distance_matrix.accessor().focus()[1]);// square matrix
    }

    scitbx::af::shared<std::size_t>
    get_rho(){
      if (!rho_set) {
        rho = scitbx::af::shared<std::size_t>(NN);
        const double* Dij = distance_matrix.begin();
        for (int i=0; i < NN; ++i){
          for (int j=0; j < NN; ++j){
            if (Dij[i*NN + j] < d_c) { rho[i]+=1; }
          }
        }
        rho_set=true;
      }
      return rho;
    }

    scitbx::af::shared<double>
    get_delta(scitbx::af::shared<std::size_t> rho_order, double delta_i_max){
      //! delta_i is measured by computing the minimum distance between the point i
      /*! and any other point with higher OR EQUAL density (emphasis mine)
       */
      delta = scitbx::af::shared<double>(NN,delta_i_max);
      const double* Dij = distance_matrix.begin();
      for (int p=1; p < NN; ++p){ // first iteration through rho order
        int i = rho_order[p];
        for (int q=0; q < p; ++q){
          int j = rho_order[q];
          if (rho[j] >= rho[i] && Dij[i*NN+j] < delta[i]) {
            delta[i] = Dij[i*NN+j];
          }
        }
      }
      return delta;
    }


    public:
      scitbx::af::shared<std::size_t> rho;
      scitbx::af::shared<double> delta;

    private:
      scitbx::af::flex_double distance_matrix;
      double d_c;
      bool rho_set;
      std::size_t NN;

  };

}

using namespace boost::python;
namespace sx_clustering{
namespace boost_python { namespace {

  void
  sx_clustering_init_module() {
    using namespace boost::python;
    typedef return_value_policy<return_by_value> rbv;
    typedef default_call_policies dcp;

    class_<Rodriguez_Laio_clustering_2014>("Rodriguez_Laio_clustering_2014",no_init)
      .def(init<>())
      .def(init<double >((arg_("d_c"))))
      .def(init<scitbx::af::flex_double, double >(
          (arg_("distance_matrix"),arg_("d_c"))))
      .def("get_rho",&Rodriguez_Laio_clustering_2014::get_rho)
      .def("get_delta",&Rodriguez_Laio_clustering_2014::get_delta,(
          (arg_("rho_order"),arg_("delta_i_max"))))
    ;

    def("foo2",&sx_clustering::foo2);
  }

}
}} // namespace sx_clustering::boost_python::<anonymous>

BOOST_PYTHON_MODULE(sx_clustering_ext)
{
  sx_clustering::boost_python::sx_clustering_init_module();
}
