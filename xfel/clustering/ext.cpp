// Author:Nick Sauter
#include <cctbx/boost_python/flex_fwd.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>

namespace sx_clustering {
  void foo2(){
    std::cout<<"HELLO foo2"<<std::endl;
  }

  class Rodriguez_Laio_Clustering_2014 {
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

    class_<Rodriguez_Laio_Clustering_2014>("Rodriguez_Laio_Clustering_2014",init<>())
    ;

    def("foo2",&sx_clustering::foo2);
  }

}
}} // namespace lunus::boost_python::<anonymous>

BOOST_PYTHON_MODULE(sx_clustering_ext)
{
  sx_clustering::boost_python::sx_clustering_init_module();
}
