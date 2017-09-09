#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <cfloat>
#include <cctbx/uctbx/determine_unit_cell/NCDist.h>
#include <scitbx/array_family/tiny.h>

namespace af = scitbx::af;

using namespace boost::python;

namespace ncdist2017 {
#include <NCDist.h>
}

namespace cctbx { namespace uctbx {

  double NCDist_wrapper(af::tiny<double,6> mm1,af::tiny<double,6> mm2){
    return NCDist(&mm1[0],&mm2[0]);
  }

  double NCDist2017_wrapper(af::tiny<double,6> mm1,af::tiny<double,6> mm2){
    return ncdist2017::NCDist(&mm1[0],&mm2[0]);
  }

}}


BOOST_PYTHON_MODULE(determine_unit_cell_ext)
{
  def ("NCDist",&cctbx::uctbx::NCDist_wrapper);
  def ("NCDist2017",&cctbx::uctbx::NCDist2017_wrapper);
}
