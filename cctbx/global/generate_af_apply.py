from generate_af_algebras import *

def write_copyright():
  print \
"""/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (Ralf W. Grosse-Kunstleve)

   *****************************************************
   THIS IS AN AUTOMATICALLY GENERATED FILE. DO NOT EDIT.
   *****************************************************

   Generated by:
     %s
 */""" % (sys.argv[0],)

def one_type(array_type_name):
  f = open("%s_apply.h" % (array_type_name,), "w")
  sys.stdout = f
  write_copyright()
  include_array_type_name = array_type_name
  if (array_type_name == "ref"):
    include_array_type_name = "versa"
  generic_include = "functors"
  if (base_array_type_name(array_type_name) == "tiny"):
    generic_include = "operators"
  print """
#ifndef CCTBX_ARRAY_FAMILY_%s_APPLY_H
#define CCTBX_ARRAY_FAMILY_%s_APPLY_H

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <cctbx/array_family/type_holder.h>
#include <cctbx/array_family/%s.h>
#include <cctbx/array_family/generic_array_%s.h>

namespace cctbx { namespace af {
""" % ((array_type_name.upper(),) * 2 + (
    include_array_type_name, generic_include))

  generate_unary_apply(array_type_name)

  print """}} // namespace cctbx::af

#endif // DOXYGEN_SHOULD_SKIP_THIS

#endif // CCTBX_ARRAY_FAMILY_%s_APPLY_H""" % (array_type_name.upper(),)
  sys.stdout = sys.__stdout__
  f.close()

def run():
  for array_type_name in (
    "tiny_plain", "tiny",
    "small_plain", "small",
    "shared_plain", "shared",
    "versa_plain", "versa",
    "ref"):
    one_type(array_type_name)

if (__name__ == "__main__"):
  run()
