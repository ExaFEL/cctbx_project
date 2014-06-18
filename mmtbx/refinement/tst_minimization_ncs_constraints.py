from __future__ import division
from  iotbx.pdb.multimer_reconstruction import ncs_group_object
from iotbx.pdb.multimer_reconstruction import multimer
import mmtbx.refinement.minimization_ncs_constraints
from libtbx.test_utils import approx_equal
import mmtbx.refinement.adp_refinement
from scitbx.array_family import flex
from libtbx import adopt_init_args
import mmtbx.utils.ncs_utils as nu
import mmtbx.f_model
import mmtbx.utils
import iotbx.pdb
import os
import sys


ncs_1_copy="""\
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.496590 -0.643597  0.582393        0.00000
MTRIX2   2  0.867925  0.376088 -0.324443        0.00000
MTRIX3   2 -0.010221  0.666588  0.745356        0.00000
MTRIX1   3 -0.317946 -0.173437  0.932111        0.00000
MTRIX2   3  0.760735 -0.633422  0.141629        0.00000
MTRIX3   3  0.565855  0.754120  0.333333        0.00000
ATOM      1  N   THR A   1       8.111  11.079  10.645  1.00 20.00           N
ATOM      2  CA  THR A   1       8.000   9.721  10.125  1.00 20.00           C
ATOM      3  C   THR A   1       8.075   8.693  11.249  1.00 20.00           C
ATOM      4  O   THR A   1       8.890   8.817  12.163  1.00 20.00           O
ATOM      5  CB  THR A   1       9.101   9.420   9.092  1.00 20.00           C
ATOM      6  OG1 THR A   1       9.001  10.342   8.000  1.00 20.00           O
ATOM      7  CG2 THR A   1       8.964   8.000   8.565  1.00 20.00           C
TER
"""


class ncs_minimization_test(object):

  def __init__(self,
               n_macro_cycle,
               sites,
               u_iso,
               transformations,
               finite_grad_differences_test,
               use_geometry_restraints,
               shake_site_mean_distance = 1.5,
               d_min = 2,
               shake_angles_sigma = 0.035,
               shake_translation_sigma = 0.5):
    """ create temp test files and data for tests """
    adopt_init_args(self, locals())
    self.test_files_names = [] # collect names of files for cleanup
    # 1 NCS copy: starting template to generate whole asu; place into P1 box
    pdb_inp = iotbx.pdb.input(source_info=None, lines=ncs_1_copy)
    pdb_obj = iotbx.pdb.hierarchy.input(pdb_string=ncs_1_copy)
    mtrix_object = pdb_inp.process_mtrix_records()
    ph = pdb_inp.construct_hierarchy()
    xrs = pdb_inp.xray_structure_simple()
    xrs_one_ncs = xrs.orthorhombic_unit_cell_around_centered_scatterers(
      buffer_size=8)
    ph.adopt_xray_structure(xrs_one_ncs)
    of = open("one_ncs_in_asu.pdb", "w")
    print >> of, mtrix_object.as_pdb_string()
    print >> of, ph.as_pdb_string(crystal_symmetry=xrs_one_ncs.crystal_symmetry())
    of.close()
    # 1 NCS copy -> full asu (expand NCS). This is the answer-structure
    m = multimer(file_name="one_ncs_in_asu.pdb",
                 round_coordinates=False,
                 reconstruction_type='cau',error_handle=True,eps=1e-2)
    assert m.number_of_transforms == 2, m.number_of_transforms
    xrs_asu = m.assembled_multimer.extract_xray_structure(
      crystal_symmetry = xrs_one_ncs.crystal_symmetry())
    m.write("full_asu.pdb")
    # force ASU none-rounded coordinates into xray structure
    xrs_asu.set_sites_cart(m.sites_cart())
    assert xrs_asu.crystal_symmetry().is_similar_symmetry(
      xrs_one_ncs.crystal_symmetry())
    # Generate Fobs from answer structure
    f_obs = abs(xrs_asu.structure_factors(d_min=d_min, algorithm="direct").f_calc())
    r_free_flags = f_obs.generate_r_free_flags()
    mtz_dataset = f_obs.as_mtz_dataset(column_root_label="F-obs")
    mtz_dataset.add_miller_array(
      miller_array=r_free_flags,
      column_root_label="R-free-flags")
    mtz_object = mtz_dataset.mtz_object()
    mtz_object.write(file_name = "data.mtz")
    # Shake structure - subject to refinement input
    xrs_shaken = xrs_one_ncs.deep_copy_scatterers()
    if sites: xrs_shaken.shake_sites_in_place(
      mean_distance=shake_site_mean_distance)
    if self.u_iso:
      u_random = flex.random_double(xrs_shaken.scatterers().size())
      xrs_shaken = xrs_shaken.set_u_iso(values=u_random)
    if self.transformations:
      transforms_obj = ncs_group_object()
      transforms_obj.build_ncs_obj_from_pdb_ncs(
      transform_info = mtrix_object,
      pdb_hierarchy_inp = pdb_obj)
      x = nu.concatenate_rot_tran(transforms_obj)
      x = nu.shake_transformations(
        x = x,
        shake_angles_sigma=self.shake_angles_sigma,
        shake_translation_sigma=self.shake_translation_sigma)
      transforms_obj = nu.separate_rot_tran(x,transforms_obj)
      mtrix_object = transforms_obj.build_MTRIX_object()
    ph.adopt_xray_structure(xrs_shaken)
    of = open("one_ncs_in_asu_shaken.pdb", "w")
    print >> of, mtrix_object.as_pdb_string()
    print >> of, ph.as_pdb_string(crystal_symmetry=xrs.crystal_symmetry())
    of.close()
    self.f_obs = f_obs
    self.r_free_flags = r_free_flags
    self.xrs_one_ncs = xrs_one_ncs
    # Get restraints manager
    self.grm = None
    pdb_str = m.assembled_multimer.as_pdb_string(
      crystal_symmetry=xrs_one_ncs.crystal_symmetry())
    self.iso_restraints = None
    if(self.use_geometry_restraints):
      self.grm = nu.get_restraints_manager(pdb_string=pdb_str)
      if(self.u_iso):
        temp = mmtbx.refinement.adp_refinement.adp_restraints_master_params
        self.iso_restraints = temp.extract().iso

  def run_test(self):
    # Refinement
    params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
    params.algorithm = "direct"
    # Get the xray_structure of the shaken ASU
    m_shaken = multimer(
      file_name="one_ncs_in_asu_shaken.pdb",
      round_coordinates=False,
      reconstruction_type='cau',error_handle=True,eps=1e-2)
    xrs_shaken_asu = m_shaken.assembled_multimer.extract_xray_structure(
      crystal_symmetry=self.xrs_one_ncs.crystal_symmetry())
    # force non-rounded coordinates into xray structure
    xrs_shaken_asu.set_sites_cart(m_shaken.sites_cart())
    # Save the shaken ASU for inspection
    m_shaken.write(
      pdb_output_file_name='asu_shaken.pdb',
      crystal_symmetry=self.xrs_one_ncs.crystal_symmetry())
    self.transforms_obj = m_shaken.transforms_obj
    # Create a boolean selection string for selecting chains in NCS
    ncs_selection = self.transforms_obj.ncs_atom_selection
    self.ncs_atom_selection = self.transforms_obj.ncs_atom_selection
    assert ncs_selection.count(True) > 0
    self.fmodel = mmtbx.f_model.manager(
      f_obs                        = self.f_obs,
      r_free_flags                 = self.r_free_flags,
      xray_structure               = xrs_shaken_asu,
      sf_and_grads_accuracy_params = params,
      target_name                  = "ls_wunit_k1")
    r_start = self.fmodel.r_work()
    assert r_start > 0.1, r_start
    print "start r_factor: %6.4f" % r_start
    for macro_cycle in xrange(self.n_macro_cycle):
      if self.transformations and \
              not self.transforms_obj.transform_chain_assignment: continue
      data_weight = None
      if(self.use_geometry_restraints):
        data_weight = nu.get_weight(self)
      target_and_grads_object = mmtbx.refinement.minimization_ncs_constraints.\
        target_function_and_grads_reciprocal_space(
          fmodel                 = self.fmodel,
          transforms_obj         = self.transforms_obj,
          restraints_manager     = self.grm,
          data_weight            = data_weight,
          refine_sites           = self.sites,
          refine_u_iso           = self.u_iso,
          refine_transformations = self.transformations,
          iso_restraints         = self.iso_restraints)
      minimized = mmtbx.refinement.minimization_ncs_constraints.lbfgs(
        target_and_grads_object      = target_and_grads_object,
        xray_structure               = self.fmodel.xray_structure,
        transforms_obj               = self.transforms_obj,
        finite_grad_differences_test = self.finite_grad_differences_test,
        max_iterations               = 60,
        refine_sites                 = self.sites,
        refine_u_iso                 = self.u_iso,
        refine_transformations       = self.transformations)
      refine_type = 'adp'*self.u_iso + 'sites'*self.sites \
                    + 'transformation'*self.transformations
      outstr = "  macro_cycle {0:3} ({1})   r_factor: {2:6.4f}   " + \
            self.finite_grad_differences_test * \
            "finite_grad_difference_val: {3:.4f}"
      print outstr.format(
        macro_cycle, refine_type,self.fmodel.r_work(),
        minimized.finite_grad_difference_val)
      assert (minimized.finite_grad_difference_val < 1.0e-3)
      assert approx_equal(self.fmodel.r_work(), target_and_grads_object.fmodel.r_work())
      # break test if r_work is very small
      if target_and_grads_object.fmodel.r_work() < 1.0e-6: break
    # check results
    if(self.u_iso):
      assert approx_equal(self.fmodel.r_work(), 0, 1.e-5)
    elif(self.sites):
      if(self.use_geometry_restraints):
        assert approx_equal(self.fmodel.r_work(), 0, 0.0001)
      else:
        assert approx_equal(self.fmodel.r_work(), 0, 1.e-5)
    elif self.transformations:
        assert approx_equal(self.fmodel.r_work(), 0, 0.0001)
    else: assert 0
    # output refined model
    xrs_refined = self.fmodel.xray_structure
    m_shaken.assembled_multimer.adopt_xray_structure(self.fmodel.xray_structure)
    output_file_name = "refined_u_iso%s_sites%s.pdb"%(str(self.u_iso),
      str(self.sites))
    m_shaken.write(output_file_name)
    self.test_files_names.append(output_file_name)
    # check final model
    if(not self.use_geometry_restraints):
      # XXX fix later for case self.use_geometry_restraints=True
      pdb_inp_answer = iotbx.pdb.input(source_info=None, lines=ncs_1_copy)
      pdb_inp_refined = iotbx.pdb.input(file_name=output_file_name)
      xrs1 = pdb_inp_answer.xray_structure_simple()
      xrs2 = pdb_inp_refined.xray_structure_simple().select(ncs_selection)
      mmtbx.utils.assert_xray_structures_equal(
        x1 = xrs1,
        x2 = xrs2,
        sites = False)
      delta = flex.vec3_double([xrs1.center_of_mass()]*xrs2.scatterers().size())-\
              flex.vec3_double([xrs2.center_of_mass()]*xrs2.scatterers().size())
      xrs2.set_sites_cart(sites_cart = xrs2.sites_cart()+delta)
      mmtbx.utils.assert_xray_structures_equal(
        x1 = xrs1,
        x2 = xrs2)

  def clean_up_temp_test_files(self):
    """delete temporary test files """
    files_to_delete = ['one_ncs_in_asu.pdb','full_asu.pdb',
                       'one_ncs_in_asu_shaken.pdb',
                       'asu_shaken.pdb','data.mtz']
    files_to_delete.extend(self.test_files_names)
    for fn in files_to_delete:
      if os.path.isfile(fn): os.remove(fn)

def exercise_without_geometry_restaints():
  print 'Running ',sys._getframe().f_code.co_name
  for sites, u_iso, n_macro_cycle in [(True, False, 100), (False, True, 50)]:
    t = ncs_minimization_test(
      n_macro_cycle   = n_macro_cycle,
      sites           = sites,
      u_iso           = u_iso,
      transformations = False,
      finite_grad_differences_test = True,
      use_geometry_restraints = False,
      shake_site_mean_distance = 0.5,
      d_min = 2.0)
    t.run_test()
    t.clean_up_temp_test_files()

def exercise_site_refinement():
  print 'Running ',sys._getframe().f_code.co_name
  t = ncs_minimization_test(
    n_macro_cycle   = 60,
    sites           = True,
    u_iso           = False,
    transformations = False,
    finite_grad_differences_test = True,
    use_geometry_restraints = True,
    shake_site_mean_distance = 1.5,
    d_min = 3)
  t.run_test()
  t.clean_up_temp_test_files()

def exercise_u_iso_refinement():
  print 'Running ',sys._getframe().f_code.co_name
  t = ncs_minimization_test(
    n_macro_cycle   = 50,
    sites           = False,
    u_iso           = True,
    transformations = False,
    finite_grad_differences_test = True,
    use_geometry_restraints = True,
    shake_site_mean_distance = 1.5,
    d_min = 2)
  t.run_test()
  t.clean_up_temp_test_files()

def exercise_transformation_refinement():
  """  Test transformation refinement  """
  print 'Running ',sys._getframe().f_code.co_name
  t = ncs_minimization_test(
    n_macro_cycle   = 40,
    sites           = False,
    u_iso           = False,
    transformations = True,
    finite_grad_differences_test = True,
    use_geometry_restraints = True,
    shake_angles_sigma = 0.032,
    shake_translation_sigma = 0.5,
    d_min = 2)
  t.run_test()
  t.clean_up_temp_test_files()

if __name__ == "__main__":
  exercise_without_geometry_restaints()
  exercise_site_refinement()
  exercise_u_iso_refinement()
  exercise_transformation_refinement()
