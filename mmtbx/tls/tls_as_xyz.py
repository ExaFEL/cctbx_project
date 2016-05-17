from __future__ import division
from scitbx import matrix
import math
import random
import mmtbx.utils
from scitbx.array_family import flex
import iotbx.pdb
from libtbx.utils import Sorry
from cctbx import adptbx
import sys
from mmtbx.tls import analysis
from mmtbx.tls import tools

random.seed(2679941)

def print_step(s, log):
  n = 50-len(s)
  print >> log, s, "*"*n

def run(pdb_file_name,
        n_models,
        log,
        output_file_name_prefix,
        eps=1.e-7):
  pdb_inp = iotbx.pdb.input(file_name = pdb_file_name)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  asc = pdb_hierarchy.atom_selection_cache()
  cs = pdb_inp.crystal_symmetry_from_cryst1()
  tls_extract = mmtbx.tls.tools.tls_from_pdb_inp(
    remark_3_records = pdb_inp.extract_remark_iii_records(3),
    pdb_hierarchy    = pdb_hierarchy)
  for i_group, tls_params_one_group in enumerate(tls_extract.tls_params):
    selection = asc.selection(tls_params_one_group.selection_string)
    pdb_hierarchy_sel = pdb_hierarchy.select(selection)
    xrs = pdb_hierarchy_sel.extract_xray_structure(crystal_symmetry=cs)
    deg_to_rad_scale = math.pi/180
    # Units: T[A], L[deg**2], S[A*deg]
    T = matrix.sym(sym_mat3=tls_params_one_group.t)
    L = matrix.sym(sym_mat3=tls_params_one_group.l)
    S = matrix.sqr(tls_params_one_group.s)
    tlso = tools.tlso(
      t      = T.as_sym_mat3(),
      l      = L.as_sym_mat3(),
      s      = S,
      origin = tls_params_one_group.origin)
    # sanity check
    if(not adptbx.is_positive_definite(tls_params_one_group.t, eps)):
      raise Sorry("T matrix is not positive definite.")
    if(not adptbx.is_positive_definite(tls_params_one_group.l, eps)):
      raise Sorry("L matrix is not positive definite.")
    r = analysis.run(T=T, L=L*(deg_to_rad_scale**2), S=S*deg_to_rad_scale,
      log=log).self_check()
    ensemble_generator_obj = ensemble_generator(
      tls_from_motions_object = r,
      pdb_hierarchy           = pdb_hierarchy_sel,
      xray_structure          = xrs,
      n_models                = n_models,
      origin                  = tls_params_one_group.origin,
      log                     = log)
    ensemble_generator_obj.write_pdb_file(
      file_name=output_file_name_prefix+"_ensemble_%s.pdb"%str(i_group))
    # get U from TLS
    u_from_tls = tools.uaniso_from_tls_one_group(
      tlso          = tlso,
      sites_cart    = xrs.sites_cart(),
      zeroize_trace = False)
    # get U from ensemble
    pdb_hierarchy_from_tls = pdb_hierarchy_sel.deep_copy()
    pdb_hierarchy_from_ens = pdb_hierarchy_sel.deep_copy()
    u_from_ens = tools.u_cart_from_ensemble(
      models = ensemble_generator_obj.states.root.models())
    for i in xrange(xrs.sites_cart().size()):
      print "atom %d:"%i
      print "  Ucart(from TLS):", ["%8.5f"%u for u in u_from_tls[i]]
      print "  Ucart(from ens):", ["%8.5f"%u for u in u_from_ens[i]]
    pdb_hierarchy_from_tls.atoms().set_uij(u_from_tls)
    pdb_hierarchy_from_ens.atoms().set_uij(u_from_ens)
    pdb_hierarchy_from_tls.write_pdb_file(
      file_name = output_file_name_prefix+"_u_from_tls_%s.pdb"%str(i_group),
      crystal_symmetry = cs)
    pdb_hierarchy_from_ens.write_pdb_file(
      file_name = output_file_name_prefix+"_u_from_ensemble_%s.pdb"%str(i_group),
      crystal_symmetry = cs)
  return ensemble_generator_obj

class ensemble_generator(object):
  def __init__(self,
               tls_from_motions_object,
               pdb_hierarchy,
               xray_structure,
               n_models,
               origin,
               log=None):
    if(log is None): log = sys.stdout
    xray_structure.convert_to_isotropic()
    xray_structure = xray_structure.set_b_iso(value=0)
    sites_cart = xray_structure.sites_cart()-origin
    self.states = mmtbx.utils.states(
      xray_structure = xray_structure,
      pdb_hierarchy  = pdb_hierarchy)
    r = tls_from_motions_object
    print >> log
    print >> log, "Generating ensemble of %d models:"%n_models
    for trial in xrange(n_models):
      print >> log, "model #%d"%trial
      dx0,dy0,dz0 = step_i__get_dxdydz(r=r, log = log)
      d_r_M_V  = formula_49(r=r, log = log)
      sites_cart_new = flex.vec3_double()
      for site_cart in sites_cart:
        r_L = r.R_ML.transpose() * site_cart
        d_r_M_L = step_i__compute_delta_L_r_dp(r=r,
          r_L=r_L, dx0=dx0,dy0=dy0,dz0=dz0)
        d_r_M = d_r_M_L + d_r_M_V # (42)
        sites_cart_new.append(matrix.col(site_cart) + d_r_M)
      self.states.add(sites_cart = sites_cart_new+origin)

  def write_pdb_file(self, file_name):
    self.states.write(file_name = file_name)

def step_i__get_dxdydz(r, log, eps=1.e-7):
  """
  Generation of shifts from screw rotations.
  """
  dx0, dy0, dz0 = 0, 0, 0
  if(abs(r.dx)>eps): dx0 = random.normalvariate(0,r.dx)
  if(abs(r.dy)>eps): dy0 = random.normalvariate(0,r.dy)
  if(abs(r.dz)>eps): dz0 = random.normalvariate(0,r.dz)
  return dx0, dy0, dz0

def step_i__compute_delta_L_r_dp(r, r_L, dx0, dy0, dz0):
  sxb,syb,szb = r.sx, r.sy, r.sz
  x,y,z = r_L
  cos, sin = math.cos, math.sin
  d_lx_r_L = matrix.col((
    sxb*dx0,
    (y-r.wy_lx)*(cos(dx0)-1) - (z-r.wz_lx)*sin(dx0),
    (y-r.wy_lx)*sin(dx0)     + (z-r.wz_lx)*(cos(dx0)-1)
    ))
  d_ly_r_L = matrix.col((
    (z-r.wz_ly)*sin(dy0)     + (x-r.wx_ly)*(cos(dy0)-1),
    syb*dy0,
    (z-r.wz_ly)*(cos(dy0)-1) - (x-r.wx_ly)*sin(dy0)
    ))
  d_lz_r_L = matrix.col((
    (x-r.wx_lz)*(cos(dz0)-1) - (y-r.wy_lz)*sin(dz0),
    (x-r.wx_lz)*sin(dz0)     + (y-r.wy_lz)*(cos(dz0)-1),
    szb*dz0
    ))

  d_r_L = d_lx_r_L + d_ly_r_L + d_lz_r_L
  d_r_M = r.R_ML * d_r_L
  return d_r_M

def formula_49(r, log, eps=1.e-7):
  """
  Generate shifts from group translation.
  """
  tx0, ty0, tz0 = 0, 0, 0
  if(r.tx>eps): tx0 = random.normalvariate(0,r.tx)
  if(r.ty>eps): ty0 = random.normalvariate(0,r.ty)
  if(r.tz>eps): tz0 = random.normalvariate(0,r.tz)
  d_r_V = matrix.col((tx0,ty0,tz0))
  return r.R_MV * d_r_V
