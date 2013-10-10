
from __future__ import division
from mmtbx.validation import validation, atoms, atom_info
from libtbx import easy_run
import libtbx.load_env
import string
import sys

class clash (atoms) :
  __clash_attr__ = [
    "overlap",
    "probe_type",
    "max_b_factor",
  ]
  __slots__ = atoms.__slots__ + __clash_attr__

  @staticmethod
  def header () :
    return "%-20s %-20s  %7s" % ("Atom 1", "Atom 2", "Overlap")

  def id_str (self, spacer=" ") :
    return "%s%s%s" % (self.atoms_info[0].id_str(), spacer,
      self.atoms_info[1].id_str())

  def id_str_no_atom_name (self) :
    return "%s %s" % (self.atoms_info[0].id_str()[0:11],
      self.atoms_info[1].id_str()[0:11])

  def format_old (self) :
    return "%s :%.3f" % (self.id_str(), self.overlap)

  def as_string (self) :
    return "%-20s %-20s  %7.3f" % (self.atoms_info[0].id_str(),
      self.atoms_info[1].id_str(), self.overlap)

  def __cmp__ (self, other) : # sort in descending order
    return cmp(self.overlap, other.overlap)

class clashscore(validation):
  __slots__ = validation.__slots__ + [
    "clashscore",
    "clashscore_b_cutoff",
    "clash_dict",
    "clash_dict_b_cutoff",
    "list_dict",
    "b_factor_cutoff",
    "probe_file",
  ]
  program_description = "Analyze clashscore for protein model"

  def get_result_class (self) : return clash

  def __init__ (self,
      pdb_hierarchy,
      keep_hydrogens=True,
      nuclear=False,
      force_unique_chain_ids=False,
      time_limit=120,
      b_factor_cutoff=None,
      save_probe_unformatted_file=None,
      save_modified_hierarchy=False,
      verbose=False,
      out=sys.stdout) :
    validation.__init__(self)
    self.b_factor_cutoff = b_factor_cutoff
    self.clashscore = None
    self.clashscore_b_cutoff = None
    self.clash_dict = {}
    self.clash_dict_b_cutoff = {}
    self.list_dict = {}
    self.probe_file = None
    if (not libtbx.env.has_module(name="probe")):
      raise RuntimeError(
        "Probe could not be detected on your system.  Please make sure "+
        "Probe is in your path.\nProbe is available at "+
        "http://kinemage.biochem.duke.edu/")
    if verbose:
      if not nuclear:
        print "\nUsing electron cloud x-H distances and vdW radii"
      else:
        print "\nUsing nuclear cloud x-H distances and vdW radii"
    import iotbx.pdb.hierarchy
    from scitbx.array_family import flex
    n_models = len(pdb_hierarchy.models())
    for i_mod, model in enumerate(pdb_hierarchy.models()):
      r = iotbx.pdb.hierarchy.root()
      mdc = model.detached_copy()
      r.append_model(mdc)
      tmp_r = r
      # removed old style SEGID handling for compatibility with Probe
      # 130622 - JJH
      #from mmtbx import utils
      #bare_chains = \
      #  utils.find_bare_chains_with_segids(pdb_hierarchy=r)
      #if bare_chains:
      #  tmp_r = r.deep_copy()
      #  tmp_r.atoms().reset_i_seq()
      #  seg_dict = utils.seg_id_to_chain_id(pdb_hierarchy=tmp_r)
      #  rename_txt = utils.assign_chain_ids(pdb_hierarchy=tmp_r,
      #                                      seg_dict=seg_dict)
      #else:
      #  tmp_r = r
      #duplicate_chain_ids = \
      #  utils.check_for_duplicate_chain_ids(pdb_hierarchy=tmp_r)
      #if duplicate_chain_ids:
      #  utils.force_unique_chain_ids(pdb_hierarchy=tmp_r)
      if keep_hydrogens:
        elements = tmp_r.atoms().extract_element()
        h_count = elements.count(' H') + elements.count(' D')
        if h_count > 0:
          has_hd = True
        else:
          has_hd = False
        # if no hydrogens present, force addition for clashscore
        # calculation
        if not has_hd:
          if verbose:
            print "\nNo H/D atoms detected - forcing hydrogen addition!\n"
          keep_hydrogens = False
      input_str = tmp_r.as_pdb_string()
      occ_max = flex.max(tmp_r.atoms().extract_occ())
      pcm = probe_clashscore_manager(
        pdb_string=input_str,
        keep_hydrogens=keep_hydrogens,
        nuclear=nuclear,
        time_limit=time_limit,
        largest_occupancy=occ_max,
        b_factor_cutoff=b_factor_cutoff,
        verbose=verbose)
      if (save_modified_hierarchy) :
        self.pdb_hierarchy = pdb.hierarchy.\
          input(pdb_string=pcm.h_pdb_string).hierarchy
      self.clash_dict[model.id] = pcm.clashscore
      self.clash_dict_b_cutoff[model.id] = pcm.clashscore_b_cutoff
      self.list_dict[model.id] = pcm.bad_clashes
      if (n_models == 1) or (self.clashscore is None) :
        self.results = pcm.bad_clashes
        self.clashscore = pcm.clashscore
        self.clashscore_b_cutoff = pcm.clashscore_b_cutoff
      if (save_probe_unformatted_file is not None) and (n_models == 1) :
        open(save_probe_unformatted_file, "w").write(pcm.probe_unformatted)
        self.probe_file = save_probe_unformatted_file

  def get_clashscore(self):
    return self.clashscore

  def get_clashscore_b_cutoff(self):
    return self.clashscore_b_cutoff

  def show_old_output (self, out=sys.stdout, verbose=False) :
    self.print_clashlist_old(out=out)
    self.show_summary(out=out)

  def show_summary (self, out=sys.stdout, prefix="") :
    if (len(self.clash_dict) == 1) :
      print >> out, prefix + "clashscore = %.2f" % self.clash_dict['']
      if self.clash_dict_b_cutoff[''] is not None:
        print >> out, "clashscore (B factor cutoff = %d) = %f" % \
          (self.b_factor_cutoff,
           self.clash_dict_b_cutoff[''])
    else:
      for k in sorted(self.clash_dict.keys()) :
        print >> out, prefix + "MODEL %s clashscore = %.2f" % (k,
          self.clash_dict[k])
        if self.clash_dict_b_cutoff[k] is not None:
          print >> out, "MODEL%s clashscore (B factor cutoff = %d) = %f" % \
            (k, self.b_factor_cutoff, self.clash_dict_b_cutoff[k])

  def print_clashlist_old (self, out=sys.stdout):
    for k in self.list_dict.keys():
      if k == '':
        print >> out, "Bad Clashes >= 0.4 Angstrom:"
        for result in self.list_dict[k] :
          print >> out, result.format_old()
      else:
        print >> out, "Bad Clashes >= 0.4 Angstrom MODEL%s" % k
        for result in self.list_dict[k] :
          print >> out, result.format_old()

  def show (self, out=sys.stdout, prefix="", outliers_only=None, verbose=None) :
    if (len(self.clash_dict) == 1) :
      for result in self.list_dict[''] :
        print >> out, prefix + str(result)
    else :
      for k in self.list_dict.keys():
        for result in self.list_dict[k] :
          print >> out, prefix + str(result)
    self.show_summary(out=out, prefix=prefix)

  def as_coot_data (self) :
    data = []
    for result in self.results :
      if result.is_outlier() :
        data.append((result.atoms_info[0].id_str(),
          result.atoms_info[1].id_str(), result.overlap, result.xyz))
    return data

class probe_clashscore_manager(object):
  def __init__(self,
               pdb_string,
               keep_hydrogens=True,
               nuclear=False,
               time_limit=120,
               largest_occupancy=10,
               b_factor_cutoff=None,
               verbose=False):
    assert (libtbx.env.has_module(name="reduce") and
            libtbx.env.has_module(name="probe"))

    self.b_factor_cutoff = b_factor_cutoff
    ogt = 10
    blt = self.b_factor_cutoff
    if largest_occupancy < ogt:
      ogt = largest_occupancy

    self.trim = "phenix.reduce -quiet -trim -"
    self.probe_atom_b_factor = None
    if not nuclear:
      self.build = "phenix.reduce -oh -his -flip -pen9999" +\
                   " -keep -allalt -limit%d -" % time_limit
      self.probe_txt = \
        'phenix.probe -u -q -mc -het -once "ogt%d not water" "ogt%d" -' % \
          (ogt, ogt)
      self.probe_atom_txt = \
        'phenix.probe -q -mc -het -dumpatominfo "ogt%d not water" -' % ogt
      if blt is not None:
        self.probe_atom_b_factor = \
          'phenix.probe -q -mc -het -dumpatominfo "blt%d ogt%d not water" -' % \
            (blt, ogt)
    else: #use nuclear distances
      self.build = "phenix.reduce -oh -his -flip -pen9999" +\
                   " -keep -allalt -limit%d -nuc -" % time_limit
      self.probe_txt = \
        'phenix.probe -u -q -mc -het -once -nuclear' +\
          ' "ogt%d not water" "ogt%d" -' % (ogt, ogt)
      self.probe_atom_txt = \
        'phenix.probe -q -mc -het -dumpatominfo -nuclear' +\
          ' "ogt%d not water" -' % ogt
      if blt is not None:
        self.probe_atom_b_factor = \
          'phenix.probe -q -mc -het -dumpatominfo -nuclear' +\
            ' "blt%d ogt%d not water" -' % (blt, ogt)

    if not keep_hydrogens:
      h_pdb_string = self.run_reduce(pdb_string)
    else:
      if verbose:
        print "\nUsing input model H/D atoms...\n"
      h_pdb_string = pdb_string
    self.h_pdb_string = h_pdb_string
    self.run_probe_clashscore(self.h_pdb_string)

  def run_reduce(self, pdb_string):
    clean_out = easy_run.fully_buffered(self.trim,
                  stdin_lines=pdb_string)
    if (clean_out.return_code != 0) :
      raise RuntimeError("Reduce crashed with command '%s' - dumping stderr:\n%s"
        % (self.trim, "\n".join(clean_out.stderr_lines)))
    build_out = easy_run.fully_buffered(self.build,
                  stdin_lines=clean_out.stdout_lines)
    if (build_out.return_code != 0) :
      raise RuntimeError("Reduce crashed with command '%s' - dumping stderr:\n%s"
        % (self.build, "\n".join(build_out.stderr_lines)))
    reduce_str = string.join(build_out.stdout_lines, '\n')
    return reduce_str

  #def update_clashscore(self, pdb_string):
  #  self.run_probe_clashscore(pdb_string)

  def run_probe_clashscore(self, pdb_string):
    clash_hash = {}
    hbond_hash = {}
    clashscore = None
    probe_out = easy_run.fully_buffered(self.probe_txt,
      stdin_lines=pdb_string)
    if (probe_out.return_code != 0) :
      raise RuntimeError("Probe crashed - dumping stderr:\n%s" %
        "\n".join(probe_out.stderr_lines))
    probe_unformatted = probe_out.stdout_lines
    self.probe_unformatted = "\n".join(probe_unformatted)
    for line in probe_unformatted:
      name, pat, type, srcAtom, targAtom, min_gap, gap, \
      kissEdge2BullsEye, dot2BE, dot2SC, spike, score, stype, \
      ttype, x, y, z, sBval, tBval = line.split(":")
      atom1 = decode_atom_string(srcAtom)
      atom2 = decode_atom_string(targAtom)
      if (cmp(srcAtom,targAtom) < 0):
        atoms = [ atom1, atom2 ]
      else:
        atoms = [ atom2, atom1 ]
      gap = float(gap)
      x, y, z = float(x), float(y), float(z)
      clash_obj = clash(
        atoms_info=atoms,
        overlap=gap,
        probe_type=type,
        outlier=abs(gap) > 0.4,
        max_b_factor=max(float(sBval), float(tBval)),
        xyz=(x,y,z))
      key = clash_obj
      if (type == "so" or type == "bo"):
        if (gap <= -0.4):
          if (key in clash_hash) :
            if (gap < clash_hash[key].overlap):
              clash_hash[key] = clash_obj
          else :
            clash_hash[key] = clash_obj
      elif (type == "hb"):
        if (key in hbond_hash) :
          if (gap < hbond_hash[key].overlap):
            hbond_hash[key] = clash_obj
        else :
          hbond_hash[key] = clash_obj

    #sort the output
    temp = []
    for k in clash_hash.keys():
      if not k in hbond_hash:
        temp.append(clash_hash[k])
    self.n_clashes = len(temp)
    if self.b_factor_cutoff is not None:
      clashes_b_cutoff = 0
      for clash_obj in temp:
        if clash_obj.max_b_factor < self.b_factor_cutoff:
          clashes_b_cutoff += 1
      self.n_clashes_b_cutoff = clashes_b_cutoff
    used = []
    self.bad_clashes = []
    for clash_obj in sorted(temp) :
      test_key = clash_obj.id_str_no_atom_name()
      if test_key not in used:
        used.append(test_key)
        self.bad_clashes.append(clash_obj)
    probe_info = easy_run.fully_buffered(self.probe_atom_txt,
      stdin_lines=pdb_string).raise_if_errors().stdout_lines
    if (len(probe_info) == 0) :
      raise RuntimeError("Empty PROBE output.")
    self.n_atoms = 0
    for line in probe_info :
      dump, n_atoms = line.split(":")
    self.n_atoms = int(n_atoms)
    self.natoms_b_cutoff = None
    if self.probe_atom_b_factor is not None:
      probe_info_b_factor = easy_run.fully_buffered(self.probe_atom_b_factor,
        stdin_lines=pdb_string).raise_if_errors().stdout_lines
      for line in probe_info_b_factor :
        dump_b, natoms_b_cutoff = line.split(":")
      self.natoms_b_cutoff = int(natoms_b_cutoff)
    if self.n_atoms == 0:
      clashscore = 0.0
    else:
      clashscore = (self.n_clashes * 1000) / self.n_atoms
    self.clashscore = clashscore
    clashscore_b_cutoff = None
    if self.natoms_b_cutoff is not None and self.natoms_b_cutoff == 0:
      clashscore_b_cutoff = 0.0
    elif self.natoms_b_cutoff is not None:
      clashscore_b_cutoff = \
        (self.n_clashes_b_cutoff*1000) / self.natoms_b_cutoff
    self.clashscore_b_cutoff = clashscore_b_cutoff

def decode_atom_string (atom_str) :
  # Example:
  # ' A  49 LEU HD11B'
  return atom_info(
    chain_id=atom_str[0:2],
    resseq=atom_str[2:6],
    icode=atom_str[6],
    resname=atom_str[7:10],
    altloc=atom_str[15],
    name=atom_str[11:15])
