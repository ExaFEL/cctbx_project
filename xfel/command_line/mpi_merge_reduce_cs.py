# LIBTBX_SET_DISPATCHER_NAME mpi.cluster_merge_cs
from __future__ import division
import sys,time,os
from xfel.command_line import cxi_merge
from libtbx.utils import Usage, multi_out
from cctbx.array_family import flex
from libtbx import easy_pickle

from xfel.command_line.single_node_merge import get_observations

from xfel.command_line.single_node_merge import scaling_manager as scaling_manager_base

#MPI_Reduce Op of _add_all_frames_mpi
def mpi_merge_op(data0, data1, datatype):
  data0.n_accepted += data1.n_accepted
  data0.n_file_error += data1.n_file_error
  data0.n_low_corr += data1.n_low_corr
  data0.n_low_signal += data1.n_low_signal
  data0.n_processed += data1.n_processed
  data0.n_wrong_bravais += data1.n_wrong_bravais
  data0.n_wrong_cell += data1.n_wrong_cell

  for key in data0.failure_modes.keys():
    data0.failure_modes[key] = data0.failure_modes.get(key,0) + data1.failure_modes[key]

  for index, isigi in data1.ISIGI.iteritems() :
    if (index in data0.ISIGI):
      data0.ISIGI[index] += isigi
    else:
      data0.ISIGI[index] = isigi
  
  data0.completeness += data1.completeness
  #print "observations count increased by %d" % flex.sum(data.completeness)
  data0.completeness_predictions += data1.completeness_predictions
  #print "predictions count increased by %d" % flex.sum(data.completeness_predictions)
  data0.summed_N += data1.summed_N
  data0.summed_weight += data1.summed_weight
  data0.summed_wt_I += data1.summed_wt_I

  data0.corr_values.extend(data1.corr_values)
  data0.d_min_values.extend(data1.d_min_values)
  if not data0.params.short_circuit:
    data0.observations.extend(data1.observations)
  data0.rejected_fractions.extend(data1.rejected_fractions)
  data0.wavelength.extend(data1.wavelength)

  data0.uc_values.add_cells(data1.uc_values)
  return data0

class scaling_manager_mpi(scaling_manager_base):

  def mpi_initialize (self, file_names) :
    tar_list,self.integration_pickle_names = file_names
    self.t1 = time.time()

    assert self.params.backend == 'FS' # for prototyping rank=0 marshalling
    from xfel.cxi.merging_database_fs import manager
    db_mgr = manager(self.params)
    db_mgr.initialize_db(self.miller_set.indices())
    self.master_db_mgr = db_mgr

  def mpi_finalize (self) :
    t2 = time.time()
    print >> self.log, ""
    print >> self.log, "#" * 80
    print >> self.log, "FINISHED MERGING"
    print >> self.log, "  Elapsed time: %.1fs" % (t2 - self.t1)
    print >> self.log, "  %d of %d integration files were accepted" % (
      self.n_accepted, len(self.integration_pickle_names))
    print >> self.log, "  %d rejected due to wrong Bravais group" % \
      self.n_wrong_bravais
    print >> self.log, "  %d rejected for unit cell outliers" % \
      self.n_wrong_cell
    print >> self.log, "  %d rejected for low signal" % \
      self.n_low_signal
    print >> self.log, "  %d rejected due to up-front poor correlation under min_corr parameter" % \
      self.n_low_corr
    print >> self.log, "  %d rejected for file errors or no reindex matrix" % \
      self.n_file_error
    for key in self.failure_modes.keys():
      print >>self.log, "  %d rejected due to %s"%(self.failure_modes[key], key)

    checksum = self.n_accepted  + self.n_file_error \
               + self.n_low_corr + self.n_low_signal \
               + self.n_wrong_bravais + self.n_wrong_cell \
               + sum([val for val in self.failure_modes.itervalues()])
    assert checksum == len(self.integration_pickle_names)

    high_res_count = (self.d_min_values <= self.params.d_min).count(True)
    print >> self.log, "Of %d accepted images, %d accepted to %5.2f Angstrom resolution" % \
      (self.n_accepted, high_res_count, self.params.d_min)

    if self.params.raw_data.sdfac_refine:
      self.scale_errors()

    if self.params.raw_data.errors_from_sample_residuals:
      self.errors_from_residuals()

class run_manager(object):
 def initialize(self,args):
  import iotbx.phil
  from xfel.command_line.cxi_merge import master_phil
  if ("--help" in args) :
    iotbx.phil.parse(master_phil).show(attributes_level=2)
    return
  phil = iotbx.phil.process_command_line(args=args, master_string=master_phil).show()
  work_params = phil.work.extract()
  from xfel.merging.phil_validation import application
  application(work_params)

  if ((work_params.d_min is None) or
      (work_params.data is None) or
      ( (work_params.model is None) and work_params.scaling.algorithm != "mark1") ) :
    command_name = os.environ["LIBTBX_DISPATCHER_NAME"]
    raise Usage(command_name + " "
                "d_min=4.0 "
                "data=~/scratch/r0220/006/strong/ "
                "model=3bz1_3bz2_core.pdb")
  if ((work_params.rescale_with_average_cell) and
      (not work_params.set_average_unit_cell)) :
    raise Usage("If rescale_with_average_cell=True, you must also specify "+
      "set_average_unit_cell=True.")
  if [work_params.raw_data.sdfac_auto, work_params.raw_data.sdfac_refine, work_params.raw_data.errors_from_sample_residuals].count(True) > 1:
    raise Usage("Specify only one of sdfac_auto, sdfac_refine or errors_from_sample_residuals.")

  # Read Nat's reference model from an MTZ file.  XXX The observation
  # type is given as F, not I--should they be squared?  Check with Nat!
  log = open("%s.log" % work_params.output.prefix, "w")
  out = multi_out()
  out.register("log", log, atexit_send_to=None)
  out.register("stdout", sys.stdout)
  print >> out, "I model"
  if work_params.model is not None:
    from xfel.merging.general_fcalc import run as run_fmodel
    i_model = run_fmodel(work_params)
    work_params.target_unit_cell = i_model.unit_cell()
    work_params.target_space_group = i_model.space_group_info()
    i_model.show_summary()
  else:
    i_model = None

  print >> out, "Target unit cell and space group:"
  print >> out, "  ", work_params.target_unit_cell
  print >> out, "  ", work_params.target_space_group
  from xfel.command_line.cxi_merge import consistent_set_and_model
  self.miller_set, self.i_model = consistent_set_and_model(work_params,i_model)
  self.work_params = work_params
  self.frame_files = get_observations(work_params)
  self.out = out

 def other(self,scaler):
  out = self.out
  work_params = self.work_params
  miller_set = self.miller_set
  if scaler.n_accepted == 0:
    return None
  scaler.show_unit_cell_histograms()
  if (work_params.rescale_with_average_cell) :
    average_cell_abc = scaler.uc_values.get_average_cell_dimensions()
    average_cell = uctbx.unit_cell(list(average_cell_abc) +
      list(work_params.target_unit_cell.parameters()[3:]))
    work_params.target_unit_cell = average_cell
    print >> out, ""
    print >> out, "#" * 80
    print >> out, "RESCALING WITH NEW TARGET CELL"
    print >> out, "  average cell: %g %g %g %g %g %g" % \
      work_params.target_unit_cell.parameters()
    print >> out, ""
    assert False,"must do this step again with MPI"
    scaler.reset()
    scaler.scale_all(frame_files)
    scaler.show_unit_cell_histograms()
  print >> out, "\n"

  # Sum the observations of I and I/sig(I) for each reflection.
  sum_I = flex.double(miller_set.size(), 0.)
  sum_I_SIGI = flex.double(miller_set.size(), 0.)
  for i in xrange(miller_set.size()) :
    index = miller_set.indices()[i]
    if index in scaler.ISIGI :
      for t in scaler.ISIGI[index]:
        sum_I[i] += t[0]
        sum_I_SIGI[i] += t[1]

  miller_set_avg = miller_set.customized_copy(
    unit_cell=work_params.target_unit_cell)
  table1 = cxi_merge.show_overall_observations(
    obs=miller_set_avg,
    redundancy=scaler.completeness,
    redundancy_to_edge=scaler.completeness_predictions,
    summed_wt_I=scaler.summed_wt_I,
    summed_weight=scaler.summed_weight,
    ISIGI=scaler.ISIGI,
    n_bins=work_params.output.n_bins,
    title="Statistics for all reflections",
    out=out,
    work_params=work_params)
  print >> out, ""
  if work_params.model is not None:
    n_refl, corr = scaler.get_overall_correlation(sum_I)
  else:
    n_refl, corr = ((scaler.completeness > 0).count(True), 0)
  print >> out, "\n"
  table2 = cxi_merge.show_overall_observations(
    obs=miller_set_avg,
    redundancy=scaler.summed_N,
    redundancy_to_edge=scaler.completeness_predictions,
    summed_wt_I=scaler.summed_wt_I,
    summed_weight=scaler.summed_weight,
    ISIGI=scaler.ISIGI,
    n_bins=work_params.output.n_bins,
    title="Statistics for reflections where I > 0",
    out=out,
    work_params=work_params)
  #from libtbx import easy_pickle
  #easy_pickle.dump(file_name="stats.pickle", obj=stats)
  #stats.report(plot=work_params.plot)
  #miller_counts = miller_set_p1.array(data=stats.counts.as_double()).select(
  #  stats.counts != 0)
  #miller_counts.as_mtz_dataset(column_root_label="NOBS").mtz_object().write(
  #  file_name="nobs.mtz")
  if work_params.data_subsubsets.subsubset is not None and work_params.data_subsubsets.subsubset_total is not None:
    easy_pickle.dump("scaler_%d.pickle"%work_params.data_subsubsets.subsubset, scaler)
  explanation = """
Explanation:
Completeness       = # unique Miller indices present in data / # Miller indices theoretical in asymmetric unit
Asu. Multiplicity  = # measurements / # Miller indices theoretical in asymmetric unit
Obs. Multiplicity  = # measurements / # unique Miller indices present in data
Pred. Multiplicity = # predictions on all accepted images / # Miller indices theoretical in asymmetric unit"""
  print >> out, explanation
  mtz_file, miller_array = scaler.finalize_and_save_data()
  #table_pickle_file = "%s_graphs.pkl" % work_params.output.prefix
  #easy_pickle.dump(table_pickle_file, [table1, table2])
  loggraph_file = os.path.abspath("%s_graphs.log" % work_params.output.prefix)
  f = open(loggraph_file, "w")
  f.write(table1.format_loggraph())
  f.write("\n")
  f.write(table2.format_loggraph())
  f.close()
  result = cxi_merge.scaling_result(
    miller_array=miller_array,
    plots=scaler.get_plot_statistics(),
    mtz_file=mtz_file,
    loggraph_file=loggraph_file,
    obs_table=table1,
    all_obs_table=table2,
    n_reflections=n_refl,
    overall_correlation=corr)
  easy_pickle.dump("%s.pkl" % work_params.output.prefix, result)
  return result


if (__name__ == "__main__"):
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  size = comm.Get_size()

  merge_op = MPI.Op.Create(mpi_merge_op, commute=True)


  # set things up
  if rank == 0:
    run = run_manager()
    run.initialize(args=sys.argv[1:])
    scaler_master = scaling_manager_mpi(
      miller_set=run.miller_set,
      i_model=run.i_model,
      params=run.work_params,
      log=run.out)
    scaler_master.mpi_initialize(run.frame_files)

    transmitted_info = dict(file_names=run.frame_files,
                            miller_set=run.miller_set,
                            model = run.i_model,
                            params = run.work_params )
  else:
    transmitted_info = None

  transmitted_info = comm.bcast(transmitted_info, root = 0)
  scaler_worker = scaling_manager_mpi(transmitted_info["miller_set"],
                                        transmitted_info["model"],
                                        transmitted_info["params"],
                                        log = sys.stdout)
  # now actually do the work
  tar_file_names = transmitted_info["file_names"][0]
  assert scaler_worker.params.backend == 'FS' # only option that makes sense
  from xfel.cxi.merging_database_fs import manager2 as manager
  db_mgr = manager(scaler_worker.params)
  scaler_worker.finished_db_mgr = db_mgr

  if rank == 0:
    for ix in xrange(len(tar_file_names)):
      rankreq = comm.recv(source=MPI.ANY_SOURCE)
      comm.send(ix,dest=rankreq)
    for rankreq in range(size-1):
      rankreq = comm.recv(source=MPI.ANY_SOURCE)
      comm.send('endrun',dest=rankreq)
  else:
    while True:
      comm.send(rank,dest=0)
      idx = comm.recv(source=0)
      if idx == 'endrun': 
        scaler_worker.finished_db_mgr = db_mgr
        break
      scaler_worker.tar_to_scale_frame_adapter(tar_list=[tar_file_names[idx],], db_mgr=db_mgr)

      # might want to clean up a bit before returning
    del scaler_worker.log
#    del scaler_worker.params
    del scaler_worker.miller_set
    del scaler_worker.i_model
    del scaler_worker.reverse_lookup
  print "I AM FINISHED: RANK %d"%rank

  from time import time as tt

  print "REDUCE START RANK=%d TIME=%f"%(rank,tt())
  scaler_workers = comm.reduce(scaler_worker, op=merge_op, root=0)
  print "REDUCE STOP RANK=%d TIME=%f"%(rank,tt())
  comm.Barrier()
  MPI.Finalize()

  if rank == 0:
    print "ISIGI_LEN=%s START RANK=%d TIME=%f"%(len(scaler_master.ISIGI.values()),rank,tt())
    scaler_master._add_all_frames(scaler_workers)
    for call_instance in scaler_workers.finished_db_mgr.sequencer:
      if call_instance["call"] == "insert_frame":
        frame_id_zero_base = scaler_master.master_db_mgr.insert_frame(**call_instance["data"])
      elif call_instance["call"] == "insert_observation":
        call_instance["data"]['frame_id_0_base'] = [frame_id_zero_base] * len(call_instance["data"]['frame_id_0_base'])
        scaler_master.master_db_mgr.insert_observation(**call_instance["data"])
    print "ISIGI_LEN=%s STOP RANK=%d TIME=%f"%(len(scaler_master.ISIGI.values()),rank,tt())
    scaler_master.master_db_mgr.join() # database written, finalize the manager
    scaler_master.mpi_finalize()

    print "OK"
    run.other(scaler_master)
