# LIBTBX_SET_DISPATCHER_NAME mpi.merge
from __future__ import division
from xfel.command_line.cxi_merge import run
import sys,time
from xfel.command_line import cxi_merge
from libtbx import Auto
from scitbx import matrix

from xfel.command_line.single_node_merge import get_observations
cxi_merge.get_observations = get_observations

from xfel.command_line.single_node_merge import scaling_manager as scaling_manager_base
class scaling_manager(scaling_manager_base):

  def scale_all (self, file_names) :
    tar_list,integration_pickle_names = file_names
    t1 = time.time()
    if self.params.backend == 'MySQL':
      from xfel.cxi.merging_database import manager
    elif self.params.backend == 'SQLite':
      from xfel.cxi.merging_database_sqlite3 import manager
    else:
      from xfel.cxi.merging_database_fs import manager

    db_mgr = manager(self.params)
    db_mgr.initialize_db(self.miller_set.indices())

    # Unless the number of requested processes is greater than one,
    # try parallel multiprocessing on a parallel host.  Block until
    # all database commands have been processed.
    nproc = self.params.nproc
    if (nproc is None) or (nproc is Auto):
      import libtbx.introspection
      nproc = libtbx.introspection.number_of_processors()
    if nproc > 1:
      try :
        import multiprocessing
        self._scale_all_parallel(tar_list, db_mgr)
      except ImportError, e :
        print >> self.log, \
          "multiprocessing module not available (requires Python >= 2.6)\n" \
          "will scale frames serially"
        self._scale_all_serial(tar_list, db_mgr)
    else:
      self._scale_all_serial(tar_list, db_mgr)
    db_mgr.join()

    t2 = time.time()
    print >> self.log, ""
    print >> self.log, "#" * 80
    print >> self.log, "FINISHED MERGING"
    print >> self.log, "  Elapsed time: %.1fs" % (t2 - t1)
    print >> self.log, "  %d of %d integration files were accepted" % (
      self.n_accepted, len(integration_pickle_names))
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
    assert checksum == len(integration_pickle_names)

    high_res_count = (self.d_min_values <= self.params.d_min).count(True)
    print >> self.log, "Of %d accepted images, %d accepted to %5.2f Angstrom resolution" % \
      (self.n_accepted, high_res_count, self.params.d_min)

    if self.params.raw_data.sdfac_refine:
      self.scale_errors()

    if self.params.raw_data.errors_from_sample_residuals:
      self.errors_from_residuals()

  def _scale_all_parallel (self, file_names, db_mgr) :
    import multiprocessing
    import libtbx.introspection

    nproc = self.params.nproc
    if (nproc is None) or (nproc is Auto) :
      nproc = libtbx.introspection.number_of_processors()

    # Input files are supplied to the scaling processes on demand by
    # means of a queue.
    #
    # XXX The input queue may need to either allow non-blocking
    # put():s or run in a separate process to prevent the procedure
    # from blocking here if the list of file paths does not fit into
    # the queue's buffer.
    input_queue = multiprocessing.Manager().JoinableQueue()
    for file_name in file_names:
      print file_name
      input_queue.put(file_name)
    pool = multiprocessing.Pool(processes=nproc)
    # Each process accumulates its own statistics in serial, and the
    # grand total is eventually collected by the main process'
    # _add_all_frames() function.
    for i in xrange(nproc) :
      sm = scaling_manager(self.miller_set, self.i_model, self.params)
      pool.apply_async(
        func=sm,
        args=[input_queue, db_mgr],
        callback=self._add_all_frames)
    pool.close()
    pool.join()

    # Block until the input queue has been emptied.
    input_queue.join()

  def _scale_all_serial (self, tar_list, db_mgr) :
    """
    Scale frames sequentially (single-process).  The return value is
    picked up by the callback.
    """
    self.tar_to_scale_frame_adapter(tar_list, db_mgr)
    return (self)

  def __call__ (self, input_queue, db_mgr) :
    """ Scale frames sequentially within the current process.  The
    return value is picked up by the callback.  See also
    self.scale_all_serial()"""
    from Queue import Empty

    try :
      while True:
        try:
          file_name = input_queue.get_nowait()
        except Empty:
          return self
        self.tar_to_scale_frame_adapter(tar_list=[file_name,], db_mgr=db_mgr)
        input_queue.task_done()

    except Exception, e :
      print >> self.log, str(e)
      return None

cxi_merge.scaling_manager = scaling_manager

if (__name__ == "__main__"):
  result = run(args=sys.argv[1:])
  if result is None:
    sys.exit(1)

