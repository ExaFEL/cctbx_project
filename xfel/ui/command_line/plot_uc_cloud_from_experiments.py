from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.plot_uc_cloud_from_experiments
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from libtbx.phil import parse

help_message = """
Plot a cloud of unit cell dimensions from stills. Provide either a combined_experiments.json
file or a specify individual .json files on the command line. To generate an overlay of
multiple plots (similar to grouping by run tag in the XFEL GUI), provide multiple
combined_experiments.json files named as ${tag}_combined_experiments.json and set
extract_tags to True in the phil scope.
"""

phil_str = """
  iqr_ratio = 1.5
    .type = float
    .help = Interquartile range multiplier for outlier rejection. Use None to disable outlier rejection.
  ranges = None
    .type = floats(6)
    .help = Lower and upper bounds for the ranges to display for each of the a, b and c axes
  extract_tags = False
    .type = bool
    .help = Extract tags from the names of multiple combined_experiments.json filenames and use
    .help = these tags to label multiple groups of experiments.
"""
phil_scope = parse(phil_str)

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # The script usage
    usage = "usage: %s [options] [param.phil] filenames" % libtbx.env.dispatcher_name

    self.tag = None
    self.reference_detector = None

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      read_experiments=True,
      check_format=False,
      epilog=help_message
      )

  def run(self):
    '''Execute the script.'''
    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True)
    def get_info(experiment):
      a, b, c, alpha, beta, gamma = experiment.crystal.get_unit_cell().parameters()
      return {'a':a,
              'b':b,
              'c':c,
              'alpha':alpha,
              'beta':beta,
              'gamma':gamma,
              'n_img':0}
    experiments_list = [e.data for e in params.input.experiments]
    if params.extract_tags:
      import os
      experiments_tags = [os.path.basename(f.filename).split("_combined_experiments.json")[0] for f in params.input.experiments]
      info_list = []
      for experiments in experiments_list:
        infos = []
        for experiment in experiments:
          infos.append(get_info(experiment))
        info_list.append(infos)
    else:
      experiments_tags = [""]
      infos = []
      for experiments in experiments_list:
        for experiment in experiments:
          infos.append(get_info(experiment))
      info_list = [infos]
    import xfel.ui.components.xfel_gui_plotter as pltr
    plotter = pltr.PopUpCharts()
    plotter.plot_uc_histogram(info_list=info_list, legend_list=experiments_tags, iqr_ratio = params.iqr_ratio, \
      ranges = params.ranges)
    plotter.plt.show()

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
