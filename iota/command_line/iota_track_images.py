from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME iota.track_images
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

'''
Author      : Lyubimov, A.Y.
Created     : 07/21/2017
Last Changed: 07/21/2017
Description : IOTA image-tracking GUI launcher
'''

import wx
from iota import iota_version
from iota.components import iota_tracker as trk

class MainApp(wx.App):
  ''' App for the main GUI window  '''

  def OnInit(self):
    self.frame = trk.TrackerWindow(None, -1, title='IOTA IMAGE TRACKER v.{}'
                                               ''.format(iota_version))
    self.frame.SetMinSize(self.frame.GetEffectiveMinSize())
    self.frame.SetPosition((150, 150))
    self.frame.Show(True)
    self.frame.Layout()
    self.SetTopWindow(self.frame)
    return True

# ---------------------------------------------------------------------------- #

if __name__ == '__main__':
  app = MainApp(0)
  app.MainLoop()

