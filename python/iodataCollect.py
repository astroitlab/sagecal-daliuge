import os,logging
from dlg.drop import BarrierAppDROP

logger = logging.getLogger(__name__)

class IODataCollectApp(BarrierAppDROP):

    def initialize(self, **kwargs):
        super(IODataCollectApp, self).initialize(**kwargs)

    def run(self):
        if len(self.inputs) < 1:
            raise Exception("at least one input is expected by this application")
        if len(self.outputs) != 1:
            raise Exception("at least and only one output is expected by this application")

        with open(self.outputs[0].path, 'w') as f:
            f.writelines([fd.path for fd in self.inputs])
