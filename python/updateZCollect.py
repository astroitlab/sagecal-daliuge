import os,logging
from dfms.drop import BarrierAppDROP
import struct

logger = logging.getLogger(__name__)

class UpdateZCollectApp(BarrierAppDROP):

    def initialize(self, **kwargs):
        super(UpdateZCollectApp, self).initialize(**kwargs)

    def run(self):
        if len(self.inputs) < 1:
            raise Exception("at least one input is expected by this application")
        if len(self.outputs) != 1:
            raise Exception("at least and only one output is expected by this application")

        with open(self.outputs[0].path, 'w') as f:
            for d in self.inputs:
                with open(d.path, 'r') as df:
                    tagId = struct.unpack("i", df.read(4))
                    if tagId == 0:
                        f.writelines([df.path,'\n'])