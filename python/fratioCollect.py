import os,logging
<<<<<<< HEAD
from dfms.drop import BarrierAppDROP
=======
from dlg.drop import BarrierAppDROP
>>>>>>> 7e1b3c310aa6385f4401327fe70bad900190208e
import struct

logger = logging.getLogger(__name__)

class FratioCollectApp(BarrierAppDROP):

    def initialize(self, **kwargs):
        super(FratioCollectApp, self).initialize(**kwargs)

    def run(self):
        if len(self.inputs) < 1:
            raise Exception("at least one input is expected by this application")
        if len(self.outputs) != 1:
            raise Exception("at least and only one output is expected by this application")

        with open(self.outputs[0].path, 'w') as f:
            for d in self.inputs:
                if os.path.exists(d.path):
                    with open(d.path, 'r') as df:
                        tagId = struct.unpack("i", df.read(4))[0]
                        if tagId == 0:
                            f.writelines([d.path,'\n'])
