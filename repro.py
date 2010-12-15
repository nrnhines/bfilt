import copy
import time
import noisegen
from neuron import h
import os
import os.path

class Repro(object):
    def __init__(self, model, parm=None, channels=0, seed=0):
        # parm = None: saves parm to what's saved in session on first datagen
        self.model = model
        self.parm = parm
        self.channels = channels
        self.seed = seed
        self.timestamp = time.time()

    def retimestamp(self):
        self.timestamp = time.time()

    def datagen(self):
        ses_timestamp = os.path.getmtime(self.model)
        assert(ses_timestamp < self.timestamp)
        if self.channels == 0:
            h.load_file(self.model)  # load session file
            mrflist = h.List("MulRunFitter")
            mrf = mrflist.o(int(mrflist.count())-1)
            N = mrf.p.pf.generatorlist.o(0).gen.po
            if self.parm == None:
                self.parm = N.getParm()
            else:
                N.setParm(self.parm)
            G = noisegen.Gen(N)
            G.reseed(self.seed)
            Data = G.datasim()
            return Data
        else:
            h.load_file(self.model)  # e.g. self.model = exper_channel.hoc
            h.load_file("exper_data.hoc")
            h.set_true_parameters()
            self.parm = h.true_parameters()
            Data = h.experimentalDataGenerator(self.channels, self.seed)
            return Data