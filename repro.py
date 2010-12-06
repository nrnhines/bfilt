import copy
import time
import noisegen
from neuron import h
import os
import os.path

class Repro(object):
    def __init__(self, modelses, parm=None, channels=0, seed=0):
        # parm = None: saves parm to what's saved in session on first datagen
        self.modelses = modelses
        self.parm = parm
        self.channels = channels
        self.seed = seed
        self.timestamp = time.time()

    def retimestamp(self):
        self.timestamp = time.time()

    def datagen(self):
        ses_timestamp = os.path.getmtime(self.modelses)
        assert(ses_timestamp < self.timestamp)
        h.load_file(self.modelses)  # load session file
        mrflist = h.List("MulRunFitter")
        mrf = mrflist.o(int(mrflist.count())-1)
        N = mrf.p.pf.generatorlist.o(0).gen.po
        if self.parm == None:
            self.parm = N.getParm()
        else:
            N.setParm(self.parm)
        if self.channels == 0:
            G = noisegen.Gen(N)
            G.reseed(self.seed)
            Data = G.datasim()
            return Data
        else:
            print "Channel noise not implemented"