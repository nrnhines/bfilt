import copy
import time
import noisegen
from neuron import h
import os
import os.path

class Repro(object):
    def __init__(self, model, channels=0, seed=0, tvec=None, parm=None):
        # tvec == None: uses CollectionTimes
        # not tvec == None and channels == 0: asserts tvec = CollectionTimes
        # parm = None: saves parm to what's saved in session on first datagen
        self.model = model
        self.parm = parm
        self.tvec = tvec
        self.channels = channels
        self.seed = seed
        self.timestamp = time.time()

    def retimestamp(self):
        self.timestamp = time.time()

    def collectionTimes(self, N):
        if self.tvec == None:
            tvec = N.Eve.collectionTimes
        else:
            tvec = self.tvec
        return tvec

    def postProcess(self,Data, tvec):
        temp = 1

    def datagen(self):
        ses_timestamp = os.path.getmtime(self.model)
        assert(ses_timestamp < self.timestamp)
        if self.channels == 0:
            h.load_file(self.model)  # load session file for channels ==0
            mrflist = h.List("MulRunFitter")
            mrf = mrflist.o(int(mrflist.count())-1)
            N = mrf.p.pf.generatorlist.o(0).gen.po
            if self.parm == None:
                self.parm = N.getParm()
            else:  # these parms may have changed.  Reset.
                N.setParm(self.parm)
            self.collectionTimes(N)
            #else:  # If you have saved a tvec, make sure it is still valid
            #   assert(len(self.tvec) == len(N.Eve.CollectionTimes))
            #   for i in range(len(self.tvec)):
            #       assert(self.tvec[i] ==N.Eve.CollectionTimes[i])
            G = noisegen.Gen(N)
            G.reseed(self.seed)
            Data = G.datasim()
            return Data
        else:
            h.load_file(self.model)  # e.g. self.model = exper_channel.hoc
            h.load_file("exper_data.hoc")
            # m = h.ExperimentalModelContext(self.channels)
            mrflist = h.List("MulRunFitter")
            mrf = mrflist.o(int(mrflist.count())-1)
            N = mrf.p.pf.generatorlist.o(0).gen.po
            h.set_true_parameters()
            if self.parm == None:
                self.parm = h.true_parameters
            else:  # Saved params should match true params, unlike above
                assert(len(self.parm) == len(h.true_parameters))
                for i in range(len(self.parm)):
                    assert(self.parm[i] == h.true_parameters[i])
            self.tvec = self.collectionTimes(N)
            # Don't always want tvec=CollectionTimes (ie when filling fitness)
            Data = h.experimentalDataGenerator(self.channels, self.seed, self.tvec)
            self.postProcess(N, Data, tvec)
            return Data

class ReproFitness(Repro):
    def collectionTimes(self, N):
        ff = N.rf.fitnesslist.o(0)
        xmd = ff.xdat.c
        return xmd

    def postProcess(self, N, ymd, xmd):
        ff = N.rf.fitnesslist.o(0)
        sav = ff.n_chdat
        ff.set_data(xmd,ymd)
        h.execute("setxy()", ff)
        ff.n_chdat = sav  # Comment Out???

