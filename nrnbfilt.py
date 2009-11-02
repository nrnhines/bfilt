from neuron import h

class NrnBFilt(object):
  def __init__(self, ho):
    self.rf = ho

  def likelihood(self):
    print 'likelihood', self, self.rf
    for var in self.rf.yvarlist:
      print var
    for data in self.rf.fitnesslist:
      print data
    return 1.0
