from neuron import h
import testcr
testcr.cvodewrap.fs.use_fixed_step = 1.0
testcr.parrun(nchannels=100)
