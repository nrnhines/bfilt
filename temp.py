from runTheFit import *
M.Obs.C[0].sigma = 0
P.B = numpy.matrix([0,0,0])
M.change(P)
Data = M.sim()

from data2nrn import data2nrn
p = data2nrn(0, Data, M.FitEvents)
print p

from neuron import h
g = h.Graph()
p[1].line(g, p[0])
g.exec_menu("View = plot")

