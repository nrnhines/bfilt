import numpy

bounds = []
bounds.append([numpy.matrix([[0,1,0,0]]),0])
bounds.append([numpy.matrix([[0,-1,0,0]]),-1])
bounds.append([numpy.matrix([[0,0,1,0]]),0])
bounds.append([numpy.matrix([[0,0,-1,0]]),-1])
bounds.append([numpy.matrix([[0,0,0,1]]),0])
bounds.append([numpy.matrix([[0,0,0,-1]]),-1])
