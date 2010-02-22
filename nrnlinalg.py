from neuron import h
import numpy

#Determinant
def det(m):
    hm = h.Matrix(m.shape[0], m.shape[1])
    for i in range(int(hm.nrow())):
        hm.setcol(i, h.Vector(m[:,i]))
    e = h.ref(0)  
    d = hm.det(e)  
    return d*10.0**e[0]

def expm(m):
    hm = h.Matrix(m.shape[0], m.shape[1])
    for i in range(int(hm.nrow())):
        hm.setcol(i, h.Vector(m[:,i]))
    hm = hm.exp()
    mo = numpy.matrix(m)
    print m
    
    for i in range(int(hm.nrow())):
        for j in range(int(hm.nrow())):
            mo[i,j] = hm.getval(i, j)
    return mo
