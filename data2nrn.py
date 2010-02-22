from neuron import h
import numpy

def data2nrn(i, data, table):
    d = h.Vector()
    t = h.Vector()
    for j in range(0, len(table)) :
        pt = table[j]
        oitems = pt[1]
        for index in range(0, len(oitems)):
            if oitems[index] == i:
                t.append(pt[0][-1])
                x = data[j][index, 0]
                d.append(float(x))
    return (t, d)
