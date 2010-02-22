import math
import numpy

sqrtEps = math.sqrt(numpy.finfo(numpy.double).eps

def jac3(fun,value,x,y,z):
    DFx = pylab.zeros((len(value),len(x)))
    k = 0
    while k<len(x):
        temp = x[k]
        if abs(temp) > 1
            h = sqrtEps*abs(temp)
        else
            h = sqrtEps
        x[k] = temp + h
        h = x[k] - temp  # trick to reduce finite precision error
        df = fun(x,y,z)
        x[k] = temp
        DFx[:,k] = (df - value)/h;
        k += 1
    return DFx
