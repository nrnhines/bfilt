import numpy
def ctest(a0, recurse=False):
    if recurse:
        pa = ctest(a0,False)
    else:
        a = numpy.matrix(a0)
        pa = []
        pa.append(a)
        pa.append(a)
    print 'pa & pa[-1] & pa[0] & pa[1]', pa, pa[-1],  pa[0], pa[1]
    return pa
