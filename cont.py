import numpy
import math
import copy

sqrteps = math.sqrt(numpy.finfo(numpy.double).eps)
maxiter = 20
tolf = 1e-5
tolx = numpy.inf
minstep = 1e-5

def jac(FUN,args,n,fv=None):
  # print 'jac args:', args
  # print 'jac n:', n
  # print 'jac FUN', FUN
  if fv == None:
    fv = FUN(*args)
    print 'j -1: fv', fv
  DF = numpy.matrix(numpy.zeros([1,len(args[n])]))  # nongeneral should be len(fv)
  for j in range(len(args[n])):
    temp = args[n][0,j]
    if abs(temp) > 1:
      h = sqrteps*abs(temp)
    else:
      h=sqrteps
    # print 'args[n]', args[n]
    x = args[n].copy()
    # if type(x) is numpy.core.defmatrix.matrix:
    # x[0,j] = temp + h
    # h = x[0,j] + h
    # else:
    x[0,j] = temp+h
    h = x[0,j] - temp  # trick to reduce finite precision error
    arglist = list(args)
    arglist[n] = x
    argtuple = tuple(arglist)
    print 'ARGTUPLE', argtuple
    f = FUN(*argtuple)
    print 'j', j, 'f', f
    # print argtuple #for debugging remove later
    DF[:,j] = (f - fv)/h
  return DF

def jac1(FUN,args,n,ic,fv=None):
  # print 'jac 1', [args[n][ic]]
  argsn = args[n]
  if type(argsn) is numpy.core.defmatrix.matrix:
    v = argsn[0,ic]
  elif type(argsn) is list:
    v = argsn[ic]
  else:
    v = argsn
  cargs = (numpy.matrix(v),FUN,args,n,ic)
  DF = jac(cFUN,cargs,0,fv)  # derivative wrt 0th comp of args
  return DF

def cFUN(aic,FUN,args,n,ic):
  print 'CFUN args', args
  print 'CFUN n', n
  print 'CFUN FUN', FUN
  print 'CFUN aic', aic

  argscopy = copy.deepcopy(args)
  if type(argscopy[n]) is numpy.core.defmatrix.matrix:
    temp = argscopy[n][0,ic]
  elif type(argscopy[n]) is list:
    temp = argscopy[n][ic]
  else:
    temp= argscopy
  if type(aic) is list:
    argscopy[n][ic] = aic[0]
  else:
    if type(argscopy[n]) is numpy.core.defmatrix.matrix:
      argscopy[n][0,ic] = aic
    else:
      argscopy[n][ic] = aic
  fx = FUN(*argscopy)

  return fx

def newt1(FUN, xx, args, n, ix):
  print 'newt1 FUN', FUN
  print 'newt1 xx', xx
  print 'newt1 args', args
  print 'newt1 n', n
  print 'newt1 ix', ix
  (x,fx,its) = newtn(cFUN, xx, (FUN, args, n, ix))
  return (x,fx,its)

def newtn(FUN, x, args): # returns (x, fx, its)
  dx = numpy.inf
  global maxiter, tolf, tolx

  for i in range(maxiter):
    fx = FUN(x,*args)
    normf = math.sqrt(fx*fx)  # nongeneral: should be dot product
    normdx = math.sqrt(dx*dx)  # nongeneral: should be dot product
    if dx == numpy.inf:
      print ' Iteration #', i, 'norm(F)', normf
    else:
      print ' Iteration #', i, 'norm(F)', normf, 'norm(dx)', normdx
    # print 'normf', normf, 'tolf', tolf, 'normdx', normdx, 'tolx', tolx
    # print (normf <= tolf and normdx <= tolx)
    if (normf <= tolf and normdx <= tolx):
      print ' Convergence!'
      return (x,fx,i)

    Df = jac(FUN,(x,)+args,0,fx)
    dx = -Df.I*fx
    x = x + dx
    # print 'Df', Df
    # print 'Df.I', Df.I
    # print 'fx', fx
    # print 'dx', dx

  i = numpy.inf
  print ' Maximum number of iterations exceeded'
  return (x,fx,i)

def cont(FUN, px, pa, args, ic, aicf, step):
  global minstep
  a0 = numpy.matrix(pa[-1])
  (x,fx,its)  = newtn(FUN, numpy.matrix(px[-1]), (a0,)+args)
  assert(its < numpy.inf)

  a = numpy.matrix(pa[-1])

  # print 'px', px
  # print 'pa', pa

  if aicf == a[0,ic]:
    return (px,pa)
  elif (aicf - a[0,ic])*step < 0:
    step = step * (-1)
  atend = False

  # print 'x', x
  DxF = jac(FUN,(x,a)+args,0)
  DaF = jac1(FUN,(x,a)+args,1,ic)
  # print 'DxF', DxF
  # print 'DaF', DaF
  DaX = -DxF.I*DaF
  atbegin = True

  while not atend:
    if not atbegin:
      print 'px[-1]', px[-1]
      print 'px[-2]', px[-2]
      print 'pa[-1][0,ic]', pa[-1][0,ic]
      print 'pa[-2][0,ic]', pa[-2][0,ic]
      DaX = (px[-1] - px[-2])/(pa[-1][0,ic] - pa[-2][0,ic])

    print 'DaX', DaX
    print 'pa[-1]', pa[-1]
    print 'px[-1]', px[-1]
    nt = math.sqrt((DaX*DaX.T)[0,0] + 1)
    tx = DaX/nt
    ta = 1/nt

    # print 'px & px[-1]', px, px[-1]
    # print 'pa & pa[-1]', pa, pa[-1]

    x0 = px[-1].copy()
    a0 = pa[-1].copy()
    its = numpy.inf

    while its == numpy.inf:
      x = x0 + step*tx
      a = a0.copy()
      a[0,ic] = a[0,ic] + step*ta
      print 'a', a
      print 'a0', a0
      print 'aicf', aicf
      if (aicf - a[0,ic])*step <= 0:
        atend = True
        laststep = step
        step = (aicf - a0[0,ic])/ta
        x = x0 + step*tx
        a[ic] = aicf

      (x,fx,its) = newtn(FUN, x, (a,)+args)

      print 'step =', step, 'its =', its, 'a(ic) = ', a[ic]

      if its == numpy.inf:
        atend = False

      if its > 10:
        step = step * 0.5
      elif its < 4:
        step = step * 1.3

      if numpy.abs(step) < minstep:
        print 'STEP less than MINSTEP; exiting ...'
        return (px,pa)

      if numpy.abs(step*ta) < minstep:
        print 'Increment to contin variable less than MINSTEP; exiting ...'
        return (px,pa)

    px.append(x.copy())
    pa.append(a.copy())
    atbegin = False
    # print 'add px', px
    # print 'add pa', pa

  return (px,pa)

def cont1(FUN, px, args, ix, ia, aicf, step):
  global minstep
  x0 = numpy.matrix(px[-1])
  (xix,fx,its)  = newt1(FUN, x0[0,ix], (x0,)+args, 0, ix)
  assert(its < numpy.inf)

  x = numpy.matrix(x0)
  x0[0,ix] = xix
  a = numpy.matrix(px[-1][0,ia])

  if aicf == a:
    return (px,pa)
  elif (aicf - a)*step < 0:
    step = step * (-1)
  atend = False

  # print 'x', x
  DxF = jac1(FUN,(x,)+args,0,ix)
  DaF = jac1(FUN,(x,)+args,0,ia)
  # print 'DxF', DxF
  # print 'DaF', DaF
  DaX = -DxF.I*DaF
  atbegin = True

  while not atend:
    if not atbegin:
      DaX = numpy.matrix(px[-1][0,ix] - px[-2][0,ix])/(px[-1][0,ia] - px[-2][0,ia])

    print 'DaX', DaX
    print 'px[-1]', px[-1]
    nt = math.sqrt((DaX*DaX.T)[0,0] + 1)
    tx = DaX/nt
    ta = 1/nt

    # print 'px & px[-1]', px, px[-1]
    # print 'pa & pa[-1]', pa, pa[-1]

    x0 = px[-1][0,ix]
    a0 = px[-1][0,ia]
    its = numpy.inf

    while its == numpy.inf:
      x = x0 + step*tx
      a = a0
      a = a + step*ta
      print 'a', a
      print 'a0', a0
      print 'aicf', aicf
      if (aicf - a)*step <= 0:
        atend = True
        laststep = step
        step = (aicf - a0)/ta
        x = x0 + step*tx
        a = aicf

      (x,fx,its) = newt1(FUN, x, (px[-1],)+args, 0, ix)

      print 'step =', step, 'its =', its

      if its == numpy.inf:
        atend = False

      if its > 10:
        step = step * 0.5
      elif its < 4:
        step = step * 1.3

      if numpy.abs(step) < minstep:
        print 'STEP less than MINSTEP; exiting ...'
        return px

      if numpy.abs(step*ta) < minstep:
        print 'Increment to contin variable less than MINSTEP; exiting ...'
        return px

    newx = px[-1].copy()
    newx[0,ix] = x[0,0]
    newx[0,ia] = a
    px.append(newx)
    atbegin = False
    # print 'add px', px
    # print 'add pa', pa

  return px

def xypoints(px,pa):
  xy = []
  for i in range(len(px)):
    xy.append([px[i][0,0],pa[i][0,0]])
  return xy

def testFUN(x,a,val):
  # print 'test x', x
  # print 'test a', a
  # print 'test val', val
  return x[0]*x[0] + a[0]*a[0] - val

def testFUN1(x,val):
  print 'testFUN1 x', x
  print 'testFUN1 v', val
  f = x[0,0]*x[0,0] + x[0,1]*x[0,1] - val
  print 'testFUN1 f', f
  return f

def test():
  global testFUN
  (px,pa) = cont(testFUN,[numpy.matrix([[1.0]])],[numpy.matrix([[1.0]])],(2,),0,-1.0,0.01)
  return (px, pa)

def test1():
  global testFUN
  px = cont1(testFUN1, [numpy.matrix([[1.0, 1.0]])], (2,), 0, 1, -1, 0.01)