import numpy
import math

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
  DF = numpy.matrix(numpy.zeros([1,len(args[n])]))  # nongeneral should be len(fv)
  for j in range(len(args[n])):
    temp = args[n][j]
    if abs(temp) > 1:
      h = sqrteps*abs(temp)
    else:
      h=sqrteps
    # print 'args[n]', args[n]
    x = args[n].copy()
    x[j] = temp+h
    h = x[j] - temp  # trick to reduce finite precision error
    arglist = list(args)
    arglist[n] = x
    argtuple = tuple(arglist)
    f = FUN(*argtuple)
    # print argtuple #for debugging remove later
    DF[:,j] = (f - fv)/h
  return DF

def jac1(FUN,args,n,ic,fv=None):
  # print 'jac 1', [args[n][ic]]
  cargs = (numpy.matrix([[args[n][0,ic]]]),FUN,args,n,ic)
  DF = jac(cFUN,cargs,0,fv)  # derivative wrt 0th comp of args
  return DF

def cFUN(aic,FUN,args,n,ic):
  temp = args[n][ic]
  args[n][ic] = aic[0]
  # print 'cFUN temp', temp
  # print 'cFUN args', args
  # print 'cFUN aic', aic
  fx = FUN(*args)
  args[n][ic] = temp
  # print 'cFUN args[n][ic]', args[n][ic]
  # print 'cFUN fx', fx
  return fx

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


def cont(FUN, x0, a0, args, ic, aicf, step):
  global minstep
  (x,fx,its)  = newtn(FUN, numpy.matrix(x0), (a0,)+args)
  assert(its < numpy.inf)

  a = numpy.matrix(a0)
  px = [numpy.matrix(x)]
  pa = [numpy.matrix(a0)]

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
      DaX = (px[-1] - px[-2])/(pa[-1][ic] - pa[-2][ic])

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
    # print 'add px', px
    # print 'add pa', pa

  return (px,pa)

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

def test():
  global testFUN
  (px,pa) = cont(testFUN,[1.0],[1.0],(2,),0,0,0.01)
  return (px, pa)