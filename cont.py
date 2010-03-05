import numpy
import math

sqrteps = math.sqrt(numpy.finfo(numpy.double).eps)
maxiter = 20
tolf = 1e-10
tolx = numpy.inf
minstep = 1e-5

def jac(FUN,args,n,fv=None):
  if fv == None:
    fv = FUN(*args)
  DF = numpy.zeros([1,len(args[n])])  # nongeneral should be len(fv)
  for j in range(len(args[n])):
    temp = args[n][j]
    if abs(temp) > 1:
      h = sqrteps*abs(temp)
    else:
      h=sqrteps
    x = args[n][:]
    x[j] = temp+h
    h = x[j] - temp  # trick to reduce finite precision error
    arglist = list(args)
    arglist[n] = x
    argtuple = tuple(arglist)
    f = FUN(*argtuple)
    print argtuple #for debugging remove later
    DF[:,j] = (f - fv)/h
  return DF

def jac1(FUN,args,n,ic,fv):

def cFUN(aic,FUN,args,n,ic):
  temp = args[n][ic]
  args[n][ic] = aic
  fx = FUN(*args)
  args[n][ic] = temp
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
    if (normf <= tolf and normdx <= tolx):
      print ' Convergence!'
      return (x,fx,i)

    Df = jac(FUN,x,args,fx)
    dx = - Df.I*fx
    x = x + dx

  i = numpy.inf
  print ' Maximum number of iterations exceeded'


def cont(FUN, x0, a0, args, ic, aicf, step):
  global minstep
  (x,fx,its)  = newtn(FUN, x0, (a0,)+args)
  assert(its < numpy.inf)

  a = a0
  px = [x]
  pa = [a0]

  if aicf == a[ic]:
    return (px,pa)
  elif (aicf - a[ic])*step < 0:
    step = step * (-1)
  atend = False

  print 'x', x
  DxF = jac(FUN,(x,a)+args,0)
  DaF = jac1(FUN,(x,a)+args,1,ic)
  DaX = -DxF.I*DaF
  atbegin = True

  while not atend:
    if not atbegin:
      DaX = (px[-1] - px[-2])/(pa[-1][ic] - pa[-2][ic])

    nt = math.sqrt((DaX*DaX.T)[0,0] + 1)
    tx = DaX/nt
    ta = 1/nt

    x0 = px[-1]
    a0 = pa[-1][ic]
    its = numpy.inf

    while its == numpy.inf:
      x = x0 + step*tx
      a[ic] = a0 + step*ta

      if (aicf - a[ic])*step <= 0:
        atend = True
        laststep = step
        step = (aicf - a0)/ta
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

    px.append(x)
    pa.append(a)

  return (px,pa)

def testFUN(x,a,val):
  print 'test x', x
  print 'test a', a
  print 'test val', val
  return x[0]*x[0] + a[0]*a[0] - val

def test():
  global testFUN
  (px,pa) = cont(testFUN,[1.0],[1.0],(2,),0,0,0.01)
  return (px, pa)