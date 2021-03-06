from neuron import h
import numpy
import scipy.stats as stats
import pickle

def read(fname):
  f = open(fname, "r")
  n = pickle.load(f)
  print 'n=',n
  rl = []
  for i in range(n):
    try:
      r = pickle.load(f)
      rl.append(r)
    except:
      break
  print '#items read = ', len(rl)
  f.close()
  return rl

def csvfile(fname, differDOF=None):
  rl = read(fname+'.dat')
  f = open(fname+'.csv','w')
  f.write('Seed, Time, pVal, ml, otml\n')
  for r in rl:
    assert(r[0][3] >= 3)
    v1 = str(r[0][1])
    v2 = str(r[0][0])
    v3 = str(chisqprob(r, differDOF))
    v4 = str(r[3][1])
    v5 = str(r[1][1])
    c = ', '
    n = '\n'
    f.write(str(v1+c+v2+c+v3+c+v4+c+v5+n))

def chisqprob(r, differDOF=None):
  r3 = r[3]
  if differDOF == None:
    n = len(r3[0])
  else:
    n = differDOF
  ml = r3[1]
  otml = r[1][1]
  cs = 2.0*(otml - ml)
  pval = stats.chisqprob(cs, n)
  return pval

def within95(r, differDOF=None):
  alpha = 0.05
  pval = chisqprob(r, differDOF)
  return pval >= alpha

def n_within(rl, differDOF=None):
  if rl[0][3] < 3 :
    print "file does not contain data needed by n_within"
    return
  i = 0
  for r in rl:
    if (within95(r, differDOF)):
      i += 1
  return i

def plot_parms(rl):
  if rl[0][3] < 3 :
    print "file does not contain data needed by plot_parms"
    return
  g = h.Graph(0)
  g.view(2)
  g.size(0,10,0,10)
  g.label(.5,.9,'#channels '+str(rl[0][1]),2,1,.5, 0, 1)
  for r in rl:
    parm = r[3][0]
    g.mark(parm[0], parm[1], "O", 5, 1, 1)
  return g

def doit(fname):
  rl = read(fname)
  g = plot_parms(rl)
  n = len(rl)
  percent = 100*n_within(rl)/n
  print "%d%% of %d instances of %d channels have (2,4) within 95%% confidence region"%(percent,n,rl[0][0][2])
  return g

if 1 and __name__ == '__main__':
  a = []
  #a.append(doit('results_100.512'))
  #a.append(doit('results_1000.512'))
  #a.append(doit('results.fixed.1.128'))
  a.append(doit('test256.dat'))
