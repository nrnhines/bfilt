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

def within95(r):
  alpha = 0.05
  n = len(r[2])
  ml = r[5]
  otml = r[3]
  cs = 2.0*(otml - ml)
  pval = stats.chisqprob(cs, n)
  return pval >= alpha

def n_within(rl):
  i = 0
  for r in rl:
    if (within95(r)):
      i += 1
  return i

def plot_parms(rl):
  g = h.Graph(0)
  g.view(2)
  g.size(0,10,0,10)
  g.label(.5,.9,'#channels '+str(rl[0][1][1]),2,1,.5, 0, 1)
  for r in rl:
    parm = r[4]
    g.mark(parm[0], parm[1], "O", 5, 1, 1)
  return g

def doit(fname):
  rl = read(fname)
  g = plot_parms(rl)
  n = len(rl)
  percent = 100*n_within(rl)/n
  print "%d%% of %d instances of %d channels have (2,4) within 95%% confidence region"%(percent,n,rl[0][1][1])
  return g

if __name__ == '__main__':
  a = []
  a.append(doit('results_100.512'))
  a.append(doit('results_1000.512'))
  a.append(doit('results.fixed.1.128'))
