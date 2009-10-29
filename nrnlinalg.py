from neuron import h

#Determinant
def det(m):
  hm = h.Matrix(m.shape[0], m.shape[1])
  for i in range(int(hm.nrow())):
    hm.setcol(i, h.Vector(m[:,i]))
  e = h.ref(0)  
  d = hm.det(e)  
  return d*10.0**e[0]
