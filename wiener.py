import random
import math
import pylab

class Wiener(random.Random):
  # "Private" Variables -- need to be changed together
  __seed = 0
  __jump = 0
  __dt = 0.001
  __tMax = 0
  __changed = True 
  
  # "Public" Variables
  W = pylab.zeros(1)  # Don't change this yourself
  recalc = False      # When True, recalcs W when parameters change
                      # When False, just notes that parameters have changed

  # Functions to set private variable
  #   then either recalcuate W or note parameters have changed
  
  def setseed(self,s): # seed parameter for random generator 
    self.__seed = s
    self.__changed = True
    if self.recalc:
      self.calcW() 
      self.__changed = False
    print self.__seed
 
  def setjump(self,j): # jump parameter for random generator
    self.__jump = j
    self.__changed = True
    if self.recalc:
      self.calcW() 
    print self.__jump
 
  def setdt(self,dt):
    self.__dt = dt
    self.__changed = True
    if self.recalc:
      self.calcW() 
    print self.__dt

  def settMax(self,tMax):
    self.__tMax = tMax
    self.__changed = True
    if self.recalc:
      self.calcW() 
    print self.__tMax

  # Functions to return values of private variables

  def getseed(self):
    return self.__seed

  def getjump(self):
    return self.__jump

  def getdt(self):
    return self.__dt

  def gettMax(self):
    return self.__tMax

  def hasChanged(self):
    return self.__changed


  # Calculate Wiener Process
  def calcW(self):
    print "Recalculating..."
    self.seed(self.__seed)
    self.jumpahead(self.__jump)
    nMax = int(math.floor(self.__tMax/self.__dt))
    self.W = pylab.zeros(nMax+1)
    k = 0
    while k < nMax:
      k = k + 1;
      self.W[k] = self.W[k-1] + self.normalvariate(0,1)*math.sqrt(self.__dt)
    self.__changed = False
    print "...Done!"

  # Plot Wiener Process
  def plotW(self):
    if self.__changed:
      self.calcW()   
    print "Plotting ... (Kill Plot Window to Return to Python)"
    pylab.plot(pylab.arange(0,self.__tMax+self.__dt,self.__dt),self.W)
    pylab.show()

  # Evaluate Wiener Process a specified time
  def evalW(self,t):
    if self.__changed:
      self.calcW()   
    if t > self.__tMax:
      t = self.__tMax
      print "Warning: t > tMax"
    return self.W[int(round(t/self.__dt))]
