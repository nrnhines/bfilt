import random
import math
import pylab

class Wiener():
  # "Private" Variables -- need to be changed together
  Proc = random.Random()
  Meas = random.Random()  
  __seed = 0
  __jumpW = 0
  __jumpM = 10000
  __Wdt = 0.001
  __tMax = 0
  __changed = True 
  
  # "Public" Variables
  W = pylab.zeros(1)  # Don't change this yourself
  M = pylab.zeros(1)  # This either
  recalc = False      # When True, recalcs W when parameters change
                      # When False, just notes that parameters have changed
  savefname = '/home/sean/Desktop/PyPlots/Wplot.pdf'

  # Functions to set private variable
  #   then either recalcuate W or note parameters have changed
  
  def setseed(self,s): # seed parameter for random generator 
    self.__seed = s
    self.__changed = True
    if self.recalc:
      self.calcMW() 
      self.__changed = False
    print self.__seed
 
  def setjump(self,j,m): # jump parameter for random generator
    self.__jumpW = j
    self.__jumpM = m
    self.__changed = True
    if self.recalc:
      self.calcMW() 
    print self.__jumpW
    print self.__jumpM
 
  def setWdt(self,Wdt):
    self.__Wdt = Wdt
    self.__changed = True
    if self.recalc:
      self.calcMW() 
    print self.__Wdt

  def settMax(self,tMax):
    self.__tMax = tMax
    self.__changed = True
    if self.recalc:
      self.calcMW() 
    print self.__tMax

  # Functions to return values of private variables

  def getseed(self):
    return self.__seed

  def getjumpW(self):
    return self.__jumpW

  def getjumpM(self):
    return self.__jumpM

  def getWdt(self):
    return self.__Wdt

  def gettMax(self):
    return self.__tMax

  def hasChanged(self):
    return self.__changed

  # Calculate Wiener Process
  def calcMW(self):
    print "Recalculating..."
    self.Proc.seed(self.__seed)
    self.Proc.jumpahead(self.__jumpW)
    self.Meas.seed(self.__seed)
    self.Meas.jumpahead(self.__jumpM)
    nMax = int(math.floor(self.__tMax/self.__Wdt))
    self.W = pylab.zeros(nMax+1)
    self.M = pylab.zeros(nMax+1)
    k = 0
    while k < nMax:
      k = k + 1;
      self.W[k]=self.W[k-1]+self.Proc.normalvariate(0,1)*math.sqrt(self.__Wdt)
      self.M[k]=self.Meas.normalvariate(0,1)
    self.__changed = False
    print "...Done!"

  # Plot Wiener Process
  def plotW(self):
    if self.__changed:
      self.calcMW()   
    print "Plotting ... (Kill Plot Window to Return to Python)"
    pylab.plot(pylab.arange(0,self.__tMax+self.__Wdt,self.__Wdt),self.W)
    pylab.show()
    # pylab.savefig(self.savefname)

  # Evaluate Wiener Process at a specified time
  def evalW(self,t):
    if self.__changed:
      self.calcMW()   
    if t > self.__tMax:
      t = self.__tMax
      print "Warning: t > tMax"
    return self.W[int(round(t/self.__Wdt))]

  # Evaluate Measurement Noise at a specified time
  def evalM(self,t):
    if self.__changed:
      self.calcMW()   
    if t > self.__tMax:
      t = self.__tMax
      print "Warning: t > tMax"
    return self.M[int(round(t/self.__Wdt))]

