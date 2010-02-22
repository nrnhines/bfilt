import random
import math
import pylab
import wiener
import numpy
from neuron import h
cvode = h.CVode()
# for i in  range(len(yvec)):
#  print i, yvec[i], ydotvec[i]

def astate(state,time,noise)

def anoise(noise,time,state)

def hstate(state,time,noise)

def hnoise(noise,time,state)


class SDE(wiener.Wiener):
    
    a = -1.0
    b = 0.1
    Y = pylab.zeros((1,1)) 
    plotIndex = 0
    __Ydt = 1
    
    # Set Initial Condition
    def setic(self,nMax):
        yvec = h.Vector()
        cvode.states(yvec)
        self.Y = pylab.zeros((len(numpy.array(yvec)),nMax+1))
        self.Y[:,0] = numpy.array(yvec)
    
    # Test function
    def fos(self,t,y):
        return self.a*y
    
    # Real function
    def afun(self,t,y):
        yvec = h.Vector(y) 
        ydotvec = h.Vector()
        cvode.f(t, yvec, ydotvec)
        return numpy.array(ydotvec) 
    
    # Return Private Variable
    def getYdt(self):
        return self.__Ydt
    
    # Euler's Method
    def euler(self,dt=-1.):
        if dt < 0:
            dt = self.getdt()
        if (dt/self.getdt())%1==0.0:  # Wdt divides dt
            print 'Solving with time step, dt = "%s"' % dt
        else:
            print "Solving but Wdt does not divide dt..."
        tMax = self.gettMax()
        nMax = int(math.floor(tMax/dt))
        self.__Ydt = dt
        k = 0
        self.setic(nMax)
        while k < nMax:
            k = k + 1;
            self.Y[:,k] = self.Y[:,k-1] + (
            self.afun((k-1)*dt,self.Y[:,k-1])*dt +
            self.b*(self.evalW(k*dt) - self.evalW((k-1)*dt)))
    
    # Plot Solution Process
    def plotY(self):
        print "Plotting Y ... (Kill Plot Window to Return to Python)"
        tMax = self.gettMax()
        pylab.plot(pylab.arange(0,tMax+self.__Ydt,self.__Ydt),self.Y[self.plotIndex])
        pylab.show()

