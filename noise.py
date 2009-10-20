import random
import math

class NoiseParams:
    def __init__(self):
        self.dt = 1
        self.tstop = 5
        self.seed = 0
    def inc(self):
        self.seed += 1
    def clock(self):
        self.seed = int(math.floor(100*time.time()))

# ZeroNoise: Noise sequence of random variable, all zero with certainty.
# ZeroNoise is the parent class to the other noise classes.
class ZeroNoise:
    # Noise classes instantiated with two arguments: 
    #    P (parameters; of class NoiseParams)
    #    I  (a unique positive integer to make seeds unique)
    # where
    #    P should be the same for all noise processes in a model
    #    I should be different for each noise process in a model
    def __init__(self, p, i):
        self.R = random.Random()
        # save local copies of P and I
        self.P = p
        self.I = i
        # calculate noise
        self.X = self.calc()
    
    # just change P not I
    def change(self, p):
        self.__init__(p, self.I)
        
    # Return Zero process
    def clear(self):
        # set seed for new noise generation
        self.R.seed(self.P.seed + 1e6*self.I)
        # how many numbers, given dt and tstop?
        kstop = 1+int(math.floor(self.P.tstop/self.P.dt))
        return [0.0]*kstop
        
    # Calculate (zero) process: subclasses will overload this function
    def calc(self):
        return self.clear()
        
    # Eavaluate at a specified time
    def eval(self,t):
        if t > self.P.tstop:
            t = self.P.tstop
            print "Warning evaluating Noise at time > tstop"
        return self.X[int(round(t/self.P.dt))]
        
# Gauss: a sequence of standard Gaussian random variables
# Gauss is suitable for measurement noise
class Gauss(ZeroNoise):
    def calc(self):
        X = self.clear() # seed & initialize list for efficient processing
        k = 0
        while k < len(X):
            X[k]=self.R.normalvariate(0,1)
            k += 1; # increment k 2nd because X[0] is Gaussian
        return X

# Wiener: discrete approximation to a Wiener process
# Wiener is suitable for process noise
class Wiener(ZeroNoise):
    def calc(self):
        X = self.clear() # seed & initialize list for efficient processing
        k = 0
        while k<len(X)-1:
            k += 1;  # increment k 1st because X[0] = 0, by def'n
            X[k]=X[k-1]+self.R.normalvariate(0,1)*math.sqrt(self.P.dt)
        return X

# Other vector classes are subclasses of ZeroVector
class ZeroVector:
    def __init__(self, p, i0, d):
        self.P = p
        self.I = i0
        self.D = d
        self.C = []
        k = 0
        while k < d:
            self.C.append(self.scalar(p, i0+k))
            k += 1
            
    # Function scalar overloaded in subclasses
    def scalar(self, p, i):
        return ZeroNoise(p, i)
        
    # Change parameters
    def change(self, p):
        i0 = self.I
        d = self.D
        self.__init__(p, i0, d)
        
    # Evaluate the vector at a given time
    def eval(self, t):
        k=0
        E = []
        while k<self.D:
            E.append(self.C[k].eval(t))
            k += 1
        return E

class WienerVector(ZeroVector):
    def scalar(self, p, i):
        return Wiener(p, i)
        
class GaussVector(ZeroVector):
    def scalar(self, p, i):
        return Gauss(p, i)
        
# Noise class gives you:
# a vector of Wiener (process noise), and
# a vector of Gauss (measurement noise)
class Noise:
    def __init__(self, p, pdim, mdim):
        self.P = p
        self.PDim = pdim
        self.MDim = mdim
        self.W = WienerVector(p, 0, pdim)  #jumpahead: 0..pdim-1
        self.G = GaussVector(p, pdim, mdim)  #: pdim..pdim+mdim-1
    
    # Change parameters
    def change(self, p):
        pdim = self.PDim
        mdim = self.MDim
        self.__init__(p, pdim, mdim)

