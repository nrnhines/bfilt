from neuron import h
class Demo(object):
  def __init__(self):
    self.x = 2
    self.cbstate = False
    self.sbut = True
    self.ss = 'A variable label'
    self.build()
    self.map()

  def build(self):
    self.box = h.HBox()
    self.box.intercept(1)
    self.box.ref(self)
    h.xpanel("")
    h.xbutton("Button 1", (self.bact, 1))
    h.xbutton("Button 2", (self.bact, (2,)))
    h.xbutton("Button 3", self.bact_noarg)
    h.xbutton("Button 4", (self.bact_2arg, ("hello", 4)))
    for i in range(3):
      h.xradiobutton("Radio "+str(i), (self.bact, i))
    h.xmenu("Menu")
    for i in range(3):
      h.xbutton("Item "+str(i), (self.bact, i))
    for i in range(3):
      h.xradiobutton("Radio "+str(i), (self.bact, i))
    h.xcheckbox('checkbox', (self, 'cbstate'), self.cb)
    h.xstatebutton('state', (self, 'sbut'), self.st)
    h.xmenu()
    h.xpanel()
    self.g = h.Graph()
    self.g.menu_tool('graph menu tool 1', self.gcb, (self.gsel, 1))
    h.xpanel("")
    h.xvalue("x", (self, "x"), 1, self.chgx)
    h.xcheckbox('checkbox', (self, 'cbstate'), self.cb)
    h.xstatebutton('state', (self, 'sbut'), self.st)
    h.xlabel('fixed label')
    h.xvarlabel((self, 'ss'))
    h.xslider((self, 'x'), 0, 20, self.slide)
    self.g.menu_tool('graph menu tool 2', self.gcb, (self.gsel, 2))
    h.xpanel()
    self.box.intercept(0)

  def map(self):
    self.box.map("Demo")

  def bact(self, i):
    print "bact", i

  def bact_noarg(self):
    print "bact_noarg"

  def bact_2arg(self, a1, a2):
    print "bact_2arg", a1, a2

  def chgx(self):
    print "x =", self.x
    self.ss = 'x is ' + str(self.x)

  def cb(self):
    print "cbstate=", self.cbstate

  def st(self):
    print 'sbut=', self.sbut

  def slide(self):
    print 'slide x=', self.x

  def gcb(self, mode, x, y, key):
    print mode, x, y, key

  def gsel(self, i):
    print 'menu tool', i, 'selected'

Demo()
