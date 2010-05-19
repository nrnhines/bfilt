from neuron import h
class TestGUI(object):
    def __init__(self):
        self.Vg = False
        self.mg = False
        self.ng = False
        self.hg = False
        self.Cag = False
        self.Og = False
        self.C1g = False
        self.C2g = False
        self.Vl = False
        self.ml = False
        self.nl = False
        self.hl = False
        self.Cal = False
        self.Ol = False
        self.C1l = False
        self.C2l = False
        self.Vs0 = False
        self.ms0 = False
        self.ns0 = False
        self.hs0 = False
        self.Cas0 = False
        self.Os0 = False
        self.C1s0 = False
        self.C2s0 = False
        self.Vs1 = False
        self.ms1 = False
        self.ns1 = False
        self.hs1 = False
        self.Cas1 = False
        self.Os1 = False
        self.C1s1 = False
        self.C2s1 = False
        self.ineq = True
        self.build()
        self.map()

    def build(self):
        self.box = h.HBox()
        self.box.intercept(1)
        self.box.ref(self)
        h.xpanel("")
        h.xlabel('0<=')
        h.xcheckbox('',(self,'Vg'),self.button)
        h.xcheckbox('',(self,'mg'),self.button)
        h.xcheckbox('',(self,'ng'),self.button)
        h.xcheckbox('',(self,'hg'),self.button)
        h.xcheckbox('',(self,'Cag'),self.button)
        h.xcheckbox('',(self,'Og'),self.button)
        h.xcheckbox('',(self,'C1g'),self.button)
        h.xcheckbox('',(self,'C2g'),self.button)
        h.xpanel()
        h.xpanel("")
        h.xlabel('<=1')
        h.xcheckbox('V',(self,'Vl'),self.button)
        h.xcheckbox('m',(self,'ml'),self.button)
        h.xcheckbox('n',(self,'nl'),self.button)
        h.xcheckbox('h',(self,'hl'),self.button)
        h.xcheckbox('[Ca2+]',(self,'Cal'),self.button)
        h.xcheckbox('O',(self,'Ol'),self.button)
        h.xcheckbox('C1',(self,'C1l'),self.button)
        h.xcheckbox('C2',(self,'C2l'),self.button)
        h.xpanel()
        h.xpanel("")
        h.xlabel('S0')
        h.xcheckbox('',(self,'Vs0'),self.button)
        h.xcheckbox('',(self,'ms0'),self.button)
        h.xcheckbox('',(self,'ns0'),self.button)
        h.xcheckbox('',(self,'hs0'),self.button)
        h.xcheckbox('',(self,'Cas0'),self.button)
        h.xcheckbox('',(self,'Os0'),self.button)
        h.xcheckbox('',(self,'C1s0'),self.button)
        h.xcheckbox('',(self,'C2s0'),self.button)
        h.xpanel()
        h.xpanel("")
        h.xlabel('S1')
        h.xcheckbox('V',(self,'Vs1'),self.button)
        h.xcheckbox('m',(self,'ms1'),self.button)
        h.xcheckbox('n',(self,'ns1'),self.button)
        h.xcheckbox('h',(self,'hs1'),self.button)
        h.xcheckbox('[Ca2+]',(self,'Cas1'),self.button)
        h.xcheckbox('O',(self,'Os1'),self.button)
        h.xcheckbox('C1',(self,'C1s1'),self.button)
        h.xcheckbox('C2',(self,'C2s1'),self.button)
        h.xpanel()
        h.xpanel("")
        h.xbutton("Add S2", self.button)
        h.xbutton("Remove Empty", self.button)
        h.xbutton("Close", self.button)
        h.xpanel()
        self.box.intercept(0)

    def map(self):
        self.box.map("TestGUI")

    def button(self):
        print 'button pressed'

TestGUI()