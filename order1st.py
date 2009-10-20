import wiener

initial_Wdt = 0.001
initial_Ndt = 0.01
initial_Mdt = 0.1
initial_tMax = 5
initial_seed = 0
initial_jumpW = 0
initial_jumpM = 10000

class Order1st(wiener.Wiener):
    class Param:
        Ndt = 0.1
        Mdt = 0.01
        Wdt = 0.001
        tMax = 5
        a = 0.5
        sigp = 0.01
        sigm = 0.001
        
        
