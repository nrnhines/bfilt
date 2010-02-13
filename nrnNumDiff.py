import numdifftools as nd

def evalFun(nbf, p):
    saveParm = nbf.getParm()
    newParm = nbf.getParm()
    for i in range(len(p)):
        newParm.x[i] = p[i]
    nbf.setParm(newParm)
    L = nbf.likelihood()
    nbf.setParm(saveParm)
    return L

def Hessian(nbf):
    # Read current values of parameters (hopefully MLEs)
    parm = nbf.getParm()
    parmList = []
    for i in range(int(parm.size())):
        parmList.append(parm[i])
    # Create (likelihood) inline function
    LamFun = lambda p: evalFun(nbf,p)
    # Create Hessian (of likelihood) inline function
    HessFun = nd.Hessian(LamFun)
    # Evaluate Hessian and return
    return HessFun(parmList)

def Gradient(nbf):
    # Read current values of parameters (hopefully MLEs)
    parm = nbf.getParm()
    parmList = []
    for i in range(int(parm.size())):
        parmList.append(parm[i])
    # Create (likelihood) inline function
    LamFun = lambda p: evalFun(nbf,p)
    # Create Hessian (of likelihood) inline function
    GradFun = nd.Gradient(LamFun)
    # Evaluate Hessian and return
    return GradFun(parmList)

def LikePerturbed(nbf,perturbList,delta):
    saveParm = nbf.getParm()
    newParm = nbf.getParm()
    for i in range(len(perturbList)):
        newParm.x[i] = newParm[i] + delta*perturbList[i]
    nbf.setParm(newParm)
    L = nbf.likelihood()
    nbf.setParm(saveParm)
    return L