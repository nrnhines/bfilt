import numdifftools as nd

def evalFun(fun, args, p)
    # ######### TO_DO write setMRFparam
    setMRFparam(p)
    # ######### TO_DO check syntax
    fun(args)

def HessianMRF(fun,args):
    # Read current values of parameters (hopefully MLEs)
    # ######### TO_DO write getMRFparam
    param = getMRFparam();
    # Create (likelihood) inline function
    LamFun = lambda p: evalFun(fun,args,p)
    # Create Hessian (of likelihood) inline function
    HessFun = nd.Hessian(LamFun)
    # Evaluate Hessian and return
    return likeHess(param)
