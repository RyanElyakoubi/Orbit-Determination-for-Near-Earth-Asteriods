import math

def solvekep(M,e):
    Eguess = M
    Eguess = Eguess - e*math.sin(Eguess)
    prevE = M
    while abs(Eguess - prevE) > 1e-4:
        prevE = Eguess
        Eguess = M + e*math.sin(Eguess)
    return Eguess