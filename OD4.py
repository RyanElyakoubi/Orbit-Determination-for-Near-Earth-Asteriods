
import math
import numpy as np



def deg(x):
    x = float(x)
    x = math.degrees(x)
    return x % 360

def rad(x):
    x = float(x)
    return x * math.pi / 180



def solvekep(M, e):
    Eguess = M
    Eguess = Eguess - e * math.sin(Eguess)
    orginalE = M

    while abs(Eguess - orginalE) > 1e-4:
        prevE = Eguess
        Eguess = M + e * math.sin(Eguess)
    return Eguess

epsilon = rad(23.4358)
mu = .01720209895



def ephemeris(a, e, I, h, w, m0, t0, t, R):
    n = mu / math.sqrt(a ** 3)

    I = rad(I)
    h = rad(h)

    w = rad(w)
    m0 = rad(m0)
    m = m0 + n * (t - t0)
    E = solvekep(m, e)
    x = a * (math.cos(E) - e)

    y = a * (math.sqrt(1 - e ** 2) * math.sin(E))
    z = 0
    vecxyz = np.mat([[x], [y], [z]])

    r1 = np.mat([[math.cos(w), -math.sin(w), 0],
                 [math.sin(w), math.cos(w), 0],
                 [0, 0, 1]])

    r2 = np.mat([[1, 0, 0],
                 [0, math.cos(I), -math.sin(I)],
                 [0, math.sin(I), math.cos(I)]])

    r3 = np.mat([[math.cos(h), -math.sin(h), 0],
                 [math.sin(h), math.cos(h), 0],
                 [0, 0, 1]])

    vecXYZ = r3 * r2 * r1 * vecxyz
    vecSun = R

    vecRho = (vecSun + vecXYZ)

    vecRhoHat = vecRho / np.linalg.norm(vecRho)

    rX = vecRhoHat[0][0]
    rY = vecRhoHat[1][0] * math.cos(epsilon) - vecRhoHat[2][0] * math.sin(epsilon)
    rZ = vecRhoHat[1][0] * math.sin(epsilon) + vecRhoHat[2][0] * math.cos(epsilon)

    dFinal = math.asin(rZ)
    aFinal = math.atan2(rY, rX)
    return (deg(dFinal)), (deg(aFinal % (math.pi * 2)))