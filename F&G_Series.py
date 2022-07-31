from cmath import sqrt, tau
from math import cos, sin
import numpy as np
from ODLIB import percent_error, magnitude


mu = 1.000000000
c = 3 * 10 ** 8
k = .01720209895

tau1 = -0.3261857571141891
tau3 = 0.05084081855693949
r2 = [0.26662393644794813, -1.381475976476564, -0.5048589337503169]
r2dot = [0.8442117090940343, -0.39728396707075087, 0.14202728258915864]


def f1_g1(tau1, r2, r2dot, flag):
    u = mu / np.linalg.norm(r2) ** 3
    z = np.dot(r2, r2dot) / np.linalg.norm(r2) ** 2
    q = np.dot(r2dot, r2dot) / np.linalg.norm(r2) ** 2 - u

    if flag == 4:
        f1 = 1 - ((1 / 2) * u * tau1 ** 2) + (1 / 2) * u * z * tau1 ** 3 + (1 / 24) * (
                    3 * u * q - 15 * u * z ** 2 + u ** 2) * tau1 ** 4
        g1 = tau1 - (1 / 6) * u * tau1 ** 3 + (1 / 4) * u * z * tau1 ** 4
        return f1, g1
    if flag == 3:
        f1 = 1 - ((1 / 2) * u * tau1 ** 2) + (1 / 2) * u * z * tau1 ** 3
        g1 = tau1 - (1 / 6) * u * tau1 ** 3
        return f1, g1
    return


def f3_g3(tau3, r2, r2dot, flag):
    u = mu / np.linalg.norm(r2) ** 3
    z = np.dot(r2, r2dot) / np.linalg.norm(r2) ** 2
    q = np.dot(r2dot, r2dot) / np.linalg.norm(r2) ** 2 - u

    if flag == 4:
        f3 = 1 - ((1 / 2) * u * tau3 ** 2) + (1 / 2) * u * z * tau3 ** 3 + (1 / 24) * (
                    3 * u * q - 15 * u * z ** 2 + u ** 2) * tau3 ** 4
        g3 = tau3 - (1 / 6) * u * tau3 ** 3 + (1 / 4) * u * z * tau3 ** 4
        return f3, g3
    elif flag == 3:
        f3 = 1 - ((1 / 2) * u * tau3 ** 2) + (1 / 2) * u * z * tau3 ** 3
        g3 = tau3 - (1 / 6) * u * tau3 ** 3
        return f3, g3
    return


f1 = f1_g1(tau1, r2, r2dot, 4)[0]
g1 = f1_g1(tau1, r2, r2dot, 4)[1]
f3 = f3_g3(tau3, r2, r2dot, 4)[0]
g3 = f3_g3(tau3, r2, r2dot, 4)[1]

f1test = 0.9823149707782799
f3test = 0.99961917185922
g1test = -0.32418770657924106
g3test = 0.05083441832100904

print(percent_error(f1test, f1))
print(percent_error(g1test, g1))
print(percent_error(f3test, f3))
print(percent_error(g3test, g3))

# c1 = g3/(f1 * g3 - g1 * f3)
# c3 = -g3/(f1 * g3 - g1 * f3)

# r1 = f1 * r2 + g1 * r2dot
# r3 = f3 * r2 + g3 * r2dot

# rho1 = R1 + r1
# rho2 = R2 + r2
# rho3 = R3 + r3

# def light_correction(rho, og_time, c):
#     time = og_time - rho/c
#     return time

h = np.cross(r2, r2dot)
r_s = magnitude(r2)
v_s = magnitude(r2dot)

a = (1) / ((2 / r_s - v_s ** 2))
e = (sqrt(1 - (magnitude(h) ** 2) / a)).real

n = magnitude(h) / (a ** 2 * sqrt(1 - e ** 2)).real


def deltaE(r2, r2dot, tau):
    num = (((np.dot(r2, r2dot)) / (n * a ** 2)) * (cos(n * tau - np.dot(r2, r2dot) / (n * a ** 2)))) + (
                (1 - magnitude(r2) / a) * sin(n * tau - ((np.dot(r2, r2dot)) / (n * a ** 2))))
    if num > 0:
        sign = 0
    elif num < 0:
        sign = 1
    if sign == 0:
        deltaE = n * tau + .85 * e - (np.dot(r2, r2dot)) / (n * a ** 2)
    elif sign == 1:
        deltaE = n * tau - .85 * e - ((np.dot(r2, r2dot)) / (n * a ** 2))

    return deltaE


def g_f_functions(deltaE, taui, tau):
    deltaE = deltaE(r2, r2dot, tau)
    oldDelta = 0
    while abs(deltaE - oldDelta) > 10 ** -12:
        fequation = deltaE - (1 - magnitude(r2) / a) * sin(deltaE) + ((np.dot(r2, r2dot)) / (n * a ** 2)) * (
                    1 - cos(deltaE)) - n * taui
        fprime = 1 - (1 - magnitude(r2) / a) * cos(deltaE) + ((np.dot(r2, r2dot)) / (n * a ** 2)) * sin(deltaE)
        oldDelta = deltaE
        deltaE = oldDelta - fequation / fprime

    gi = taui + (sin(deltaE) - deltaE) / n
    fi = 1 - (a / magnitude(r2)) * (1 - cos(deltaE))
    return gi, fi


g1 = g_f_functions(deltaE, tau1, tau)[0]
f1 = g_f_functions(deltaE, tau1, tau)[1]
g3 = g_f_functions(deltaE, tau3, tau)[0]
f3 = g_f_functions(deltaE, tau3, tau)[1]
print(g1, f1, g3, f3)
