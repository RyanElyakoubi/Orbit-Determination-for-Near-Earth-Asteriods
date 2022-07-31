from ODLIB import magnitude, quadrant, find_orbital_elements, percent_error
from math import cos, sin, sqrt, asin, radians, degrees
import numpy as np

import NewtonRaphson

Rvec = np.array([float(-6.573688002856684E-01), float(7.092600133268118E-01), float(3.074357332141535E-01)])

k = .01720209894
mu = 1.00000000
epsilon = radians(23.439281)

a = (find_orbital_elements("ElyakoubiInput.txt")[0])
e = (find_orbital_elements("ElyakoubiInput.txt")[1])
i = radians(find_orbital_elements("ElyakoubiInput.txt")[2])
omega = radians(find_orbital_elements("ElyakoubiInput.txt")[3])
w = radians(find_orbital_elements("ElyakoubiInput.txt")[4])

t = 2458333.500000000
T = 2458158.719533933792

M = k*(sqrt(mu/(a**3)))*(t-T)

E = NewtonRaphson.solvekep(M,e)

r = np.array([a*cos(E)-(a*e), a*sqrt(1-e**2)*sin(E), 0])

eclipticmatrix1 = np.array([[cos(omega), -sin(omega), 0],
                            [sin(omega), cos(omega), 0],
                            [0, 0, 1]])

eclipticmatrix2 = np.array([[1, 0, 0],
                            [0, cos(i), -sin(i)],
                            [0, sin(i), cos(i)]])

eclipticmatrix3 = np.array([[cos(w), -sin(w), 0],
                            [sin(w), cos(w), 0],
                            [0, 0, 1]])

recliptic = (eclipticmatrix1 @ eclipticmatrix2 @ eclipticmatrix3 @ r)

equatorialmatrix =  np.array([[1, 0, 0],
                              [0, cos(epsilon), -sin(epsilon)],
                              [0, sin(epsilon), cos(epsilon)]])


requatorial = equatorialmatrix @ recliptic

rho_vec = requatorial + Rvec

rho_unit_vec = rho_vec/magnitude(rho_vec)

DEC = (asin((rho_unit_vec[2])))
RA = (quadrant((rho_unit_vec[0])/cos(DEC), (((rho_unit_vec[1])/cos(DEC)))))
DEC = degrees(DEC)

print(DEC, RA)

print("Declination Percent Error:", percent_error(float(31.87410), DEC), "%")
print(" Right Ascension Percent Error:", percent_error(float(265.58799), RA), "%")
