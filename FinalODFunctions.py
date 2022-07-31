import numpy as np
from math import sqrt, cos, sin, radians as rad
import numpy.polynomial.polynomial as poly

mu = 1
k = .01720209895
light = 173.144633
#Functions
#converts a single civil date to JD
def JD(years,months,days,hours,minutes,seconds):
    J_0 = 367 * years - int((7/4)*(years + int((months+9)/12))) + int(275 * months / 9) + days + 1721013.5
    UT = hours / 24 + minutes / (24 * 60) + seconds / (24 * 60 * 60)
    JD = J_0 + UT
    return JD

#converts a list of civil dates to JD
def JDf(years,months,days,hours,minutes,seconds):
    JD = []
    for i in range(3):
        J_0 = 367 * years[i] - int((7/4)*(years[i] + int((months[i]+9)/12))) + int(275 * months[i] / 9) + days[i] + 1721013.5
        UT = hours[i] / 24 + minutes[i] / (24 * 60) + seconds[i] / (24 * 60 * 60)
        JD.append(J_0 + UT)
    return JD

#converts RA in hms to radians
def alpha_radf(alpha_hr,alpha_min,alpha_sec):
    alpha_hr = alpha_hr + alpha_min / 60 + alpha_sec / 3600
    alpha_deg = alpha_hr * 15
    alpha_rad = rad(alpha_deg)
    return alpha_rad

#converts dec in dms to radians
def delta_radf(delta_deg,delta_min,delta_sec):
    delta_hr = delta_deg + delta_min / 60 + delta_sec / 3600
    delta_rad = rad(delta_hr)
    return delta_rad

#puts sun vectors into a list
def get_R_vectorsf(R_x,R_y,R_z):
    R_vectors = []
    for i in range(3):
        R = np.array([R_x[i],R_y[i],R_z[i]])
        R_vectors.append(R)
    return R_vectors

#puts sun vector mags into a list
def get_R_normsf(R_vectors):
    R_norms = []
    for i in range(3):
        R_norms.append(np.linalg.norm(R_vectors[i]))
    return R_norms

#converts lists of RA in hms to radians
def get_alphasf(alpha_hours,alpha_minutes,alpha_seconds):
    alphas = []
    for i in range(3):
        alphas.append(alpha_radf(alpha_hours[i],alpha_minutes[i],alpha_seconds[i]))
    return alphas

#covnerts lits of dec in dms to radians
def get_deltasf(delta_deg,delta_minutes,delta_seconds):
    deltas = []
    for i in range(3):
        deltas.append(delta_radf(delta_deg[i],delta_minutes[i],delta_seconds[i]))
    return deltas

#creates list of rho hat vectors from the RA's and decs
def get_rho_hats(alphas,deltas):
    rho_hats = []
    for i in range(3):
        rho_hat = np.array([cos(alphas[i])*cos(deltas[i]),sin(alphas[i])*cos(deltas[i]),sin(deltas[i])])
        rho_hats.append(rho_hat)
    return rho_hats

#D constants
def get_D(rho_hats,R_vectors):
    D_0 = np.dot(rho_hats[0],np.cross(rho_hats[1],rho_hats[2]))
    D_1i = []
    for i in range(3):
        D_1i.append(np.dot(np.cross(R_vectors[i],rho_hats[1]),rho_hats[2]))
    D_2i = []
    for i in range(3):
        D_2i.append(np.dot(np.cross(rho_hats[0],R_vectors[i]),rho_hats[2]))
    D_3i = []
    for i in range(3):
        D_3i.append(np.dot(rho_hats[0],np.cross(rho_hats[1],R_vectors[i])))
    return D_0,D_1i[0],D_1i[1],D_1i[2],D_2i[0],D_2i[1],D_2i[2],D_3i[0],D_3i[1],D_3i[2],

#finds gaussian time intervals given julian times
def get_taus(t):
    t_1,t_2,t_3 = t[0],t[1],t[2]
    tau_3 = k * (t_3 - t_2)
    tau_1 = k * (t_1 - t_2)
    tau = tau_3 - tau_1
    return tau_1,tau_3,tau


#Gets c constants from f and g series
def get_c(f,g):
    f_1 = f[0]
    f_3 = f[1]
    g_1 = g[0]
    g_3 = g[1]
    c_1 = g_3 / (f_1 * g_3 - g_1 * f_3)
    c_2 = -1
    c_3 = -g_1 / (f_1 * g_3 - g_1 * f_3)
    return c_1,c_2,c_3

#Gets r vectors given R and rho vectors
def get_r(rho_norms,rho_hats,R):
    r = []
    for i in range(3):
        r.append(rho_norms[i]*rho_hats[i] - R[i])
    return r

def arcsec_to_rad(arcsecs):
    rads = []
    for i in range(3):
        rads.append((4.84814e-6) * arcsecs[i])
    return rads
