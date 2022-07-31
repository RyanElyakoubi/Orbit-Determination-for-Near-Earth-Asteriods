
from MethodOfGauss import *
from FinalODFunctions import *
from OD2MOD import *

#Constants
k = .01720209895
mu = 1
epsilon = rad(-23.4352)
t_0 = JD(2018,7,22,6,0,0)

def mog(alphas,deltas,t,R_x,R_y,R_z,t_0):
    #get the sun vectors and their norms
    R = get_R_vectorsf(R_x,R_y,R_z)
    R_1 = R[0]
    R_2 = R[1]
    R_3 = R[2]

    R_norms = get_R_normsf(R)

    R_1_norm = R_norms[0]
    R_2_norm = R_norms[1]
    R_3_norm = R_norms[2]

    #find rho_hats:
    rho_hats = get_rho_hats(alphas,deltas)

    rho_1_hat = rho_hats[0]
    rho_2_hat = rho_hats[1]
    rho_3_hat = rho_hats[2]

    #find D constants:
    D_constants = get_D(rho_hats,R)
    D_0 = get_D(rho_hats,R)[0]
    D_11,D_12,D_13 = get_D(rho_hats,R)[1],get_D(rho_hats,R)[2],get_D(rho_hats,R)[3]
    D_21,D_22,D_23 = get_D(rho_hats,R)[4],get_D(rho_hats,R)[5],get_D(rho_hats,R)[6]
    D_31,D_32,D_33 = get_D(rho_hats,R)[7],get_D(rho_hats,R)[8],get_D(rho_hats,R)[9]

    #find Gaussian time intervals:
    t_1 = t[0]
    t_2 = t[1]
    t_3 = t[2]
    taus = get_taus(t)
    tau_1,tau_3,tau = get_taus(t)[0],get_taus(t)[1],get_taus(t)[2]

    #find initial r_2 value using scalar equation of Lagrange:
    A_1 = tau_3 / tau
    B_1 = (A_1 / 6) * (tau**2 - tau_3**2)
    A_3 = -tau_1 / tau
    B_3 = (A_3 / 6) * (tau**2 - tau_1**2)
    A = (A_1 * D_21 - D_22 + A_3 * D_23) / (-D_0)
    B = (B_1 * D_21 + B_3 * D_23) / (-D_0)
    E = -2*(np.dot(rho_2_hat, R_2))
    F = R_2_norm**2
    a = -(A**2 + A*E + F)
    b = -mu*(2*A*B + B*E)
    c = -mu**2 * B**2

    #solve scalar equation of lagrange to find all roots
    r_2_norms = poly.polyroots([c,0,0,b,0,0,a,0,1])
    r_2_norms.tolist()

    #Check all roots which are real, positive, and result in positive rho mags
    for i in range(len(r_2_norms)):
        #Calculate rho norms for each root (this means preliminarily calculating f and g series and c constants)
        f_1 = 1 - (mu * tau_1**2) / (2 * r_2_norms[i]**3)
        f_3 = 1 - (mu * tau_3**2) / (2 * r_2_norms[i]**3)
        f = f_1,f_3

        g_1 = tau_1 - (mu * tau_1**3) / (6 * r_2_norms[i]**3)
        g_3 = tau_3 - (mu * tau_3**3) / (6 * r_2_norms[i]**3)
        g = g_1,g_3

        c_1 = get_c(f,g)[0]
        c_2 = get_c(f,g)[1]
        c_3 = get_c(f,g)[2]

        rho_1_norm = (c_1 * D_11 + c_2 * D_12 + c_3 * D_13) / (c_1 * D_0)
        rho_2_norm = (c_1 * D_21 + c_2 * D_22 + c_3 * D_23) / (c_2 * D_0)
        rho_3_norm = (c_1 * D_31 + c_2 * D_32 + c_3 * D_33) / (c_3 * D_0)


        r_2_norms_valid = []
        #Check if root is positive, real, and returns positive rho norms
        if np.isreal(r_2_norms[i]) and r_2_norms[i] > 0 and rho_1_norm > 0 and rho_2_norm > 0 and rho_3_norm > 0:
            r_2_norms_valid.append(np.real(r_2_norms[i]))
            break

    #loop through all valid roots
    for r_2_norm in r_2_norms_valid:
        #find c_1, c_3 using very truncated f, g series:
        f_1 = 1 - (mu * tau_1**2) / (2 * r_2_norm**3)
        f_2 = 1
        f_3 = 1 - (mu * tau_3**2) / (2 * r_2_norm**3)
        g_1 = tau_1 - (mu * tau_1**3) / (6 * r_2_norm**3)
        g_2 = 0
        g_3 = tau_3 - (mu * tau_3**3) / (6 * r_2_norm**3)
        c_1 = g_3 / (f_1 * g_3 - g_1 * f_3)
        c_2 = -1
        c_3 = -g_1 / (f_1 * g_3 - g_1 * f_3)

        #find rho_1_norm,rho_2_norm,rho_3_norm:
        rho_1_norm = (c_1 * D_11 + c_2 * D_12 + c_3 * D_13) / (c_1 * D_0)
        rho_2_norm = (c_1 * D_21 + c_2 * D_22 + c_3 * D_23) / (c_2 * D_0)
        rho_3_norm = (c_1 * D_31 + c_2 * D_32 + c_3 * D_33) / (c_3 * D_0)
        rho_norms = rho_1_norm,rho_2_norm,rho_3_norm

        #Get initial r vectors
        r_1 = get_r(rho_norms,rho_hats,R)[0]
        r_2 = get_r(rho_norms,rho_hats,R)[1]
        r_3 = get_r(rho_norms,rho_hats,R)[2]

        #Get initial velocity vector
        d_1 = -f_3 / (f_1 * g_3 - f_3 * g_1)
        d_3 = f_1 / (f_1 * g_3 - f_3 * g_1)
        r_2_dot = d_1 * r_1 + d_3 * r_3
        r_2_norm_new = np.linalg.norm(r_2)

        #Iterate through to converge to final r_2,r_2_norm
        while abs(r_2_norm - r_2_norm_new) > 10**-12:
            #Set the old position magnitude to new
            r_2_norm = r_2_norm_new

            #Calculate f,g series to the fourth order
            u = mu / r_2_norm**3
            z = np.dot(r_2,r_2_dot) / r_2_norm**2
            q = np.dot(r_2_dot,r_2_dot) / r_2_norm**2 - u

            f_1 = 1 - (1/2)*u*tau_1**2 + (1/2)*u*z*tau_1**3 + (1/24)*(3*u*q - 15*u*z**2 + u**2)*tau_1**4
            f_3 = 1 - (1/2)*u*tau_3**2 + (1/2)*u*z*tau_3**3 + (1/24)*(3*u*q - 15*u*z**2 + u**2)*tau_3**4
            f = f_1,f_3

            g_1 = tau_1 - (1/6)*u*tau_1**3 + (1/4)*u*z*tau_1**4
            g_3 = tau_3 - (1/6)*u*tau_3**3 + (1/4)*u*z*tau_3**4
            g = g_1,g_3

            c_1 = get_c(f,g)[0]
            c_2 = get_c(f,g)[1]
            c_3 = get_c(f,g)[2]


            #Calculate the new rho mags
            rho_1_norm = (c_1 * D_11 + c_2 * D_12 + c_3 * D_13) / (c_1 * D_0)
            rho_2_norm = (c_1 * D_21 + c_2 * D_22 + c_3 * D_23) / (c_2 * D_0)
            rho_3_norm = (c_1 * D_31 + c_2 * D_32 + c_3 * D_33) / (c_3 * D_0)

            #Calculate the new r vectors
            r_1 = rho_1_norm * rho_1_hat - R_1
            r_2 = rho_2_norm * rho_2_hat - R_2
            r_3 = rho_3_norm * rho_3_hat - R_3

            #Calculate the new velocity vector
            d_1 = -f_3 / (f_1 * g_3 - f_3 * g_1)
            d_3 = f_1 / (f_1 * g_3 - f_3 * g_1)
            r_2_dot = d_1 * r_1 + d_3 * r_3

            r_2_norm = r_2_norm_new
            r_2_norm_new = np.linalg.norm(r_2)




            #Light travel correction
            t_1 = t[0] - rho_1_norm / light
            t_2 = t[1] - rho_2_norm / light
            t_3 = t[2] - rho_3_norm / light
            tau_3 = k * (t_3 - t_2)
            tau_1 = k * (t_1 - t_2)
            tau = tau_3 - tau_1

        #Final position and velocity vectors after iteration
        r_2 = np.dot(epsilonrot,r_2)
        r_2_dot = np.dot(epsilonrot,r_2_dot)

        #Get orbital elements
        a,e,I,Omega,omega,M = getOrbitaElements(r_2,r_2_dot,t_2,t_0)

        return a,e,I,Omega,omega,M