


from OD2MOD import *
from FinalODFunctions import *
import OD2MOD
import FinalODFunctions

# CONSTANT VALUES
k = .01720209895
mu = 1
epsilon = rad(-23.4352)
t_0 = FinalODFunctions.JD(2021, 6, 25, 18, 0, 0 )


def mog(t_0, alphas, deltas, t, R_x, R_y, R_z):
    # SUN VECTORS
    R = FinalODFunctions.getRVectorsF(R_x,R_y,R_z)
    R_1 = R[0]
    R_2 = R[1]
    R_3 = R[2]

    R_norms = FinalODFunctions.getRNormsF(R)
    R_1_norm = R_norms[0]
    R_2_norm = R_norms[1]
    R_3_norm = R_norms[2]

    # RHO HATS
    rho_hats = FinalODFunctions.getRhoHats(alphas, deltas)

    rho_1_hat = rho_hats[0]
    rho_2_hat = rho_hats[1]
    rho_3_hat = rho_hats[2]

    D_constants = FinalODFunctions.getD(rho_hats,R)
    D_0 = FinalODFunctions.getD(rho_hats,R)[0]
    D_11,D_12,D_13 = FinalODFunctions.getD(rho_hats,R)[1],FinalODFunctions.getD(rho_hats,R)[2],FinalODFunctions.getD(rho_hats,R)[3]
    D_21,D_22,D_23 = FinalODFunctions.getD(rho_hats,R)[4],FinalODFunctions.getD(rho_hats,R)[5],FinalODFunctions.getD(rho_hats,R)[6]
    D_31,D_32,D_33 = FinalODFunctions.getD(rho_hats,R)[7],FinalODFunctions.getD(rho_hats,R)[8],FinalODFunctions.getD(rho_hats,R)[9]

    # GAUSSIAN
    t_1 = t[0]
    t_2 = t[1]
    t_3 = t[2]
    taus = FinalODFunctions.get_taus(t)
    tau_1,tau_3,tau = FinalODFunctions.get_taus(t)[0],FinalODFunctions.get_taus(t)[1],FinalODFunctions.get_taus(t)[2]


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

    # LAGRANGE
    r_2_norms = poly.polyroots([c,0,0,b,0,0,a,0,1])
    r_2_norms.tolist()

 # F&G Series
    for i in range(len(r_2_norms)):

        f_1 = 1 - (mu * tau_1**2) / (2 * r_2_norms[i]**3)
        f_3 = 1 - (mu * tau_3**2) / (2 * r_2_norms[i]**3)
        f = f_1,f_3

        g_1 = tau_1 - (mu * tau_1**3) / (6 * r_2_norms[i]**3)
        g_3 = tau_3 - (mu * tau_3**3) / (6 * r_2_norms[i]**3)
        g = g_1,g_3

        c_1 = getC(f,g)[0]
        c_2 = getC(f,g)[1]
        c_3 = getC(f,g)[2]

        rho_1_norm = (c_1 * D_11 + c_2 * D_12 + c_3 * D_13) / (c_1 * D_0)
        rho_2_norm = (c_1 * D_21 + c_2 * D_22 + c_3 * D_23) / (c_2 * D_0)
        rho_3_norm = (c_1 * D_31 + c_2 * D_32 + c_3 * D_33) / (c_3 * D_0)


        r_2_norms_valid = []
        # VALIDATION
        if np.isreal(r_2_norms[i]) and r_2_norms[i] > 0 and rho_1_norm > 0 and rho_2_norm > 0 and rho_3_norm > 0:
            r_2_norms_valid.append(np.real(r_2_norms[i]))
            break


    for r_2_norm in r_2_norms_valid:
        # C1 & C3
        f_1 = 1 - (mu * tau_1**2) / (2 * r_2_norm**3)
        f_2 = 1
        f_3 = 1 - (mu * tau_3**2) / (2 * r_2_norm**3)
        g_1 = tau_1 - (mu * tau_1**3) / (6 * r_2_norm**3)
        g_2 = 0
        g_3 = tau_3 - (mu * tau_3**3) / (6 * r_2_norm**3)
        c_1 = g_3 / (f_1 * g_3 - g_1 * f_3)
        c_2 = -1
        c_3 = -g_1 / (f_1 * g_3 - g_1 * f_3)

        #find RHO'S (3)
        rho_1_norm = (c_1 * D_11 + c_2 * D_12 + c_3 * D_13) / (c_1 * D_0)
        rho_2_norm = (c_1 * D_21 + c_2 * D_22 + c_3 * D_23) / (c_2 * D_0)
        rho_3_norm = (c_1 * D_31 + c_2 * D_32 + c_3 * D_33) / (c_3 * D_0)
        rho_norms = [rho_1_norm, rho_2_norm, rho_3_norm]

        # INITIAL R VECTORS
        r_1 = FinalODFunctions.getR(rho_norms,rho_hats,R)[0]
        r_2 = FinalODFunctions.getR(rho_norms,rho_hats,R)[1]
        r_3 = FinalODFunctions.getR(rho_norms,rho_hats,R)[2]

        # INITIAL VELOCITY VECTOR
        d_1 = -f_3 / (f_1 * g_3 - f_3 * g_1)
        d_3 = f_1 / (f_1 * g_3 - f_3 * g_1)
        r_2_dot = d_1 * r_1 + d_3 * r_3
        r_2_norm_new = np.linalg.norm(r_2)

        # FINAL R2 & R2 NORM
        while abs(r_2_norm - r_2_norm_new) > 10**-12:

            r_2_norm = r_2_norm_new

            # F&G SERIES 4TH ORDER
            u = mu / r_2_norm**3
            z = np.dot(r_2,r_2_dot) / r_2_norm**2
            q = np.dot(r_2_dot,r_2_dot) / r_2_norm**2 - u

            f_1 = 1 - (1/2)*u*tau_1**2 + (1/2)*u*z*tau_1**3 + (1/24)*(3*u*q - 15*u*z**2 + u**2)*tau_1**4
            f_3 = 1 - (1/2)*u*tau_3**2 + (1/2)*u*z*tau_3**3 + (1/24)*(3*u*q - 15*u*z**2 + u**2)*tau_3**4
            f = f_1,f_3

            g_1 = tau_1 - (1/6)*u*tau_1**3 + (1/4)*u*z*tau_1**4
            g_3 = tau_3 - (1/6)*u*tau_3**3 + (1/4)*u*z*tau_3**4
            g = g_1,g_3

            c_1 = getC(f,g)[0]
            c_2 = getC(f,g)[1]
            c_3 = getC(f,g)[2]

            # NEW RHO MAGS
            rho_1_norm = (c_1 * D_11 + c_2 * D_12 + c_3 * D_13) / (c_1 * D_0)
            rho_2_norm = (c_1 * D_21 + c_2 * D_22 + c_3 * D_23) / (c_2 * D_0)
            rho_3_norm = (c_1 * D_31 + c_2 * D_32 + c_3 * D_33) / (c_3 * D_0)

            # NEW R VECTORS
            r_1 = rho_1_norm * rho_1_hat - R_1
            r_2 = rho_2_norm * rho_2_hat - R_2
            r_3 = rho_3_norm * rho_3_hat - R_3

            # NEW VELOCITY VECTOR
            d_1 = -f_3 / (f_1 * g_3 - f_3 * g_1)
            d_3 = f_1 / (f_1 * g_3 - f_3 * g_1)
            r_2_dot = d_1 * r_1 + d_3 * r_3

            r_2_norm = r_2_norm_new
            r_2_norm_new = np.linalg.norm(r_2)

            # LIGHT TRAVEL CORRECTION
            t_1 = t[0] - rho_1_norm / light
            t_2 = t[1] - rho_2_norm / light
            t_3 = t[2] - rho_3_norm / light
            tau_3 = k * (t_3 - t_2)
            tau_1 = k * (t_1 - t_2)
            tau = tau_3 - tau_1

        # Final POSITION & VELOCITY
        r_2 = np.dot(OD2MOD.epsilonrot,r_2)
        r_2_dot = np.dot(OD2MOD.epsilonrot,r_2_dot)

        # ORBITAL ELEMENT$

        a,e,I,Omega,omega,M = OD2MOD.getOrbitalElements(r_2,r_2_dot,t_2,t_0)


        return a,e,I,Omega,omega,M


#print(mog(t_0,2459390.500000000, 276.28516667, -17.44480556,[-0.206885, 0.913277, 0.39588], [0.368894, 0.869094, 0.376727], [-0.520476, 0.8004, 0.346949]))

print(mog(t_0,[-0.206885, 0.913277, 0.39588],[-0.520476, 0.8004, 0.346949], -17.44480556,[-0.206885, 0.913277, 0.39588], [0.368894, 0.869094, 0.376727], [-0.520476, 0.8004, 0.346949]))
