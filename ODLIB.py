from tkinter import W
import numpy as np
import math as m

k = .01720209895


def get_vectors(filename):
    with open(filename) as data:
        lines = data.readlines()

        rx = float(lines[0])
        ry = float(lines[1])
        rz = float(lines[2])
        vx = float(lines[3])
        vy = float(lines[4])
        vz = float(lines[5])

    return rx, ry, rz, vx, vy, vz


rx = get_vectors("ElyakoubiInput.txt")[0]
ry = get_vectors("ElyakoubiInput.txt")[1]
rz = get_vectors("ElyakoubiInput.txt")[2]
vx = get_vectors("ElyakoubiInput.txt")[3] / k
vy = get_vectors("ElyakoubiInput.txt")[4] / k
vz = get_vectors("ElyakoubiInput.txt")[5] / k


def quadrant(cos_input, sin_input):
    if cos_input > 1 or cos_input < -1 or sin_input > 1 or sin_input < -1:
        return ('invalid input')
    if round(cos_input * cos_input + sin_input * sin_input, 1) != 1:
        return ('invalid input')
    if cos_input > 0 and sin_input > 0:
        if m.degrees(m.asin(sin_input)) >= 0 and m.degrees(m.asin(sin_input)) <= 90:
            return (m.degrees(m.asin(sin_input)))
        else:
            return m.degrees((m.acos(cos_input)))
    elif cos_input < 0 and sin_input > 0:
        if m.degrees(m.asin(sin_input)) >= 90 and m.degrees(m.asin(sin_input)) <= 180:
            return m.degrees((m.asin(sin_input)))
        else:
            return m.degrees((m.acos(cos_input)))
    elif cos_input < 0 and sin_input < 0:
        # Read about this once
        return (180 - m.degrees(m.asin(sin_input)))
    elif cos_input > 0 and sin_input < 0:
        # Read about this once
        return (360 + m.degrees(m.acos(cos_input)))
    # return(alt_deg,azi_deg1,azi_deg2)


def magnitude(vector):
    scalar = m.sqrt((vector[0] ** 2) + (vector[1] ** 2) + (vector[2] ** 2))
    return scalar


r = magnitude(np.array([rx, ry, rz]))
v = magnitude(np.array([vx, vy, vz]))
rvec = (np.array([rx, ry, rz]))
vvec = (np.array([vx, vy, vz]))


def angular_momentum(filename):
    with open(filename) as data:
        lines = data.readlines()

        position = np.array([float(lines[0]), float(lines[1]), float(lines[2])])
        velocity = np.array([float(lines[3]), float(lines[4]), float(lines[5])])

        h = np.cross(position, velocity)
        h = h * (1 / k)

    return h


h = angular_momentum("BradleyInput.txt")


def semi_major_axis(filename):
    with open(filename) as data:
        lines = data.readlines()

        position = np.array([float(lines[0]), float(lines[1]), float(lines[2])])
        velocity = (np.array([float(lines[3]), float(lines[4]), float(lines[5])])) / k

        r_s = magnitude(position)
        v_s = magnitude(velocity)

        a = (1) / ((2 / r_s) - (v_s ** 2))

    return a


a = float(semi_major_axis("BradleyInput.txt"))


def eccentricity(h, a):
    h = magnitude(h)
    e = m.sqrt(1 - ((h ** 2)) / (a))

    return e


e = eccentricity(h, a)


def inclination_angle(h):
    hz = abs(h[2])
    h = magnitude(h)
    i = m.degrees(m.acos(hz / h))

    return i


i = inclination_angle(h)


def longitude_of_the_ascending_node(h, i):
    hy = h[1]
    hx = h[0]
    h = magnitude(h)

    sininput = (hx / (h * m.sin((m.radians(i)))))
    cosinput = (-hy / (h * m.sin((m.radians(i)))))
    omega = (quadrant(cosinput, sininput))
    return omega


omega = longitude_of_the_ascending_node(h, i)


def u(omega, r, rx, ry, rz, i):
    cosu = ((rx * m.cos(m.radians(omega))) + (ry * m.sin(m.radians(omega)))) / r
    sinu = (rz) / (r * m.sin(m.radians(i)))
    u = quadrant(cosu, sinu)
    return u


u = u(omega, r, rx, ry, rz, i)


def nu(a, e, h, rvec, vvec, r):
    posvecdot = np.dot(rvec, vvec)
    sinnu = (((a * (1 - e ** 2) * posvecdot / r)) / (e * magnitude(h)))
    cosnu = ((a * (1 - e ** 2)) / (r) - 1) / (e)
    nu = quadrant(cosnu, sinnu)

    return nu


nu = nu(a, e, h, rvec, vvec, r)


def normalize(var):
    while var < 0:
        var = var + 360
    while var > 360:
        var = var - 360
    return var


def argument_of_perihelion(nu, u):
    w = u - nu
    w = normalize(w)
    return w


w = argument_of_perihelion(nu, u)


def eccentric_anomaly(a, r, e):
    E = m.acos((1 - (r / a)) / e)
    return E


E = eccentric_anomaly(a, r, e)


def mean_anomaly(E, e):
    M = m.degrees(E - e * m.sin(E))
    return M


julian_date = 2458313.500000000
M = mean_anomaly(E, e)


def time_of_perihelion_passage(julian_date, M, k, a):
    T = (julian_date - (m.radians(M) / (m.sqrt(k ** 2 / a ** 3))))
    return T


def percent_error(actual, calculated):
    percent_error = abs((calculated - actual) / actual) * 100

    return percent_error


T = time_of_perihelion_passage(julian_date, M, k, a)


def find_orbital_elements(filename):
    a = semi_major_axis(filename)
    e = eccentricity(h, a)
    i = inclination_angle(h)
    omega = longitude_of_the_ascending_node(h, i)
    w = argument_of_perihelion(nu, u)
    M = mean_anomaly(E, e)

    return a, e, i, omega, w, M
