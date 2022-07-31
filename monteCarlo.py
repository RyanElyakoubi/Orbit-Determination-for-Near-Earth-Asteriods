from math import sqrt, pi
import numpy as np
from OD4 import Obse as obs
from OD4 import emph
from statistics import stdev


def percentError(actual, calculated):
    percentError = abs((calculated - actual) / actual) * 100

    return percentError


obs.lon = "-81:24:52.9"  # longitude of obs.
obs.lat = "36:15:09.7"  # lat of obs.
obs.elev = 922  # elevation of obs, in meters (not a string!)
obs.date = "2020/07/05 07:08:00"  # (UTC date/time of observation)

line = "2003GE42,e,27.8543085,165.6469288,101.353156,2.63408,,0.375418379,8.4377868,07/05.29722/2020,2000,,"
asteroid = emph.readdb(line)
asteroid.compute(obs.date)
print(asteroid.a_ra, asteroid.a_dec)


def dist(x1, y1, x2, y2):
    dist = sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
    return dist


circle_count = 0
total_count = 0
center = np.array([1, 1])
list = []
total_trials = 100

for t in range(total_trials):
    for i in range(10000):
        x_coord = 2 * np.random.random()
        y_coord = 2 * np.random.random()

        if dist(x_coord, y_coord, center[0], center[1]) < 1:
            circle_count = circle_count + 1
        total_count = total_count + 1
    pi = (circle_count / total_count) * 4
    list.append(pi)

print("Approximate value of pi: ", np.mean(list),)
print("Standard error of approximation: ", stdev(list) / sqrt(total_trials))
print("Percent Error of approximation: ", percentError(pi, np.mean(list)), "%", sep="")
