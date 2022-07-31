
from ODLIB import find_orbital_elements, percent_error

a = find_orbital_elements("ElyakoubiInput.txt")[0]
e = find_orbital_elements("ElyakoubiInput.txt")[1]
i = find_orbital_elements("ElyakoubiInput.txt")[2]
omega = find_orbital_elements("ElyakoubiInput.txt")[3]
w = find_orbital_elements("ElyakoubiInput.txt")[4]

print("a:", a)
print("e:", e)
print("i:", i)
print("omega:", omega)
print("w:", w)

# Got this to run correctly, should return proper values for all variables in order a, e, i, omega, w.

# Semi major axis                   - Expected = 1.056800055744403  Calculated = 1.05680005574156 
# Eccentricity                      - Expected = .3442331093278357  Calculated = .344233109330765
# Inclination Angle                 - Expected = 25.15525166563580  Calculated = 25.15525166563579
# Longitude of the Ascending Node   - Expected = 236.2379806531768  Calculated = 236.2379806531768
# Argument of the Perihelion        - Expected = 255.5046142809498  Calculated = 255.50461428065267

# None of these returned percent errors of any significance (< .001%)

print("Percent Error of Semi Major Axis:", percent_error(1.056800055744403, 1.05680005574156), "%")
print("Percent Error of Eccentricity:", percent_error(.3442331093278357, .344233109330765), "%")
print("Percent Error of Inclination Angle:", percent_error(25.15525166563580, 25.15525166563579), "%")
print("Percent Error of Longitude of the Acending Node:", percent_error(236.2379806531768, 236.2379806531768), "%")
print("Percent Error of Argument of the Perihelion:", percent_error(255.5046142809498, 255.50461428065267), "%")
