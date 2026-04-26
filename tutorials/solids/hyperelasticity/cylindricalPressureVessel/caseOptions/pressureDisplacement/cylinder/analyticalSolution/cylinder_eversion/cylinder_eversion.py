#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt


################################################################################

"""
# test
A       = 1.0           # inner radius, initial configuration [m]
B       = 2.0           # outer radius, initial configuration [m]
a       = A             # inner radius, deformed configuration [m]
b       = B             # outer radius, deformed configuration [m]
lz      = 1.0           # axial stretch [-]
mu      = float(1e6)    # shear modulus [Pa]

#p0     = float(5e4)    # internal pressure [Pa]
"""

# test
A       = 8.0           # inner radius, initial configuration [m]
B       = 16.0          # outer radius, initial configuration [m]
a       = A             # inner radius, deformed configuration [m]
b       = B             # outer radius, deformed configuration [m]
lz      = 1.0           # axial stretch [-]
mu      = float(1e9) / 3.0   # shear modulus [Pa]

#p0     = float(100e6)    # internal pressure [Pa]

"""
# Bijelonja et al. (they use the Mooney-Rivlin model)
A       = 7.0           # inner radius, initial configuration [m]
B       = 18.625        # outer radius, initial configuration [m]
a       = A             # inner radius, deformed configuration [m]
b       = B             # outer radius, deformed configuration [m]
lz      = 1.0           # axial stretch [-]
mu      = float(160e6)   # shear modulus [Pa]

p0     = float(100e6)    # internal pressure [Pa]
"""

"""
# aorta
A       = 0.02          # inner radius, initial configuration [m]
B       = 0.021         # outer radius, initial configuration [m]
a       = A             # inner radius, deformed configuration [m]
b       = B             # outer radius, deformed configuration [m]
lz      = 1.0           # axial stretch [-]
mu      = float(3e6)    # shear modulus [Pa]

#p0     = float(1e4)    # internal pressure [Pa]
"""

################################################################################


# 'f(r)'
def f1(p0, a):

    la      = lz * a
    la2     = la * a
    l2a2    = la * la

    Z       = la2 + B*B - A*A
    ZB      = Z - B*B

    return p0                                                                  \
         - mu * np.log(B/A) / lz                                               \
         + 0.5 * mu * np.log(Z / la2) / lz                                     \
         + 0.5 * mu * ZB / (Z * lz)                                            \
         - 0.5 * mu * ZB / l2a2


# 'df(r)/dr'
def df1dr(p0, a):

    la      = lz * a
    la2     = la * a
    l2a3    = la * la2

    Z       = la2 + B*B - A*A
    ZB      = Z - B*B

    return mu * a / Z                                                          \
         + mu * a * B*B / (Z*Z)                                                \
         + mu * ZB / l2a3                                                      \
         - 2.0 * mu / la


# Newton-Raphson method (single variable only)
def newton(p0, f, dfdx, x0 = 0.0, tol = 1e-9,):

    x_old = x0
    x_new = x0 + 2.0 * tol

    while abs(x_new - x_old) > tol:

        x_old = x_new
        x_new = x_old - f(p0, x_old) / dfdx(p0, x_old)

        #print("x    = {:.9f}".format(x_new))

    return x_new


# main function to call from outside
def evert(p0):

    return newton(p0, f1, df1dr, x0 = A)


# deformed radius
def radius(R, r0):

    return np.sqrt(r0*r0 + (R*R - A*A) / lz)


# deformed angle
def angle(t):

    return t


# deformed length
def length(z):

    return lz * z


# radial stretch
def lambda_r(R, r0):

    return R / (lz * radius(R, r0))


# tangential stretch
def lambda_t(R, r0):

    return radius(R, r0) / R


def pressure(R, r0, p0):

    r = radius(R, r0)

    return p0                                                                  \
         - mu * np.log(R * r0 / (r * A)) / lz                                  \
         + mu * np.power(R / (lz * r), 2)                                      \
         + mu * (lz * r0*r0 - A*A) * (1.0 / (r*r) - 1.0 / (r0*r0)) / (2.0 * lz*lz)


################################################################################

#print("r0 = {:.10f}".format(evert(p0)))
"""
p0 = np.linspace(0.0, 100e6, 11)
r0 = np.zeros(11)

for i, rad in enumerate(r0):

    r0[i] = newton(p0[i], f1, df1dr, x0 = A)

plt.grid(True, 'both')
plt.plot(p0, r0 - A)
plt.show()
"""
"""
R = np.linspace(A, B, 10)
r0 = evert(p0)

print(radius(R, r0))
print(pressure(R, r0, p0))

#plt.grid(True, 'both')
#plt.plot(R, radius(R, r0))
#plt.plot(R, sigma_rr(R, r0, p0))
#plt.show()
"""
