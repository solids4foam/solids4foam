import numpy as np
import matplotlib.pyplot as plt

#inner radius
r_i=8
#outer radius
r_o=16
#number of samples
num_samples=50
#youngs
E = 1000 * 10**6
#pressure
p_i = 10 * 10**6
#poisson
nu = 0.5
#shear
G = E / (2 * (1 + nu))
#delta
delta = r_o - r_i
#current r
r = np.zeros(num_samples, dtype=float)
#radial stress
sigma_rr = np.zeros(num_samples, dtype=float)
# hoop stress
sigma_ff = np.zeros(num_samples, dtype=float)
# axial stress
sigma_z = np.zeros(num_samples, dtype=float)

# radial displacements
u_r = np.zeros(num_samples, dtype=float)
#pressure
p = np.zeros(num_samples, dtype=float)
#f(x)
y = np.zeros(num_samples, dtype=float)

for i in range(0, num_samples):
    if i==0:
        r[i] = r_i
    else:
        r[i] = r[i-1] + delta/(num_samples-1)

    r[i]=np.round(r[i], decimals=3)
    #radial stress
    sigma_rr[i] = p_i * r_i**2 /(r_o**2 - r_i**2) * (1 - r_o**2/r[i]**2)

    #hoop stress
    sigma_ff[i] = p_i * r_i**2 /(r_o**2 - r_i**2) * (1 + r_o**2/r[i]**2)

    #radial displacements
    u_r[i] = (1 + nu) / E * p_i * r_i**2 /(r_o**2 - r_i**2) * (r_o**2 / r[i] + (1 - 2*nu) * r[i])

    #axial stress
    sigma_z[i] = p_i * r_i**2 /(r_o**2 - r_i**2)
    #pressure
    p[i] = 1/3 * (sigma_rr[i] + sigma_ff[i] + sigma_z[i])



a = np.transpose(np.array([[r], [sigma_rr], [sigma_ff], [sigma_z], [u_r], [p]]))
mat = np.matrix(a)

text = np.array(["#position", "radial stress", "hoop stress", "axial stress", "displacements", "pressure"], dtype='str')
mat = np.vstack([text, mat])


# Save the results to a text file
np.savetxt('result' + str(p_i/10**6) + 'MPa.txt', mat, fmt='%s', delimiter="\t")


