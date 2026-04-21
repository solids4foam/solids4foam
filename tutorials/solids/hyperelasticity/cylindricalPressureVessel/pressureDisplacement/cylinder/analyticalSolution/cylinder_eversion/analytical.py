#!/usr/bin/env python3

import numpy as np
import sys
import os

import cylinder_eversion as cyl

def make_probes(r0, r1, num=8):
    dr = (r1 - r0) / float(num)

    r = r0
    i = 0
    probes = np.empty([num + 1, 3])

    while i <= num:
        probes[i] = np.array([r, 0.0, 0.0])

        i = i + 1
        r = r0 + float(i) * dr

    return probes


def print_vector_tabbed(vec, outfile, skip = False, fmt_string = "{:.9e}"):
    for i in vec:
        print(fmt_string.format(i), end = '\t', file = outfile)
    if skip:
        print(file = outfile)


def displacement(pos, r0):

    #return np.tensordot(F, pos, axes=1) - pos
    return np.array([cyl.radius(pos[0], r0),
                     cyl.angle(pos[1]),
                     cyl.length(pos[2])
                    ]) - pos


def cyl2cart(Q, t):
    if t.ndim == 1:
        return np.matmul(np.matmul(np.transpose(Q), t), Q)
    else:
        return np.tensordot(Q, t, axes=1)


def print_header(outfile):
    with open(outfile, "w") as of:
        #print("# x\ty\tz\tp\tsigmaEq\tUx\tUy\tUz\tsigmaxx\tsigmayy\tsigmazz\tsigmaxy\tsigmaxz\tsigmayz", file = of)
        print("# x\ty\tz\tp\tsigmaEq\tUr\tUt\tUz\tsigmarr\tsigmatt\tsigmazz\tsigmart\tsigmarz\tsigmatz", file = of)


def print_output(pos, r0, p0, outfile):

    with open(outfile, "a") as of:

        R = np.sqrt(pos[0]*pos[0] + pos[1]*pos[1])
        T = np.arctan2(pos[1], pos[0])
        Z = pos[2]
        pos_cyl = np.array([R, T, Z])

        # cylindrical to cartesian
        Q = np.array([np.array([ np.cos(T), -np.sin(T),  0.0 ]),
                      np.array([ np.sin(T),  np.cos(T),  0.0 ]),
                      np.array([ 0.0,        0.0,        1.0 ])
                     ])

        stretch = np.array([cyl.lambda_r(R, r0),
                            cyl.lambda_t(R, r0),
                            cyl.lz])
        I = np.identity(3)
        F = np.diag(stretch)
        b = np.matmul(F, np.transpose(F))
        p = cyl.pressure(R, r0, p0)

        sigma_cyl = cyl.mu*b - p*I
        U_cyl = displacement(pos_cyl, r0)

        #sigma = cyl2cart(Q, sigma_cyl)
        #U = cyl2cart(Q, U_cyl)

        #sigmaEq = np.sqrt(1.5 * np.sum(np.square((sigma - 1.0/3.0 * np.trace(sigma) * I))))
        sigmaEq = np.sqrt(1.5 * np.sum(np.square((sigma_cyl - 1.0/3.0 * np.trace(sigma_cyl) * I))))

        # print xyz
        print_vector_tabbed(pos, of)
        # print p
        print("{:.9e}\t".format(p), end = "", file = of)
        # print sigmaEq
        print("{:.9e}\t".format(sigmaEq), end = "", file = of)
        # print U(x,y,z)
        #print_vector_tabbed(U, of)
        print_vector_tabbed(U_cyl, of)
        # print sigma(xx,yy,zz)
        #print_vector_tabbed(sigma.diagonal(), of)
        print_vector_tabbed(sigma_cyl.diagonal(), of)
        # print sigma(xy,xz,yz)
        #print_vector_tabbed(np.array([sigma[0, 1], sigma[0, 2], sigma[1,2]]), of, True)
        print_vector_tabbed(np.array([sigma_cyl[0, 1], sigma_cyl[0, 2], sigma_cyl[1,2]]), of, True)


# check if we're running from case root dir
# check cmdline args
if len(sys.argv) < 2:
    print("pressure value not supplied")
    exit(1)
elif len(sys.argv) > 3:
    print ("Too many arguments")
    exit(1)
elif not isinstance(float(sys.argv[1]), float):
    print ("Supply floating point value for pressure")
    exit(1)


# check output directory
if not os.path.isdir("results"):
    os.mkdir("results")

p0 = float(sys.argv[1])

r0 = cyl.evert(p0)

#print("a: ", r0)

probe_positions = make_probes(8.0, 16.0, num=32)

# compute and print output
outfile = "results/analytical"
print_header(outfile)

for pos in probe_positions:
    print_output(pos, r0, p0, outfile)
