

# Translating the provided Julia code into Python

import numpy as np

def myrandom(r, a, b):
    r = (r * 41475557) % 1
    return r, a + r * (b - a)

n = 99
dens = 0.5
dvert = 2
myseed = 1

for igen in range(5):
    r = (4 * myseed + 1) / 16384 / 16384
    myseed += 1

    Fc = np.zeros((n+1, n+1))
    Fl = np.zeros((n+1, n+1))
    F = np.zeros((n+1, n+1))

    for i in range(n):
        for j in range(i+1, n+1):
            r, num = myrandom(r, 0, 1)
            if num < dens:
                r, num = myrandom(r, 0, 10)
                Fc[i, j] = num
                Fc[j, i] = num
            else:
                r, num = myrandom(r, -10, 0)
                Fc[i, j] = num
                Fc[j, i] = num

    for i in range(n+1):
        r, num = myrandom(r, 0, dvert)
        Fl[i, i] = num

    for i in range(n):
        for j in range(i+1, n+1):
            Fl[i, j] = 0.5 * (Fl[i, i] + Fl[j, j])
            Fl[j, i] = Fl[i, j]

    for i in range(n+1):
        for j in range(i, n+1):
            F[i, j] = Fl[i, j] - Fc[i, j]
            F[j, i] = F[i, j]

    filename = f"Problem_{n+1}x{n+1}({dens})_{igen+1}.txt"
    np.savetxt(filename, F, fmt='%.15f')
    
# This code will generate 5 files with matrices. The matrices are not displayed here to save space.
