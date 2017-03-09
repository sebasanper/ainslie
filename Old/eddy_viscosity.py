__author__ = 'sebasanper'

from math import sqrt, exp
import time
# from 2Dintegrate_simpson import  simpson_integrate

def ainslie(ct, u0, distance_parallel, distance_perpendicular):

    # centreline = open('centreline.dat', 'w')
    # velocity = open('velocity.dat', 'w')

    h = 0.1
    L = distance_parallel
    n = int(L / h) + 1
    Uc1 = [0.0 for x in range(n)]
    d1 = [0.0 for x in range(n)]
    I0 = 8.0  # Must be given in percent.
    k = 0.41  # von Karman constant
    Ct = ct  # Thrust coefficient
    U0 = u0
    # dr = 0.1
    Y = distance_perpendicular
    # m = int(Y / dr)

    Dmi = Ct - 0.05 - (16.0 * Ct - 0.5) * I0 / 1000.0

    def b(deficit):  # Wake width measure
        return (3.56 * Ct / (8.0 * deficit * (1.0 - 0.5 * deficit))) ** 0.5

    def F(x):  # Factor for near and far wake
        if x >= 5.5:
            return 1.0
        if x < 5.5:
            if x >= 4.5:
                return 0.65 + ((x - 4.5) / 23.32) ** (1.0 / 3.0)
            else:
                return 0.65 - ((-x + 4.5) / 23.32) ** (1.0 / 3.0)

    def E(x1, Uf, Ud, Dm):  # Eddy viscosity term
        return F(x1) * ((0.015 * b(Dm) * (Uf - Ud)) + (k ** 2.0) * I0 / 100.0)

    # Uc = U0 * (1.0 - Dmi)  # Boundary condition at x = 2.0
    # d = Dmi
    Uc1[0] = U0 * (1.0 - Dmi)  # Boundary condition at x = 2.0
    d1[0] = Dmi
    for i in range(1, n):  # For all positions in the wake centreline direction. Recursive. Whole grid
        Uc1[i] = Uc1[i - 1] + (h * 16.0 * E(i * h, U0, Uc1[i - 1], d1[i - 1]) * (Uc1[i - 1] ** 3.0 - U0 * Uc1[i - 1] ** 2.0 - Uc1[i - 1] * U0 ** 2.0 + U0 ** 3.0) / (Uc1[i - 1] * Ct * U0 ** 2.0))
        d1[i] = 1.0 - Uc1[i] / U0
        # centreline.write('{0:f}\t{1:f}\t{2:f}\n'.format(i * h + 2.0, Uc1[i], d1[i]))
        # for j in range(m):  # For all positions in the perpendicular direction.
        #     U = U0 * (1.0 - d1[i] * exp(-3.56 * (j * dr / b(d1[i])) ** 2.0))
        #     velocity.write('{0:f}\t{1:f}\t{2:f}\n'.format(i * h + 2.0, j * dr, U))
        # velocity.write('\n')

    ########### Code to calculate wake deficit at a specific point instead of the whole grid.

    U = U0 * (1.0 - d1[n - 1] * exp(- 3.56 * (Y / b(d1[n - 1])) ** 2.0))
    return U

    # centreline.close()
    # velocity.close()

if __name__ == '__main__':
    with open('centreline.dat', 'w') as out:
        for i in range(2000):
            out.write('{0}\n'.format(ainslie(0.79, 8.5, i * 0.01, 0.8)))

