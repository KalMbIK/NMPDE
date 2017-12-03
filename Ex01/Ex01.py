import numpy as np
from math import factorial, pow, fabs


def myDot(a, b):
    if len(a) != len(b):
        raise Exception
    sum = 0.
    for i in range(len(a)):
        sum += a[i] * b[i]
    return sum


def taylor(k_, m_):
    return pow(k_, m_) / factorial(m_)


def generateTaylorRow(gr, i):
    temp = []
    for p in gr:
        temp.append(taylor(p, i))
    return temp


def generateTaylorTable(gr):
    temp = []
    for i in range(len(gr)):
        temp.append(generateTaylorRow(gr, i))
    return np.array(temp)


def generateCoeffs(gr, m_):
    tTaylor = generateTaylorTable(gr)
    b = np.zeros(len(tTaylor))
    b[m_] = 1
    return np.linalg.solve(tTaylor, b)


def getAccuracy(coeffs, gr, m_, tol):
    tR = len(gr)
    tTR = generateTaylorRow(gr, tR)
    dt = fabs(myDot(tTR, coeffs))
    while dt <= tol:
        # print dt
        tR += 1
        dt = fabs(myDot(generateTaylorRow(gr, tR), coeffs))
    return tR - m_


# program parameters
p = 1
q = 3
m = 2
grid = np.arange(-p, q + 1)
size = np.size(grid)

if m > size:
    raise AttributeError("input parameters p, q should "
                         "correspond to sufficiently many points to support an approximation")

coeffs = generateCoeffs(grid, m)
accuracy = getAccuracy(coeffs, grid, m, tol=10 ** -8)

print "Coeffs: " + str(coeffs)
print "Accuracy: " + str(accuracy)
