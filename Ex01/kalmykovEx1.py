import numpy as np
from math import factorial, pow

class DiffApproxGenerator(object):
    def __init__(self, m, p, q):
        self.leftBorder = p
        self.rightBorder = q
        self.order = m
        self.grid = np.arange(-p, q+1)
        self.size = np.size(self.grid)
        if m >= self.size:
            raise AttributeError("input parameters p, q should "
                                 "correspond to sufficiently many points to support an approximation")
        self.vTaylor = np.vectorize(self.taylor)
        self.tTaylor = self.generateTaylorTable()
        self.coeffs = self.generateCoeffs()

    def getGrid(self):
        return self.grid

    def taylor(self, k_, m_):
        return pow(k_,m_)/factorial(m_)

    def generateTaylorTable(self):
        temp, gr = [], self.getGrid()
        for i in range(self.size):
            temp.append(self.vTaylor(gr,i))
        return np.array(temp)

    def generateCoeffs(self):
        b = np.zeros(self.size)
        b[self.order] = 1
        return np.linalg.solve(self.tTaylor,b)

    def getAccuracy(self):
        if self.order%2 == 0:
            return self.size+1-self.order
        else:
            return self.size-self.order

symmCase = DiffApproxGenerator(1,0,3)
print symmCase.size
print symmCase.getGrid()
print symmCase.tTaylor
print "Coeffs: " + str(symmCase.coeffs)
print "Accuracy: " + str(symmCase.getAccuracy())