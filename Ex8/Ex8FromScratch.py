import numpy as np
import math as m

class FEM2D(object):
    def __init__(self, filename, exSol, rhs):
        self.parseFile(filename)
        self.exSolFun = exSol
        self.rhsFun = rhs
        self.exSolOnGrid = np.array(map(exSol,self.grid[:self.nodesInt]))
        self.rhsOnGrid = np.array(map(rhs,self.grid[:self.nodesInt]))
        self.initEngine()
        self.genMatrix()
        self.genRhsV()
        self.numSolOnGrid = np.linalg.solve(self.A, self.rhsV)
        self.findL2Norm(self.numSolOnGrid-self.exSolOnGrid)

    def getMatrix(self):
        return self.A
    def getRhsV(self):
        return self.rhsV
    def getExactSol(self):
        return self.exSolOnGrid
    def getNumSol(self):
        return self.numSolOnGrid
    def getErrors(self):
        return self.errors

    def initEngine(self):
        self.B = np.zeros((self.numTrgs,2,2))
        self.trigDiams = np.zeros(self.numTrgs)
        self.trigIters = [0,1,2]
        p0 = [-1,-1]
        p1 = [1, 0]
        p2 = [0, 1]
        self.p = [p0,p1,p2]

    def setBk(self, k):
        x0 = self.grid[self.npt[k][0]]
        x1 = self.grid[self.npt[k][1]]
        x2 = self.grid[self.npt[k][2]]
        self.B[k][0][0] = x2[1] - x0[1]
        self.B[k][0][1] = x0[1] - x1[1]
        self.B[k][1][0] = x0[0] - x2[0]
        self.B[k][1][1] = x1[0] - x0[0]

    def genMatrix(self):
        self.A = np.zeros((self.nodesInt,self.nodesInt))
        for k in range(self.numTrgs):
            self.setBk(k)
            self.trigDiams[k] = np.linalg.det(self.B[k]) / 2
            b = self.B[k]
            for l in self.trigIters:
                i = self.npt[k][l]
                if i >= self.nodesInt:
                    continue
                for m in range(0,3):
                # for m in range(l,3):
                    j = self.npt[k][m]
                    if j >= self.nodesInt:
                        continue
                    alpha = ((b[0][0]*self.p[l][0]+b[0][1]*self.p[l][1])*(b[0][0]*self.p[m][0]+b[0][1]*self.p[m][1])+(b[1][0]*self.p[l][0]+b[1][1]*self.p[l][1])*(b[1][0]*self.p[m][0]+b[1][1]*self.p[m][1]))/4/self.trigDiams[k]
                    self.A[i][j] += alpha
                    # self.A[j][i] += alpha

    def integrateFun(self, fValues, storage):
        for i in range(self.numTrgs):
            for k in self.trigIters:
                ind = self.npt[i][k]
                if ind >= self.nodesInt:
                    continue
                storage[ind] += self.trigDiams[i] / 3 * fValues[ind]

    def genRhsV(self):
        self.rhsV = np.zeros(self.nodesInt)
        self.integrateFun(self.rhsOnGrid,self.rhsV)

    def findL2Norm(self, funValues):
        self.errors = np.zeros_like(self.rhsV)
        self.integrateFun(funValues**2,self.errors)
        intVal = np.sum(self.errors)
        self.norm = m.sqrt(intVal)

    def meshToStr(self):
        s = 'Number of nodes ' + str(self.nodes)
        s += '\r\nNumber of internal nodes ' + str(self.nodesInt)
        s += '\r\nNumber of triangles in the mesh ' + str(self.numTrgs)
        s += '\r\nGRID\r\n'
        s += str(self.grid)
        s += '\r\n______________'
        s += '\r\nNPT\r\n'
        s += str(self.npt)
        s += '\r\n______________'
        return s

    def solStr(self, sol):
        s = 'size = ' + str(len(sol)) + '\r\n'
        s += str(sol)
        s += '\r\n______________'
        return s

    def solToStr(self):
        s = 'Norm ' + str(self.norm)
        s += '\r\nEXACT\r\n'
        s += self.solStr(self.exSolOnGrid)
        s += '\r\nNUMERIC\r\n'
        s += self.solStr(self.numSolOnGrid)
        return s

    def parseFile(self,fileName):
        f = file(fileName)
        line = f.readline().split()
        self.nodes = int(line[0])
        self.nodesInt = int(line[1])
        self.grid = np.zeros((self.nodes,2), dtype=np.float64)
        for i in range(self.nodes):
            line = f.readline().split()
            self.grid[i][0]=float(line[0])
            self.grid[i][1]=float(line[1])
        self.numTrgs = int(f.readline())
        self.npt = np.zeros((self.numTrgs,3), dtype=np.uint)
        for i in range(self.numTrgs):
            line = f.readline().split()
            self.npt[i][0]=int(line[0])-1
            self.npt[i][1]=int(line[1])-1
            self.npt[i][2]=int(line[2])-1
        f.close()

fileName1 = 'MESH.str.10'
fileName2 = 'MESH.unstr.8'
exSol = lambda point: m.sin(m.pi*point[0])*m.sin(m.pi*point[1])
rhs = lambda point: exSol(point)*2.*m.pi*m.pi
fem1 = FEM2D(fileName2, exSol, rhs)
print fem1.solToStr()