import numpy as np
import math as m

def parseFile(fileName):
    f = file(fileName)
    line = f.readline().split()
    nodes = int(line[0])
    nodesInt = int(line[1])
    grid = np.zeros((nodes,2), dtype=np.float64)
    for i in range(nodes):
        line = f.readline().split()
        grid[i][0]=float(line[0])
        grid[i][1]=float(line[1])
    numTrgs = int(f.readline())
    npt = np.zeros((numTrgs,3), dtype=np.uint)
    for i in range(numTrgs):
        line = f.readline().split()
        npt[i][0]=int(line[0])-1
        npt[i][1]=int(line[1])-1
        npt[i][2]=int(line[2])-1
    f.close()
    return nodes, nodesInt, numTrgs, grid, npt

def setBk(B, grid, npt, k):
    x0 = grid[npt[k][0]]
    x1 = grid[npt[k][1]]
    x2 = grid[npt[k][2]]
    B[k][0][0] = x2[1]-x0[1]
    B[k][0][1] = x0[1]-x1[1]
    B[k][1][0] = x0[0]-x2[0]
    B[k][1][1] = x1[0]-x0[0]

def genMatrixAndRhs(N, Nint, Nt, grid, npt):
    A = np.zeros((Nint,Nint))
    B = np.zeros((Nt,2,2))
    trigDiams = np.zeros(Nt)
    trigIters = [0,1,2]
    p0 = [-1,-1]
    p1 = [1, 0]
    p2 = [0, 1]
    p = [p0,p1,p2]

    # matrix
    for k in range(Nt):
        setBk(B,grid,npt,k)
        trigDiams[k] = np.linalg.det(B[k]) / 2
        b = B[k]
        for l in trigIters:
            i = npt[k][l]
            if i >= Nint:
                continue
            for m in range(0,3):
                j = npt[k][m]
                if j >= Nint:
                    continue
                alpha = ((b[0][0]*p[l][0]+b[0][1]*p[l][1])*(b[0][0]*p[m][0]+b[0][1]*p[m][1])+(b[1][0]*p[l][0]+b[1][1]*p[l][1])*(b[1][0]*p[m][0]+b[1][1]*p[m][1]))/4/trigDiams[k]
                A[i][j] += alpha

    # f
    f = np.zeros(Nint)
    trigIters = [0,1,2]
    for i in range(Nt):
        for k in trigIters:
            ind = npt[i][k]
            if ind >= Nint:
                continue
            f[ind] += trigDiams[i] / 3 * rhs(grid[ind])

    return A, f, trigDiams

def integrate(fun, npt, Nint, td, Nt):
    f = np.zeros(Nint)
    trigIters = [0,1,2]
    for i in range(Nt):
        for k in trigIters:
            ind = npt[i][k]
            if ind >= Nint:
                continue
            f[ind] += td[i] / 3 * fun[ind]
    return f

def findL2Norm(fun, npt, Nint, td, Nt):
    ff = fun**2
    u = integrate(ff, npt, Nint, td, Nt)
    intVal = np.sum(u)
    return m.sqrt(intVal)

def solveSystem(N, Nint, Nt, grid, npt):
    A, f, td = genMatrixAndRhs(N, Nint, Nt, grid, npt)
    return np.linalg.solve(A,f), td

fileName1 = 'MESH.str.10'
fileName2 = 'MESH.unstr.8'
N, Nint, Nt, grid, npt = parseFile(fileName1)
exSol = lambda point: m.sin(m.pi*point[0])*m.sin(m.pi*point[1])
rhs = lambda point: exSol(point)*2.*m.pi*m.pi

exSolOnGrid = np.array(map(exSol,grid[:Nint]))
numSol, td = solveSystem(N, Nint, Nt, grid, npt)
print 'L2norm = ' + str(findL2Norm(numSol-exSolOnGrid, npt, Nint, td, Nt))
print(3.4E-2)