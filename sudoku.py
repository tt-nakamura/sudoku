import numpy as np
import sympy as sp

def sudoku(A):
    """ solve sudoku puzzle
    A = initial state matrix
    nrows(A) == ncols(A) == square number
    empty places must be filled with zeros
    return B = final state matrix
    referece: W. Decker and G. Pfister
      "A First Course in Computational Algebraic Geometry"
       chapter 3
    """
    d = len(A); m = int(d**0.5); n = d*d

    x = sp.symbols('x:%d'%n)
    F = [sp.prod([t-i for i in range(1,d+1)]) for t in x]

    for a in range(n):
        ia,ja = divmod(a,d)
        ka = (ia//m)*m + ja//m
        for b in range(a+1,n):
            ib,jb = divmod(b,d)
            kb = (ib//m)*m + jb//m
            if ia==ib or ja==jb or ka==kb:
                t = sp.cancel((F[a]-F[b])/(x[a]-x[b]))
                F.append(t)

    c = np.reshape(A,-1)
    for a in np.flatnonzero(c):
        F.append(x[a]-c[a])

    G = sp.groebner(F, *x, order='grevlex')
    B = [int(-sp.Poly(g).TC()) for g in G]
    return np.reshape(B,(d,d))

if __name__ == '__main__':
    A = [[0,0,0,4],
         [4,0,2,0],
         [0,3,0,1],
         [1,0,0,0]]
    print(sudoku(A))
