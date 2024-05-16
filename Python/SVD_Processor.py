# -*- coding: utf-8 -*-
"""Solve linear systems based on SVD, pseudo-inverse and optimization.

This code is based on "Linear Algebraic Equations, SVD, and the Pseudo-Inverse"
by Philip N. Sabes (october, 2001) available on Research Gate. It solves
under-constrained problems subject to constraints, for instance muscle tension
problems. The 'main' function of this module runs a simplistic musculosqueletal
model for test.
"""
# =====================================
# IMPORTS
# =====================================
import logging
import numpy as np
from scipy.optimize import minimize

# =====================================
# CONSTANTS AND GLOBAL VALUES
# =====================================
lgg = logging.getLogger(__name__)
A = np.array([[1, -1, 0,  0, 0,  0],
              [0,  0, 1, -1, 0,  0],
              [0,  0, 0,  0, 1, -1]], dtype=float)
Y = 2*np.ones(3, dtype=float)


# =====================================
# GLOBAL FUNCTION
# =====================================
def cost_function(wgt, tTilde, svdProc):
    """Test cost function for muscle tension optimization."""
    return np.linalg.norm(tTilde + svdProc.computeKernelLinComb(wgt))


# =====================================
# CLASS SVD_Processor
# =====================================
class SVD_Solver():
    """Basic features to solve under-constrained linear system.

    The under-contrained linear system is formulated as Mx=y, with:
    - M in IR^(n x m) a real matrix, n < m;
    - x in IR^m
    - y in IR^n
    """

    # --------------------------------------------------------------------------
    def __init__(self, M):
        """Read the input matrix and run SVD.

        Input
        -----
        - M: np.array( (n x m) ), n < m, the linear system
        """
        lgg.debug("Initializing solver...")
        self.n, self.m = M.shape
        self.r = np.linalg.matrix_rank(M)
        self.M = M

        self.U, self.D, self.VT = np.linalg.svd(M)
        lgg.debug("SVD decomposition:\nU =\n%s\nD =\n%s\nVT =\n%s",
                  self.U, self.D, self.VT)

        self.Mdag = np.linalg.pinv(M)
        lgg.debug("Pseudo-inverse:\nMdag =\n%s", self.Mdag)

    # --------------------------------------------------------------------------
    def computeBasicSolution(self, y):
        """Compute a solution (not unique) in the kernel of the linear system.

        Input
        -----
        - y: np.array(n), the vector to solve for

        Output
        ------
        - y0: np.array(n), one vector solution in the kernel of M

        """
        assert self.M.shape[0] == self.n, "Invalid data size !"
        lgg.debug("Computing kernel basic (least square) solution")
        xTilde = self.Mdag @ y
        lgg.debug("xTilde =\n%s", xTilde)
        lgg.debug("M @ xTilde =\n%s", self.M @ xTilde)
        lgg.debug("y =\n%s", y)
        return xTilde

    # --------------------------------------------------------------------------
    def computeKernelLinComb(self, wgt):
        """Compute a weighted linear combination of kernel base vectors.

        Input
        -----
        - wgt: np.array(m-n) of float, the vector of weights

        Output
        ------
        - x0: np.array(n), the weighted linear combination of kernel's base
                           vectors
        """
        # lgg.debug("wgt (len %d):\n%s", len(wgt), wgt)
        # lgg.debug("M ! n=%d x m=%d", self.M.shape[0], self.M.shape[1])
        p = self.m - self.n
        assert p == len(wgt), "Invalid data size !"
        v = np.zeros(self.m)

        # build column by column
        for k in range(p):
            v += wgt[k]*self.VT[self.r+k, :]
        # lgg.debug("v_kernel =\n%s", v)
        # lgg.debug("M.v_kernel =\n%s", np.dot(self.M, v))
        return v


# =====================================
# MAIN CODE
# =====================================


def main():
    lgg.debug("Test matrix:\n%s", A)
    r = np.linalg.matrix_rank(A)
    lgg.debug("Rank: %d", r)
    n, m = A.shape

    mySVD = SVD_Solver(A)
    xTilde = mySVD.computeBasicSolution(Y)
    lgg.debug("VT:\n%s", mySVD.VT)

    lgg.debug("Testing optimization...")
    res = minimize(cost_function,
                   x0=np.random.randn(m-n),
                   args=(xTilde, mySVD))
    lgg.info("Optimization result:\n%s", res)
    xOpt = xTilde + mySVD.computeKernelLinComb(res.x)
    lgg.info("Optimal solution: %s", xOpt)
    lgg.info("A.xopt =\n%s", A @ xOpt)
    lgg.info("Y =\n%s", Y)

    lgg.info("Done...")


# =====================================
if __name__ == "__main__":
    # logging stuff
    lgg.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        '%(name)s - %(levelname)s - %(filename)s:%(lineno)s\t%(message)s')
    ch.setFormatter(formatter)
    lgg.addHandler(ch)

    main()
