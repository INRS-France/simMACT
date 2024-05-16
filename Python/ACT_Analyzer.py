# -*- coding: utf-8 -*-
"""Abstract class for the simulation of maximum joint actuation torques."""

# =============================
# DEPENDENCIES
# =============================
from os.path import join, dirname, abspath
from abc import ABC, abstractmethod  # Abstract Base Class

import sys
import math
import logging
import logging.config
from scipy.optimize import minimize
from matplotlib import cm

import numpy as np
import matplotlib.pyplot as plt

# append current folder in order to import local modules
sys.path.append(dirname(abspath(__file__)))
import ZNT_Processor
import MSM_Processor
import SVD_Processor
import GUI_Manager


# =============================
# CONSTANTS
# =============================
PI = math.pi
RAD = np.radians
DEG = np.degrees
DEF_MARM_THR = 1e-4     # Moment arm actuation threshold. Muscles with lower
# value are considered as not actuating the DoF

lgg = logging.getLogger(__name__)   # if module not run as 'main'
DEF_LOG_CONF = './logConfig/ACT_log.conf'   # else
DEF_LOGGER = "ACT_logger"

# Opensim models' specific constants
# ----------------------------------
DEF_MDL_DIR = join(dirname(dirname(abspath(__file__))), "Models")
MOBL_MDL_FILE = join(DEF_MDL_DIR, "MOBL/4.1/MOBL_ARMS_fixed_41-ignoreTendonCompl.osim")

# JointTVA specific constants (for test of the module)
# ----------------------
# Default joint TVA range
ANL_DOFS = ["elbow_flexion", "pro_sup"]
DEF_VAR_DOF = "elbow_flexion"
DEF_Q_RANGE = np.arange(0., 151., 30.)   # joint angles (deg)
DEF_V_RANGE = np.arange(-150., 151., 50.)   # joint velocities (deg/s)


# ==============================================================================
# General purpose functions
# ==============================================================================

# --------------------------------------------------------------------------
def muscleTensions_cost_fun(wgt, tTilde, svdProc):
    """Compute the cost function of muscle tensions.

    Cost function for minimization of the solution of an under-constrained
    linear problem defined by a (n x q) matrix A, q < n, r=rank(A),
    and its goal y_des : A.t=y_des.

    Args:
        wgt (np.array): the vector of dim 'r', i.e the weights that build a
            linear combination of the base vectors of the kernel of a linear
            system of rank 'r' described in svdProc.
        tTilde (np.array): the least square solution of A.t=y_des of dimension
            'm'
        svdProc (SVD_Processor) : the SVD_Processor object that describes the
            linear system to solve.

    Returns:
        (float): the value of the cost function
    """
    return np.linalg.norm(compute_muscle_actuation(wgt, tTilde, svdProc))


# --------------------------------------------------------------------------
def compute_muscle_actuation(wgt, tTilde, svdProc):
    r"""Compute muscle tensions.

    Compute actuating muscles tensions from SVD-based arguments.

    Args:
        wgt (np.array): the vector of dim 'r', i.e the weights that build a
            linear combination of the base vectors of the kernel of a linear system
            of rank 'r' described in svdProc.
        tTilde (np.array): the least square solution of A.t=y_des of dimension 'm'
        svdProc (SVD_Processor): the SVD_Processor object that describes the
            linear system to solve

    Returns:
        (float): the individual muscle tensions
    """
    # tBar = tTilde + sum_i wgt_i_ V_i with V_i the basis' kernel vectors of SVD
    return tTilde + svdProc.computeKernelLinComb(wgt)


# --------------------------------------------------------------------------
def lb_constr_fun(wgt, tTilde, svdProc, lb):
    """Return the lower bound of the inequality constraint.

    Args:
        wgt (np.array): the vector of dim 'r', i.e the weights that build a
            linear combination of the base vectors of the kernel of a linear
            system of rank 'r' described in svdProc.
        tTilde (np.array): the least square solution of A.t=y_des of dimension
            'm'
        svdProc (SVD_Processor): the SVD_Processor object that describes the
            linear system to solve
        lb (np.array): the vector of dim 'm' of the lower-bound tension of each
            muscle.

    Returns:
        (np.array): lower bound of the muscle actuations that obey the linear
            system.
    """
    return compute_muscle_actuation(wgt, tTilde, svdProc) - lb


# --------------------------------------------------------------------------
def ub_constr_fun(wgt, tTilde, svdProc, ub):
    """Return the upper bound of the inequality constraint.

    Args:
        wgt (np.array): the vector of dim 'r', i.e the weights that build a
            linear combination of the base vectors of the kernel of a linear
            system of rank 'r' described in svdProc.
        tTilde (np.array): the least square solution of A.t=y_des of dimension
            'm'
    svdProc (SVD_Processor): the SVD_Processor object that describes the linear
        system to solve
    ub (np.array): the vector of dim 'm' of the upper-bound tension of each
        muscle.

    Returns:
        (np.array): upper bound of the muscle actuations that obey the linear
            system.
    """
    return ub - compute_muscle_actuation(wgt, tTilde, svdProc)


# ------------------------------------------------------------------------------
def plotWrench(w):
    """3D-plot the force and moment components of a wrench.

    Args:
        w (np.array): a wrench (dim = 6). First 3 values are the force
            component, least 3 are the moment component.

    Returns:
        ax (pyplot.Axis): 3d-axis object to access plot features
    """
    ax = plt.figure().add_subplot(projection='3d')

    ax.plot([0., w[0]], [0., w[1]], [0., w[2]], 'g-', label='force (N)')
    ax.plot([0., w[3]], [0., w[4]], [0., w[5]], 'b:', label='moment (N.m)')
    ax.plot(0., 0.,  0., 'ro')
    ax.text(0., 0., 0.05, 'P')
    ax.set_xlabel(r'X_j')
    ax.set_ylabel(r'Y_j')
    ax.set_zlabel(r'Z_j')
    ax.legend(loc="best")
    return ax


# =====================================
# CLASSES
# =====================================
class AbstractAnalysis(ABC):
    """An abstract class to handle musculoskeletal analyses.

    Derived classes must set analysed coordinates before launching the
    overloaded 'runAnalysis' method.
    """

    # --------------------------------------------------------------------------
    def __init__(self, prc=None):
        """Instantiate the abstract analysis class.

        Args:
            prc (MSM_Processor): interface to an OpenSim Model
        """
        if not prc:
            prc = MSM_Processor.MSM_Processor()
        self.prc = prc
        self.MMArm = None       # Matrix of moment arms
        self.J = None           # Jacobian matrix
        self.anlDoFsIds = None     # list of the ids of the analysed DoFS

        self.Z = None           # Zonotope

        self.ub = None          # vector of muscle tensions' upper bounds
        self.lb = None          # vector of muscle tensions' lower bounds
        lgg.debug("Empty abstract analysis created...")

    # --------------------------------------------------------------------------
    def setAnalysedDofs(self, dofNames):
        """Set the list of the analysed DoFs (determine the Zonotope' size).

        Args:
            dofIds (list): the ids of the dofs analysed
        """
        for nm in dofNames:
            assert self.prc.getCoordByName(nm) is not None, \
                "Wrong coordinate name %s" % nm
        self.anlDoFsIds = [self.prc.getCoordIdByName(nm) for nm in dofNames]
        lgg.debug("Ids of analysed DoFs : %s", self.anlDoFsIds)

    # --------------------------------------------------------------------------
    def setModelConfiguration(self, q, qdot):
        """Set model's joint angles, velocities, jacobian and moment arm matrix.

        Args:
            pos (np.array(nDof, float): joint angles of the whole model
            vel (np.array(nDof), float): joint velocities of the whole model
        """
        lgg.info("Set the model's posture...")
        nc = len(self.prc.getCoordinateNameList())
        assert isinstance(q, np.ndarray) and len(q) == nc, \
            "Invalid joint angles"
        assert isinstance(qdot, np.ndarray) and len(qdot) == nc, \
            "Invalid joint velocities"

        self.prc.setModelConfiguration(q, qdot)
        self.prc.show()
        lgg.info("Joint angles and velocities are now configured !")

        self.MMArm = self.prc.getMomentArmMatrix()
        lgg.info("Moment Arm Matrix updated !")

    # --------------------------------------------------------------------------
    def setModelPerformances(self):
        """Compute model's upper and lower bounds of muscle tensions."""
        lgg.info("Computing the hypercube of min-max muscle tensions...")
        nbm = self.prc.getNumMuscles()
        fMin = np.zeros(nbm)
        fMax = np.ones(nbm)
        for i in range(nbm):
            perf = self.prc.getMusclePerformance(self.prc.getMuscleById(i),
                                                 MSM_Processor.MAX_MSCL_ACT)
            fMax[i], fMin[i] = perf[0:2]  # other perf values are unused

        # lgg.debug("Min muscle tensions:\n%s", fMin)
        # lgg.debug("Max muscle tensions:\n%s", fMax)
        self.ub = fMax
        self.lb = fMin
        lgg.info("Muscle performances are now configured !")

    # --------------------------------------------------------------------------
    def getListOfActuatingMuscles(self):
        """Get the list of the indices of muscles which actuate any DoF.

        Returns the indices of muscles which actuate the DoFs i.e with
        "efficient" moment arms (greater than DEF_MARM_THR threshold).

        Returns:
            (list): the list of the ids of actuating muscles.
        """
        tCol = np.any(np.abs(self.MMArm[self.anlDoFsIds, :]) > DEF_MARM_THR,
                      axis=0)
        return np.nonzero(tCol)[0].tolist()

    # --------------------------------------------------------------------------
    def getTorqueZonotope(self):
        """Compute the torque zonotope.

        Compute the Zonotope determined by:
        - the selected coordinates' ;
        - the actuationg muscles in the current posture.

        Returns:
            (ZNT_Processor.Zonotope)
        """
        lgg.info("Computing the partial zonotope of needed DoFs and muscles")
        # select muscles which actuate the analysed DoFs
        selMuscles = \
            self.getListOfActuatingMuscles()
        lgg.debug("List of actuating muscles: %s", selMuscles)

        # submatrix extraction
        subM = self.MMArm[np.ix_(self.anlDoFsIds, selMuscles)]
        # lgg.debug("Moment-arm submatrix:\n%s", subM)

        selCrdNames = [self.prc.getCoordById(k).getName()
                       for k in self.anlDoFsIds]

        Z = ZNT_Processor.Zonotope(self.lb[selMuscles],
                                   self.ub[selMuscles],
                                   subM,
                                   selCrdNames)
        return Z

    # --------------------------------------------------------------------------
    def getActuationDirectionFromWrench(self, w_des):
        """Compute the unit vector which create the desired wrench.

        Args:
            w_des (np.array(6, floats): the desired wrench (force component
                first, then moment component)

        Returns:
            (np.array(floats)): the unit vector which creates a wrench
                colinear to the desired wrench.
        """
        # compute the objective set point
        # -------------------------------
        # wrench with only just force, no torque,
        # and magnitude = 100 (N)
        lgg.debug("Set-point wrench:\n%s", w_des)
        tau = self.J[:, self.anlDoFsIds].T @ w_des
        return tau/np.linalg.norm(tau)

    # --------------------------------------------------------------------------
    def getTensionsFromTorques(self, tau_des):
        """Compute muscles tensions from desired actuation torques.

        Args:
            tau_des (np.array(float)): the desired actuation torques
            withPlot (Boolean): ask for the bar plot of muscle tensions

        Returns:
            (np.array(floats)): the individual muscle tensions
        """
        # This function estimates a solution 't_' that solves the linear equation
        # \tau = N.t
        # with t_ subject to constraints : t_ > t_min et t_ < t_max,
        # based on SVD decomposition.


        # Initialize SVD problem
        selMusclIds = self.getListOfActuatingMuscles()
        subM = self.MMArm[np.ix_(self.anlDoFsIds, selMusclIds)]
        n, m = subM.shape
        mySVD = SVD_Processor.SVD_Solver(subM)

        # Round at the nearest interger to help the solving (avoid no solution)
        np.around(tau_des)

        # solve the linear equation : find one particular solution
        tTilde = mySVD.computeBasicSolution(tau_des)

        # define constraints
        lbConstr = {'type': 'ineq',
                    'fun': lb_constr_fun,
                    'args': (tTilde, mySVD, self.lb[selMusclIds])}
        ubConstr = {'type': 'ineq',
                    'fun': ub_constr_fun,
                    'args': (tTilde, mySVD, self.ub[selMusclIds])}

        # search for optimal muscle tensions
        # ----------------------------------
        lgg.debug("Trying to solve the constrained linear problem...\n")
        res = minimize(muscleTensions_cost_fun,
                       x0=np.random.randn(m-n),
                       args=(tTilde, mySVD),
                       constraints=(lbConstr, ubConstr))
        if res.success:
            lgg.info("Optimization result:\n%s", res)
            tBar = tTilde + mySVD.computeKernelLinComb(res.x)
            lgg.info("Optimal solution: %s", tBar)
            lgg.info("MMarm . tBar =\n%s", subM @ tBar)
            lgg.info("b_des =\n%s", tau_des)
        else:
            lgg.warning("No solution found for desired actuation %s", tau_des)
            tBar = None

        return tBar

    # --------------------------------------------------------------------------
    def updateModel(self, q, qdot):
        """Update the posture of the model and all associated features.
        
        Args:
            q (np.array): joint angles (dim = number of DoFs)
            qdot (np.array): joint velocities (dim = number of DoFs)
        """
        lgg.info("Updating the posture of the model...")
        self.setModelConfiguration(q, qdot)
        self.setModelPerformances()
        self.Z = self.getTorqueZonotope()

    # --------------------------------------------------------------------------
    def computeMaxT_AngleVel(self, hplan, withPlot=False):
        """Compute max actuation torque at the current posture.

        Compute max actuation torques at the current posture as the
        intersection between the torque zonotope and the experimental
        hyperplan.

        Args:
            hplan (list): description of the normalized hyperplan in the
                joint-torque space; hplan = [a_x, b_y, c_z, ... , d] with
                [a_x, b_y, c_z, ...] a unit vector normal to this hyperplan
                and d the distance to the origin.
            withPlot (boolean): ask for the plot of the zonotope and the
                intersection points

        Returns:
            (list): the intersection points.
        """
        # compute intersection with the given hyperplane
        intr = self.Z.getIntersectionPointsWithPlane(np.array(hplan))

        return intr

    # --------------------------------------------------------------------------
    def getMaxActuationInDirection(self, u):
        """Get the maximum actuation torques in the desired direction.

        Args:
            u (np.array(3)): vector giving the direction of the expected force,
                            expressed in the same frame as the jacobian
                            (point P of body j)

        Returns:
            (np.array): the coordinates of the maximum exertable torques
        """
        return self.Z.getIntersectionPointWithRay(u)

    # --------------------------------------------------------------------------
    def computeWrenchFromMuscleTensions(self, t):
        """Compute the wrench at point corresponding to the current Jacobian.

        Args:
            t (np.array): individual muscle tensions (in Newton)

        Returns:
            np.array(6): W = [ F, Gamma ] in IR ^ (6x1) the wrench (in the frame
                and at the point defined by the jacobian), F is the force
                component of W and Gamma its moment component
        """
        assert len(t) == self.prc.getNumMuscles(), "Wrong data size"
        W = np.zeros(6)
        N = self.MMArm
        y_des = N @ t
        lgg.debug("Torque actuation:\n%s (N.m)", y_des)
        Jdag = np.linalg.pinv(self.J.T)
        lgg.debug("Jdag:\n%s", Jdag)
        W = Jdag @ y_des
        lgg.debug("Wrench:\n%s", W)
        lgg.debug("Norm of resulting force  F=%.4e", np.linalg.norm(W[3:]))
        lgg.debug("Norm of resulting moment G=%.4e", np.linalg.norm(W[:-3]))
        return W

    # --------------------------------------------------------------------------
    def getNorm(self, vec):
        """Return the norm of the actuation vector.

        Return the norm of the actuation vector, or NaN if vec is None

        Args:
            vec (np.array): actuation vector

        Returns:
            (float): norm(vec) or NaN
        """
        return np.NaN if vec is None else np.linalg.norm(vec)

    # --------------------------------------------------------------------------
    @abstractmethod
    def run(self):
        """Abstract method to be overloaded.

        Implement the specific steps and computation of the analysis.
        """
        print("Run requested analysis")


# ===============================================================================
class JointTVA(AbstractAnalysis):
    """Run a 1-joint many-muscles TVA analysis."""

    # --------------------------------------------------------------------------
    def __init__(self, prc=None, cname=None):
        """Initialise the instance.

        Args:
            prc (MSM_Processor): the interface with the musculoskeletal model
            cname (string): the name of the Joint to be analyzed.
         """
        super().__init__(prc)
        if cname is None:
            cname = GUI_Manager.chooseCoordinateForTVA(self.prc)
            assert cname, "No joint selected for TVA, skipping..."
        self.cname = cname

    # --------------------------------------------------------------------------
    def run(self):
        """Compute the all-muscles 1-joint TVA."""

        # Choose the considered DOFS
        self.setAnalysedDofs(ANL_DOFS)
            # considering the elbow flexion-extension and fore-arm prono-supination
        hplan = [0., 1., 0.]
            # the projection plan is "zero prono-supination actuation torque

        # get the index of elbow angle
        k = self.prc.getCoordIdByName(self.cname)

        # initialize the array containing simulated max actuation torques
        nc = len(self.prc.getCoordinateNameList())
        qdot = np.zeros(nc)
        mact = np.zeros((len(DEF_Q_RANGE), len(DEF_V_RANGE), 2))
        # raws are for joint angle, columns are for joint velocities,
        # 3rd dimension is for agonist / antagonist max exertion

        # =================================
        # Build the matrix of max exertion
        # (loop on angles and velocities)
        # =================================
        for i, alpha in enumerate(DEF_Q_RANGE):
            q = np.zeros(nc)
            for j, nu in enumerate(DEF_V_RANGE):
                qdot = np.zeros(nc)
                q[k] = alpha
                qdot[k] = nu

                # update the Zonotope according to the current model's state
                self.updateModel(q, qdot)

                # get the max flexion/extension torques with zero pronosupination torque
                extr = np.array(self.computeMaxT_AngleVel(hplan, True))

                mact[i, j, :] = np.sort(extr[:, 0])

        # ==============================
        # Plot
        # ==============================
        # fig = plt.figure()
        # fig.suptitle("TVA surface for DoF: '%s'" % self.cname)
        # X, Y = np.meshgrid(DEF_Q_RANGE, DEF_V_RANGE, indexing='ij')

        # # Flexion
        # ax3d = fig.add_subplot(121, projection="3d")
        # ax3d.set_xlabel('Joint angle (deg)')
        # ax3d.set_ylabel('Joint velocity (deg/s)')
        # ax3d.set_zlabel('Flexion max torque (N.m)')
        # ax3d.plot_surface(X, Y, mact[:, :, 1],
                          # cmap=cm.coolwarm,
                          # linewidth=0.1, antialiased=False, alpha=0.75)
        # xmin = np.min(ax3d.axes.get_xlim3d())
        # ymax = np.max(ax3d.axes.get_ylim3d())
        # ax3d.contour(X, Y, mact[:, :, 1],
                     # zdir='x', offset=xmin, cmap='copper')
        # ax3d.contour(X, Y, mact[:, :, 1],
                     # zdir='y', offset=ymax, cmap='copper')

        # # Extension
        # ax3d = fig.add_subplot(122, projection="3d")
        # ax3d.set_xlabel('Joint angle (deg)')
        # ax3d.set_ylabel('Joint velocity (deg/s)')
        # ax3d.set_zlabel('Extension max torque (N.m)')
        # ax3d.plot_surface(X, Y, -mact[:, :, 0],
                          # cmap=cm.coolwarm,
                          # linewidth=0.1, antialiased=False, alpha=0.75)
        # xmin = np.min(ax3d.axes.get_xlim3d())
        # ymax = np.max(ax3d.axes.get_ylim3d())
        # ax3d.contour(X, Y, -mact[:, :, 0],
                     # zdir='x', offset=xmin, cmap='copper')
        # ax3d.contour(X, Y, -mact[:, :, 0],
                     # zdir='y', offset=ymax, cmap='copper')

        #=== In one subplot ===
        fig = plt.figure()
        fig.suptitle("TVA surface for DoF: '%s'" % self.cname)
        X, Y = np.meshgrid(DEF_Q_RANGE, DEF_V_RANGE, indexing='ij')

        min_, max_ = mact.min(), mact.max()

        ax3d = fig.add_subplot(111, projection="3d")
        ax3d.set_xlabel('Joint angle (deg)')
        ax3d.set_ylabel('Joint velocity (deg/s)')
        ax3d.set_zlabel('max torque (N.m)')

        flx = ax3d.plot_surface(X, Y, mact[:, :, 1],
                                cmap=cm.coolwarm,
                                vmin=min_, vmax=max_,
                                linewidth=0.1, antialiased=False, alpha=0.75)
        ext = ax3d.plot_surface(X, Y, mact[:, :, 0],
                                cmap=cm.coolwarm,
                                vmin=min_, vmax=max_,
                                linewidth=0.1, antialiased=False, alpha=0.75)
        fig.colorbar(flx, location='left').set_label('Flexion(+) and Extension (-)',
                                                     rotation=90)

        xmin = np.min(ax3d.axes.get_xlim3d())
        ymax = np.max(ax3d.axes.get_ylim3d())

        ax3d.contour(X, Y, mact[:, :, 1],
                     zdir='x', offset=xmin, cmap='copper')
        ax3d.contour(X, Y, mact[:, :, 1],
                     zdir='y', offset=ymax, cmap='copper')
        ax3d.contour(X, Y, mact[:, :, 0],
                     zdir='x', offset=xmin, cmap='copper')
        ax3d.contour(X, Y, mact[:, :, 0],
                     zdir='y', offset=ymax, cmap='copper')

        plt.show()


# =====================================
# MAIN CODE
# =====================================
def main():
    """Test the module."""
    # Initialize needed objects
    # -------------------------
    lgg.info("Test the module...")
    prc = MSM_Processor.MSM_Processor(MOBL_MDL_FILE, True)
    lgg.info("with MoBL's model")
    anl = JointTVA(prc, "elbow_flexion")
    anl.run()

    lgg.info("Done !")


# =====================================
if __name__ == "__main__":
    # logging stuff
    logging.config.fileConfig(DEF_LOG_CONF, disable_existing_loggers=False)
    # 'False' so that other loggers are still displayed

    logging.getLogger('matplotlib.font_manager').setLevel(logging.WARNING)
    # disable matplotlib font logs

    lgg = logging.getLogger(DEF_LOGGER)
    main()
