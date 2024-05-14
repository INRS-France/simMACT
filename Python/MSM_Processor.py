# -*- coding: utf-8 -*-
"""Process essential features to run simulate maximum actuation torques."""
# =============================================================================
# MODULES
# =============================================================================
from os.path import join, dirname, abspath
import tkinter.filedialog
import tkinter
import os
import logging
import tqdm  # progress bar
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps as cmap
import opensim as osim
# import pandas as pd

lgg = logging.getLogger(__name__)

# from mpl_toolkits import mplot3d
# from functools import partial  # used to run a tk callback with parameter

# =============================================================================
# CONSTANTS
# =============================================================================
DEF_MDL_DIR = join(dirname(dirname(abspath(__file__))), "Models/")
OSIM_MDL_FILE = join(DEF_MDL_DIR, "arm26/arm26.osim")
MOBL_MDL_FILE = join(DEF_MDL_DIR, "MOBL/4.1/MOBL_ARMS_fixed_41.osim")
MOBL_MDL_IGN_FILE = join(DEF_MDL_DIR, "MOBL/4.1/MOBL_ARMS_fixed_41-ignoreTendonCompl.osim")
DEF_PLOT_DIR = join(dirname(dirname(dirname(abspath(__file__)))),
                    "simMACT-WIP/Plots")

DEF_MUSCLE = "BICshort"
DEF_BODY = "humerus"
DEF_NB_JNT_ANG = 8      # Default number of joint angles to plot for each
ALL_ELB_FLX_MUSCLES = ["BIClong", "BICshort", "BRA"]
DEF_COORD = "r_elbow_flex"   # "elbow_flexion"
DEF_ELB_ANG = 90.   # default albow flexion angle
MAX_MSCL_ACT = 1.0
MIN_MSCL_ACT = 0.0  # 0.05 for some implementation of Hill's model /!\

DEF_Q_RANGE = np.arange(0., 141., 10.)   # joint angles (deg)
DEF_V_RANGE = np.arange(-200., 201., 50.)   # joint velocities (deg/s)
NULL_VEL_IND = 5    # index at null velocity

# default moment arm threshold for agonist / antagonist / indiff.
DEF_MOMENT_ARM_THR = 1e-5

RAD = np.radians
DEG = np.degrees

GUI_COL_SZ = 10  # number of item per column in GUI windows


# =============================================================================
# CLASS MSM_Processor
# =============================================================================
class MSM_Processor():
    """Basic features to process OpenSim musculoskeletal models.

    This class computes the basic features (moment arm, muscle forces, jacobian
    matrices) needed for advanced computations as Torque-Velocity-Angle
    characteristics and maximum actuation torques.
    """

    # --------------------------------------------------------------------------
    def __init__(self, osimFile="", withVisu=False):
        """Read the given OSIM model file or ask the user to choose one.

        Args:
            osimFile (string): absolute path of the OpenSim model file
                (optional)
            withVisu (boolean): whether to create an OpenSim viewer window
                (optional)

        """
        # get model path
        if not os.path.exists(osimFile):
            lgg.warning("osim file '%s' not found, asking for another !",
                        osimFile)
            root = tkinter.Tk()
            opt = {}
            opt['initialdir'] = '.'
            opt['title'] = u"Choose a valid OpenSim model (osim) :"
            opt['filetypes'] = [("OpenSim model", ("*.osim"))]
            osimFile = str(tkinter.filedialog.askopenfilename(**opt))
            root.destroy()

        # load model and initialize the class
        lgg.debug("Processing %s" % osimFile)
        self.model = osim.Model(osimFile)
        self.model.setUseVisualizer(withVisu)
        self.state = self.model.initSystem()
        if withVisu:
            self.visu = self.model.getVisualizer()
            simtkVisu = self.visu.getSimbodyVisualizer()
#            simtkVisu.setShutdownWhenDestructed(True)
            if self.model.getName() == "Right":   # MoBL
                simtkVisu.setGroundHeight(-1.0)   # ground is at CLAV height !
        else:
            self.visu = None
        self.analysedCoord = None

#        self.grvHlp = osim.GravityHelper(self.model)

        lgg.debug("Ready to process model '%s'", self.model.getName())
        lgg.debug("Nb of muscles: %d", self.model.getMuscles().getSize())
        lgg.debug("Nb of DDL:     %d", self.model.getNumCoordinates())
        lgg.debug("Nb of Joints:  %d", self.model.getNumJoints())

    # # --------------------------------------------------------------------------
    # def __del__(self):
    #     """Delete the object.
    #
    #     Shut the SimTK::Visualizer down, if any.
    #     """
    #     if self.visu is not None:
    #         lgg.debug("Trying to kill the SimTK::Visualizer")
    #         self.visu.getSimbodyVisualizer().shutdown()

# ------------------------------------------------------------------------------
    def chooseItemInModel(self, prompt, listGetter):
        """GUI dialog to choose one coordinate in the model.

        Args:
            prompt (string): the title of the window
            listGetter (function): returns the list of names of 'item' in
                       the model

        Return:
            string: name of the chosen item
        """
        lgg.debug("Asking the user to choose one item...")

        root = tkinter.Tk()
        root.title(prompt)

        frm = tkinter.LabelFrame(root, text='Available items:',
                                 padx=10, pady=10)
        frm.grid(column=0, row=0, padx=10, pady=10)

        # create choice radio buttons
        items = listGetter()
        inVar = tkinter.StringVar(root, "")

        for i, iname in enumerate(items):
            rb = tkinter.Radiobutton(frm,
                                     text=iname,
                                     variable=inVar,
                                     value=iname,
                                     tristatevalue="x")
            rb.grid(column=i // GUI_COL_SZ,
                    row=i % GUI_COL_SZ,
                    ipadx=10, ipady=10)
            # 'tristate' to prevent buttons to be activated by default

        # create OK button
        bOK = tkinter.Button(root, text="OK", command=root.destroy)
        bOK.grid(column=0, row=1, padx=10, pady=10)

        # compute the size of the main window
        x_size = 2 * ((i+1) // GUI_COL_SZ)
        # number of columns in frm1 and frm2
        y_size = GUI_COL_SZ             # max number of rows
        margin = 10
        x_scale = 300                    # approximated width  of a column
        y_scale = 50                     # approximated height of a row

        sz = ""
        sz = str(x_size*x_scale + 2*margin) + "x" + \
            str((2+y_size)*y_scale + 2*margin)
        root.geometry(sz)
        root.mainloop()

        return inVar.get()

    # --------------------------------------------------------------------------
    def getNumMuscles(self):
        """Return the number of muscles in the model."""
        return self.model.getMuscles().getSize()

    # --------------------------------------------------------------------------
    def getMusclesList(self):
        """Return the list (osim.Set) of the muscles in the model."""
        ml = self.model.getMuscleList()
        return ml

    # --------------------------------------------------------------------------
    def getJointsList(self):
        """Return the list (osim.Set) of the joints in the model."""
        return self.model.getJointList()

    # --------------------------------------------------------------------------
    def getNbCoordinates(self):
        """Return the number of coordinates (DoF) in the model."""
        return self.model.getNumCoordinates()

    # --------------------------------------------------------------------------
    def getMuscleNameList(self):
        """Return the names of the muscles in the model."""
        ret = []
        for m in self.model.getMuscles():
            ret.append(m.getName())
        return ret

    # --------------------------------------------------------------------------
    def getCoordinateNameList(self):
        """Return the names of the coordinates (DoF) in the model."""
        ret = []
        for c in self.model.getCoordinateSet():
            ret.append(c.getName())
        return ret

    # --------------------------------------------------------------------------
    def getCoordByName(self, coordName):
        """Return the osim.Coordinate object according to its name.

        Args:
            coordName (string): the name of the Coordinate
        """
        return self.model.getCoordinateSet().get(coordName)

    # --------------------------------------------------------------------------
    def getCoordById(self, c_id):
        """Return the osim.Coordinate object according to its id.

        Args:
            c_id (int): the id of the Coordinate

        Return:
            (opensim.Coordinate): the Coordinate object
        """
        return self.model.getCoordinateSet().get(c_id)

    # --------------------------------------------------------------------------
    def getCoordIdByName(self, coordName):
        """Return the Id of the osim.Coordinate object according to its name.

        Args:
            coordName : string, the name of the Coordinate

        Return:
            int: the id of the Coordinate object
        """
        return self.model.getCoordinateSet().getIndex(coordName)

    # --------------------------------------------------------------------------
    def getJointsAngles(self):
        """Return the angle values of all coordinates.

        Return:
            list: the values of all joint angles.
        """
        s = self.state
        return [DEG(c.getValue(s)) for c in self.model.getCoordinateSet()]

    # --------------------------------------------------------------------------
    def getJointsVelocities(self):
        """Return the velocity values of all coordinates..

        Return:
            list: the values of all joint angle velocities.
        """
        s = self.state
        return [DEG(c.getSpeedValue(s)) for c in self.model.getCoordinateSet()]

    # --------------------------------------------------------------------------
    def getMuscleByName(self, muscleName):
        """Return the osim.Muscle object according to its name.

        Args:
            muscleName : string, the name of the expected muscle

        Return:
            opensim.Muscle: the Muscle object.
        """
        return self.model.getMuscles().get(muscleName)

    # --------------------------------------------------------------------------
    def getMuscleById(self, muscleId):
        """Return the osim.Muscle object corresponding with its id.

        Args:
            muscleId : int, the id of the muscle.

        Return:
            opensim.Muscle: the Muscle object.
        """
        return self.model.getMuscles().get(int(muscleId))
        # cast to int in case mi is a np.int64

    # --------------------------------------------------------------------------
    def getMusclesRoleForCoord(self, coord):
        """Identify agonist, antagonist and indifferent muscles.

        Args:
            coord : osim.Coordinate, the studied joint.

        Return:
            agon : np.array, the indices of agonist muscles
            antg : np.array), the indices of antagonist muscles
            indf : np.array), the indices of non actuating muscles
        """
        mav = self.getMomentArmVector(coord)
        agon = np.where(mav > DEF_MOMENT_ARM_THR)[0]
        antg = np.where(mav < -DEF_MOMENT_ARM_THR)[0]
        indf = np.where(np.abs(mav) < DEF_MOMENT_ARM_THR)[0]

        mn = [m.getName() for m in self.model.getMuscleList()]

        lgg.debug("DoF %s:\nagonists = %s\nantagonists = %s\nindifferent = %s",
                  coord.getName(),
                  [mn[i] for i in agon],
                  [mn[i] for i in antg],
                  [mn[i] for i in indf])

        return agon, antg, indf

    # --------------------------------------------------------------------------
    def getMomentArmVector(self, coord):
        """Return the moment arms for all muscle for the given coordinate.

        Args:
            coord : osim.Coordinate, the studied joint.

        Return:
            np.array: the vector of moment arms for all muscles
        """
        msclList = self.getMusclesList()
        v = np.zeros(len(msclList))
        for i, m in enumerate(msclList):
            v[i] = m.computeMomentArm(self.state, coord)
        return v

    # --------------------------------------------------------------------------
    def getMomentArmMatrix(self):
        """Return the nxm matrix of moment arms.

        Return:
            np.array(n,m): the matrix of moment arms for all 'n' coordinates
                and all 'm' muscles of the model.
        """
        cl = self.model.getCoordinateSet()
        ml = self.model.getMuscles()
        M = np.zeros((cl.getSize(), ml.getSize()), float)
        for j, m in enumerate(ml):
            for i, c in enumerate(cl):
                M[i, j] = m.computeMomentArm(self.state, c)
        return M

    # --------------------------------------------------------------------------
    def getMusclePerformance(self, mscl, act):
        """Return the basic forces exerted by the given muscle and activation.

        Args:
            mscl : osim.Muscle, the studied muscle
            act : float in [0:1], its activation level

        Return:
            fa: float, active fiber force
            fp: float, passive fiber force
            fl: float, fiber length
            nfl: float, normalized fiber length
        """
        # Set the activation of the desired muscle, don't change its fiber
        # length as it is deduced from the current configuration (vq, vqdot)
        # of the whole model
        mscl.setActivation(self.state, act)
        self.model.equilibrateMuscles(self.state)

        fa = mscl.getActiveFiberForce(self.state)      # Active fiber force
        fp = mscl.getPassiveFiberForce(self.state)     # Passive fiber force
        fl = mscl.getFiberLength(self.state)           # Fiber length
        nfl = mscl.getNormalizedFiberLength(
            self.state)  # Normalized fiber length

        return fa, fp, fl, nfl

    # --------------------------------------------------------------------------
    def setCoordinateAngVel(self, coord, q, qdot):
        """Set the configuration (angle, velocity) of one coordinate.

        Set the configuration (angle, velocity) of one coordinate and equilibrate
        the Model.

        Args:
            coord: opensim.Coordinate: the joint to configure
            q: double, the joint angle (deg)
            qdot: double, the joint velocity (deg/s)

        Return:
            tuple, the updated angle and velocity.
        """
        coord.setValue(self.state, RAD(q))
        coord.setSpeedValue(self.state, RAD(qdot))

        # Equilibrate the muscle and tendon forces.
        self.model.equilibrateMuscles(self.state)

        return DEG(coord.getValue(self.state)), \
            DEG(coord.getSpeedValue(self.state))

    # --------------------------------------------------------------------------
    def resetModelConfiguration(self):
        """Set the angles and velocities to 0 for all joints."""
        lgg.debug("Reset joint angles and velocities of the model")
        for c in self.model.getCoordinateSet():
            c.setValue(self.state, 0.)
            c.setSpeedValue(self.state, 0.)

        # Equilibrate the muscle and tendon forces.
        self.model.equilibrateMuscles(self.state)

    # --------------------------------------------------------------------------
    def setModelConfiguration(self, vq, vqd):
        """Set the joint configuration (angles, velocities) of the model.

        Also equilibrate the system (solve muscle lengths).

        Args:
            vq: np.array, vector of doubles, joint angles (deg)
            vqd: np.array, vector of doubles, joint velocities (deg/s)
        """
        coordSet = self.model.getCoordinateSet()  # set of coordinates

        for i, c in enumerate(coordSet):
            self.setCoordinateAngVel(c, vq[i], vqd[i])
            # lgg.debug("/%s/%s: (%02.2f °;%02.2f °/s)",
            #           c.getJoint().getName(),  # joint
            #           c.getName(),             # coordinate
            #           DEG(c.getValue(self.state)),
            #           DEG(c.getSpeedValue(self.state)))

        # Equilibrate the muscle and tendon forces.
        self.model.equilibrateMuscles(self.state)

    # --------------------------------------------------------------------------
    def getBodyList(self):
        """Return the list of bodies in the model."""
        return self.model.getBodySet()

    # --------------------------------------------------------------------------
    def getBodyNameList(self):
        """Return the list of the bodies'name in the model.

        Return:
            list, body names"""
        return [b.getName() for b in self.model.getBodySet()]

    # --------------------------------------------------------------------------
    def getBodyByName(self, bdname):
        """Return the body object from its name.

        Args:
            bdname: string, the name of the body.

        Return:
            osim.Body object
        """
        return self.model.getBodySet().get(bdname)

    # --------------------------------------------------------------------------
    def getBodyIdByName(self, bdname):
        """Return the Id of the osim.Body object according to its name.

        Args:
            bdname : string, the name of the Body

        Return:
            int: the id of the Body object
        """
        return self.model.getBodySet().getIndex(bdname)


    # --------------------------------------------------------------------------
    def getBodyRotation(self, name):
        """Return the Rotation matrix of one body from its name.

        simbody types (Rotation, Mat33, Vec3) are a mess => convert to numpy
        ndarray.

        Return:
            np.ndarray, the rotation matrix.
        """
        b = self.getBodyByName(name)
        R = b.getTransformInGround(self.state).R()
        nc = R.ncol()
        nr = R.nrow()
        npR = np.zeros((nr, nc))
        for i in range(nr):
            for j in range(nc):
                npR[i, j] = R.get(i, j)
        return npR

    # # --------------------------------------------------------------------------
    # def getGravityForces(self):
        # """Return the gravity force/torque applied on each body.

        # This function relies on the SimTK::Force::Gravity::getBodyForces()
        # WHICH BY DEFAULT IS NOT EXPOSED by SWIG (2023-03-02 : pending pulling
        # request #3416) so it may require a re-compilation of OpenSim's Python
        # wrapping or th access to a git artifact.

        # Return:
            # simTk.Vector_<simTk.SpatialVec>: vector of gravity forces
                # on all bodies, indexed by simTk.MobilizedBodyIndex.
        # """
        # return self.grvHlp.getGravityForces(self.state)

    # # --------------------------------------------------------------------------
    # def getGravityForcesOnBody(self, bdname):
        # """Return the gravity force/torque applied on the specified body.

        # The returned value is a pair of 3D-vectors, first the gravity moments
        # about the **body origin** (not its CoM), then the resulting gravity
        # forces.

        # Return:
            # simTk.SpatialVec>: vector of gravity forces applied on the
                # desired body, None if body not found.
        # """
        # spvec = None
        # try:
            # idx = self.getBodyIdByName(bdname)
            # lgg.debug("Body %s is at index %d", bdname, idx)
            # spvec = self.grvHlp.getGravityForces(self.state).get(idx)
        # except Exception as exc:
            # lgg.warning("Body %s not found, %s", bdname, exc)
        # return spvec

    # --------------------------------------------------------------------------
    def getJacobian(self, pos, bname):
        """Compute the jacobian matrix at position 'pos' of body 'bname'.

        Args:
            pos: np.array(3), cartesian position (in the current body frame)
                where the jacobian is to be computed
            bname: string, name of the body segment

        Return:
            jac: np.array(6, nddl), the jacobian matrix. First 3 raws
                correspond to translation (force part of the wrench), last 3
                raws correspond to rotation (torque part of the wrench)

        TODO : test if the following instructions lead to the same result
        jm = osim.Matrix()
        ms = osim.Model.getMatterSubsystem()
        ms.calcSystemJacobian(state, jm)
        jm =
        """
        bs = self.getBodyList()
        assert bname in [b.getName() for b in bs], "Unknown body %s !" % bname

        nddl = self.model.getNumCoordinates()
        cs = self.model.getCoordinateSet()
        bd = bs.get(bname)

        v3Velocity = osim.Vec3(*np.zeros(3))  # '*' to pass parameters as list
        v3AngulVel = osim.Vec3(*np.zeros(3))  # idem
        v3Endeffec = osim.Vec3(*pos)          # idem
        # simTkEngine uses Vec3 type...

        jac = np.zeros((6, nddl))

        # Loop on DDLs to 'activate' only one DDL at a time
        # (build jacobian matrix columnwise)
        for i in range(nddl):
            qd = np.zeros(nddl)
            qd[i] = 1.

            # update speed of each joint
            for j in range(nddl):
                cs.get(j).setSpeedValue(self.state, qd[j])

            # solve kinematics state
            self.model.realizeVelocity(self.state)

            # compute velocities
            simTkEngine = self.model.getSimbodyEngine()
            simTkEngine.getVelocity(self.state, bd, v3Endeffec, (v3Velocity))
            simTkEngine.getAngularVelocity(self.state, bd, v3AngulVel)

            tVel = np.array(
                [v3Velocity.get(0), v3Velocity.get(1), v3Velocity.get(2)])
            rVel = np.array(
                [v3AngulVel.get(0), v3AngulVel.get(1), v3AngulVel.get(2)])
            jac[:, i] = np.concatenate((tVel, rVel))

            # https://simbody.github.io/simbody-3.6-doxygen/api/classSimTK_1_1SimbodyMatterSubsystem.html#aca444046ca1a2cae0c80b7846c29abf4
            # j_alt = osim.Matrix()
            # ms = self.model.getMatterSubsystem()
            # ms.calcStationJacobian(self.state,
            #                        bs.getIndex(bname),
            #                        osim.Vec3(*pos),
            #                        j_alt)
            #
            # ou bien ms.calcSystemJacobian(self.state, j_alt)) ??
            #
            # lgg.debug("simTK-calculated Jacobian:\n%s", j_alt.toString())
            # lgg.debug("NRg algorithm for Jacobian:\n%s", jac)

        return jac

    # --------------------------------------------------------------------------
    def plotMuscleCharacteristics(self, mName, cName):
        """Plot the main characteristics of 1 muscle and 1 joint.

        Args:
            mName: string, the name of one muscle
            cName: string, the name of one coordinate

        Return:
            plt.Axes, handle on the figure
        """
        assert mName in self.getMuscleNameList(), \
            "Invalid muscle name '%s'" % mName
        assert cName in self.getCoordinateNameList(), \
            "Invalid coordinate name ''%s'" % cName
        lgg.debug("Plotting characteristics of '%s' linked to joint '%s'",
                  mName, cName)

        # Initializations
        coord = self.getCoordByName(cName)
        coordId = self.getCoordIdByName(cName)
        mscl = self.getMuscleByName(mName)
        jntNb = self.getNbCoordinates()
        ja_lb = DEG(coord.getRangeMin())   # lower bound on joint angle
        ja_ub = DEG(coord.getRangeMax())   # upper bound on joint angle
        angles = np.linspace(ja_lb, ja_ub, DEF_NB_JNT_ANG, endpoint=True)
        X, Y = np.meshgrid(angles, DEF_V_RANGE, indexing='ij')

        ref_jntAng = np.zeros(jntNb)
        ref_jntAng[coordId] = DEF_ELB_ANG
        ref_jntVel = np.zeros(jntNb)

        mArm = np.zeros(DEF_NB_JNT_ANG)
        mslFtr = np.zeros((DEF_NB_JNT_ANG,
                          len(DEF_V_RANGE),
                          4))      # 4 values returned
        fbr_vel = np.zeros((DEF_NB_JNT_ANG, len(DEF_V_RANGE)))
        # fiber velocity
        fvm = np.zeros((DEF_NB_JNT_ANG, len(DEF_V_RANGE)))
        # force-velocity multiplier
        flm = np.zeros(DEF_NB_JNT_ANG)
        # force-length multiplier

        # try:
            # m2012eq = osim.Millard2012EquilibriumMuscle.safeDownCast(mscl)
            # mscl = m2012eq
        # except:
            # lgg.warning("%s is not a Millard2012Equilibrium instance", mName)
        # finally:
            # lgg.warning("%s is now a %s instance", mName, type(mscl))

        # =============================
        # Compute characteristics
        # =============================
#        for i, alpha in enumerate(angles):
        for i in tqdm.tqdm(range(DEF_NB_JNT_ANG),
                           desc="Muscle characteristics progression"):
            alpha = angles[i]

            # moment arm
            mArm[i] = mscl.computeMomentArm(self.state, coord)     # cm

            # force features
            flm[i] = mscl.getActiveForceLengthMultiplier(self.state)
            for j, vel in enumerate(DEF_V_RANGE):
                q = ref_jntAng
                qdot = ref_jntVel
                q[coordId] = alpha
                qdot[coordId] = vel
                self.setModelConfiguration(q, qdot)
                mslFtr[i, j, :] = self.getMusclePerformance(mscl, MAX_MSCL_ACT)
                fvm[i, j] = mscl.getForceVelocityMultiplier(self.state)
                fbr_vel[i, j] = mscl.getFiberVelocity(self.state)
            # visu
            self.show()

        AF = mslFtr[:, :, 0] + mslFtr[:, :, 1]      # active force
        MM = np.vstack((mArm,)*len(DEF_V_RANGE))    # duplicated moment arms
        T = AF * MM.T                               # moment arm

        # =============================
        # Plot characteristics
        # =============================
        mycm = cmap["coolwarm"]
        f = plt.figure(figsize=[19.2, 10.97])
        f.suptitle("Model '%s'\nCharacteristics of muscle '%s', joint '%s'" %
                   (self.model.getName(), mName, cName))

        # 2D plot of active, passive and total force vs q
        ax = f.add_subplot(231)
        ax.grid("on")
        ax.set_xlabel("Joint angle (°)")
        ax.set_ylabel("Isometric force (N)")
        ax.plot(angles, mslFtr[:, NULL_VEL_IND, 0], label="Active")
        ax.plot(angles, mslFtr[:, NULL_VEL_IND, 1], label="Passive")
        ax.plot(angles, mslFtr[:, NULL_VEL_IND, 0]
                + mslFtr[:, NULL_VEL_IND, 1], label="Total")
        ax.legend(loc="best")

        # 2D plot of force-length relation
        ax = f.add_subplot(232)
        ax.grid("on")
        ax.set_xlabel("Joint angle (°)")
        ax.set_ylabel("f-l multipl.")
        ax.plot(angles, flm)

        # 3D plot of force-velocity relation
        ax = f.add_subplot(233, projection="3d")
        ax.grid("on")
        ax.set_xlabel("Joint angle (°)")
        ax.set_ylabel("Joint velocity (°/s)")
        ax.set_zlabel("f-v multipl.")
        ax.plot_surface(X, Y, fvm,
                        cmap=mycm,
                        linewidth=0.1, antialiased=False, alpha=0.75)
        xmin = np.min(ax.axes.get_xlim3d())
        ymax = np.max(ax.axes.get_ylim3d())
        ax.contour(X, Y, fvm, zdir='x', offset=xmin, cmap='copper')
        ax.contour(X, Y, fvm, zdir='y', offset=ymax, cmap='copper')

        # 2D plot of normalized fiber length vs joint angle
        ax = f.add_subplot(234)
        ax.grid("on")
        ax.set_xlabel("Joint angle (°)")
        ax.set_ylabel("Normalized fiber length (%)")
        ax.plot(angles, mslFtr[:, NULL_VEL_IND, 3])

        # 2D plot of moment arm vs joint angle
        ax = f.add_subplot(235)
        ax.grid("on")
        ax.set_xlabel("Joint angle (°)")
        ax.set_ylabel("Muscle moment arm (cm)")
        ax.plot(angles, 100.*np.abs(mArm))

        # 3D plot of total torque vs q, qdot
        ax = f.add_subplot(236, projection='3d')
        ax.grid("on")
        ax.set_xlabel("Joint angle (°)")
        ax.set_ylabel("Joint velocity (°/s)")
        ax.set_zlabel("Joint torque (N.m)")
        ax.plot_surface(X, Y, T,
                        cmap=mycm,
                        linewidth=0.1, antialiased=False, alpha=0.75)
        xmin = np.min(ax.axes.get_xlim3d())
        ymax = np.max(ax.axes.get_ylim3d())
        ax.contour(X, Y, T, zdir='x', offset=xmin, cmap='copper')
        ax.contour(X, Y, T, zdir='y', offset=ymax, cmap='copper')

        return f

    # --------------------------------------------------------------------------
    def plotMuscleTensions(self,
                           tensions,
                           mIndList=None,
                           ax=None):
        """Horizontal bars plot of muscle tensions.

        Args:
            tensions: list(float) of size p <= nb of muscles in the model,
                    values of muscle tensions to be plotted
            mIndList: list(int) (optional) of size p<=nb of muscles in the model,
                    indices of muscle tensions to be plotted
            ax: plt.Axes, where to plot the Zonotope. If None, create the figure

        Return:
            plt.Axes (if figure was created here)
        """
        if mIndList is None:
            mIndList = range(self.getNumMuscles())
        if ax is None:
            f = plt.figure()
            f.suptitle("2D actuation zonotope")
            ax = f.add_subplot(111)
        else:
            assert isinstance(ax, plt.Axes), "Invalid axes"

        # preliminary tests
        nbm = len(mIndList)
        assert nbm <= self.getNumMuscles(), "Invalid list of muscles"
        assert max(mIndList) <= self.getNumMuscles(), "Invalid index of muscle"
        # assert nbm == len(tensions), "Invalid input data size"

        # build list of muscle names (for legend)
        mNames = [self.getMuscleById(k).getName() for k in mIndList]

        # compute tension upper and lower bound for each muscle
        fMin = np.zeros(nbm)
        fMax = np.ones(nbm)
        for i, k in enumerate(mIndList):
            fMax[i], fMin[i] = \
                self.getMusclePerformance(self.getMuscleById(k),
                                          MAX_MSCL_ACT)[0:2]

        # plot features
        ax.set_xlabel('Tension (N)')    # , fontsize=15)
        ax.set_ylabel('Muscle')         # , fontsize=15
        ax.grid("on")

        # use error bars to plot min/max tension
        err = np.abs(np.stack((tensions-fMin, fMax-tensions)))
            # there might be neglectable rounding errors which lead to negative differences
        ax.barh(mNames, tensions-fMin, xerr=err, left=fMin, capsize=2)

        return ax

    # --------------------------------------------------------------------------
    def show(self):
        """Show the model in the OpenSim visualisation window (if any)."""
        if self.visu:
            self.visu.show(self.state)


# === MAIN CODE ===============================================================
def main():
    """Test the class MSM_Processor."""
    lgg.info("Testing MSM_Processor class with default model")

    # initialize
    # ------------
#    myPrc = MSM_Processor(MOBL_MDL_IGN_FILE, True)
    myPrc = MSM_Processor(withVisu=True)

    # set the model to a test configuration
    # ---------------------------------------
    v_q = np.zeros(myPrc.getNbCoordinates())
    v_qdot = np.zeros(myPrc.getNbCoordinates())
    myPrc.setModelConfiguration(v_q, v_qdot)
    myPrc.show()

    # get gravity effects
    # geff = myPrc.getGravityForcesOnBody(DEF_BODY)
    # print(geff)



    # Plot Muscle characteristics (force/length, force-velocity)
    # for one muscle and one joint
    # --------------------------------------------------------
    mName = "BICshort"
    cName = "r_elbow_flex"          # Arm26
#    cName = "elbow_flexion"         # MOBL
    # mName = myPrc.chooseItemInModel("Choose muscle to study:",
    #                                 myPrc.getMuscleNameList)
    # cName = myPrc.chooseItemInModel(
    #     "Choose associated coordinate:", myPrc.getCoordinateNameList)
    myPrc.plotMuscleCharacteristics(mName, cName)
    # plt.savefig(join(DEF_PLOT_DIR,
    #                  "%s-%s-%s.png" % (myPrc.model.getName(), mName, cName)))
    plt.show()
    lgg.info("Test finished...")


# =============================================================================
if __name__ == "__main__":
    # logging stuff
    lgg.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        '%(name)s - %(levelname)s - %(filename)s:%(lineno)s\n%(message)s')
    ch.setFormatter(formatter)
    lgg.addHandler(ch)

    main()
