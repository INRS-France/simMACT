# -*- coding: utf-8 -*-
"""Basic freatures for zonotope processing in the field of biomechanics.

This class uses Krut & Goutterfarde's works and the  Hyperplane Shifting Method
to build zonotopes in their H-representation. It also uses the `polytope` python
package to get the V-representation of a zonotope (other package would do as well, 
for instance `pycddlib`, but the associated license implies more constraints upon 
simMACT diffusion.

JSN, 2022-04-15 - création
"""

# Several other Python modules for polytopes/zonotopes processing
#   -> https://www.geeksforgeeks.org/polytopes-in-python/
#   -> https://pypi.org/project/polytope/
#   -> https://github.com/sadraddini/pypolycontain
#   -> https://pypi.org/project/pycddlib/
# 
# More doc about the `polytope` package:
#   https://tulip-control.github.io/polytope/
#   https://tulip-control.sourceforge.io/doc/_modules/polytope/polytope.html

# plot 3D : only vertices and edges :
# https://stackoverflow.com/questions/27270477/3d-convex-hull-from-point-cloud


# =============================
# DEPENDENCIES
# =============================
import itertools        # automatic computation of combinations
import math
import logging
import logging.config
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

import polytope as pt   # primitives for polytope processing
import scipy.spatial as sp

logger = logging.getLogger(__name__)
# =============================
# CONSTANTS
# =============================
# default threshold for colinear unit vectors (dot product)
EPS = 1e-3
DEF_DOT_TOL = 1 - EPS

# /!\ For test purpose /!\
# A 2D-zonotope
DEF_M2D = np.array([[-0.0374, -0.,     -0.,      0.0216,  0.0074, -0.],
                    [-0.0212, -0.0212, -0.0212,  0.0405,  0.0405,  0.0167]])
    # moment arm values for model arm26.osim (2 joints, 6 muscles)
    # with q(°)=[0, 60] and qdot(°/s)= [0, 0].

# A 3D-zonotope
DEF_CUB = np.eye(3)
DEF_M3D = np.array([[1., 0.1, 0.1, -1., -0.1, 0.],
                    [0., 1.,  0.2,   0., -1.,  0.],
                    [0., 0.,  1.,   0.,  0., -1.]])

# Misc
PI = math.pi
DEF_PHI = PI/6  # angle in the horizontal plane
DEF_THETA = 1.1*PI/3  # angle with the vertical direction

# =============================
# GLOBAL FUNCTIONS
# =============================
def get_angle(v1, v2):
    """Compute the planar angle between 2 N-D vectors.

    Args:
        v1(np.array(): the 1st N-D vector
        v2(np.array(): the 2nd N-D vector

    Returns:
        (float) the value (in rad) of the planar angle between the 2 vectors.
    """
    u1 = v1 / np.linalg.norm(v1)
    u2 = v2 / np.linalg.norm(v2)

    y = u1 - u2
    x = u1 + u2

    a0 = 2 * np.arctan(np.linalg.norm(y) / np.linalg.norm(x))

    if (not np.signbit(a0)) or np.signbit(pi - a0):
        return a0
    elif np.signbit(a0):
        return 0.
    else:
        return np.pi

# --------------------------------------------------------------------------
def orderPointsInPlane(pts, n):
    """Order the points as path of consecutive points.
    
    Args:
        pts (np.array(nbP,N)): npP coplanar points in a N-D space
        n (np.array(N)): the N-D vector normal to the plane

    Returns:
        (list): the ordered list, made of the sorted nbP indices of lines
                of the 'pts' array. 
    """
    ordered = []
    nbP = pts.shape[0]
    angles = np.zeros(nbP)

    bc = np.mean(pts, axis=0)
    rays = pts - bc
    u = rays[0,:]
 
    # compute angle between rays[1:] and ray[0]
    for k in range(nbP-1):
        v = rays[k+1,:]
        dt = np.dot(u,v)
        cp = np.cross(u,v)
        at = get_angle(u,v)

        if np.signbit(np.dot(cp, n)):
            # the rotation angle is backward
            ang = np.degrees(2*np.pi - at)
        else:
            ang = np.degrees(at)
        angles[k+1] = ang

    # order vertices indices according to angles
    ordered = np.argsort(angles).tolist()
        # the sequence of position in 'angles' that sorts it in increasing angle order
    return ordered

# ==============================================================================
class N_Segment():
    """Basic features to deal with N-D segments."""

    # --------------------------------------------------------------------------
    def __init__(self, points):
        """Create the N-D segment.

        Args:
            points (list of 2 points): extremities of the segment
        """
        assert (points[0].ndim == 1) and (
            points[1].ndim == 1), "Invalid format for point"
        assert points[0].size == points[1].size, "Invalid format for point"
        self.pts = points
        self.dim = points[0].size
        # logger.debug("%dD segment [%s - %s] created", self.dim, *points)

    # --------------------------------------------------------------------------
    def __str__(self):
        """Return string representation of the segment."""
        return "[%s\n %s]" % (self.pts[0], self.pts[1])

    # --------------------------------------------------------------------------
    def getIntersectionWithRay(self, u):
        """Compute the intersection with the ray of direction u (from origin).

        Return the coordinates of the point where the ray intersects the segment.

        Args:
            u (np.array(N)): the(unit) direction vector of the ray

        Returns:
            float: distance from the origin of the ray to the intersection point
            with the segment (np.nan if no intersection)
        """
        # logger.debug("Compute %d-D segment intersection with ray %s",
        #              self.dim, u)
        d = np.nan

        p1 = self.pts[0]
        p2 = self.pts[1]

        v1 = - p1
        v2 = p2 - p1
        n = np.array([-u[1], u[0]])    # normal to u

        dv2n = np.dot(v2, n)
        if dv2n:
            t1 = np.cross(v2, v1) / dv2n
            t2 = np.dot(v1, n) / dv2n
            if 0. <= t2 <= 1.:
                d = t1
        else:
            if not p2[1]:
                d = p2[0]
        return d

    # --------------------------------------------------------------------------
    def getIntersectionWithPlane(self, nf):
        """Compute the intersection of the N-segment with a plane.

        Return the coordinates of the instersection point.

        Args:
            nf (np.array(N+1)): the normal form of hyperplane, [n0, d] with n0
                                the unit N-D vector normal to the hyperplane
                                and d the shortest distance to the origin.
        Returns:
            list: the coordinates of the intersection(s) point(s).
        """
        inter = []

        P1 = self.pts[0]
        P2 = self.pts[1]
        u = P2 - P1
        n0 = nf[:-1]
        d = nf[-1]/np.linalg.norm(n0)
        dot = np.dot(n0, u)

        if abs(dot) > EPS:
            # u and n0 not normal => find 1 point belonging to the hyperplane
            M = d*n0
            # test if the projection of this point on the line (p1,p2)
            # is inside the segment [p1 p2]
            w = P1 - M

            alpha = -np.dot(n0, w)
            s = alpha / dot

            # test if intersection is strictly inside the segment
            if 0 < s < 1:
                inter.append(P1 + s*u)

        # test if the extremities of the segment belong to the hyperplane
        if abs(np.dot(n0, P1) - d) < EPS:
            inter.append(P1)
        if abs(np.dot(n0, P2) - d) < EPS:
            inter.append(P2)

#        if len(inter):
#            logger.debug("Found intersection point(s): %s", inter)
        return inter


# ==============================================================================
class N_Triangle():
    """Basic features to deal with triangles in dimension N."""

    # --------------------------------------------------------------------------
    def __init__(self, tri):
        """Create the N-dimension triangle.

        Args:
            tri (np.array((3, N))): the linewise array of the N_D coordinates
                                    of the 3 triangle's vertices
        """
        self.tri = tri
        self.dim = tri.shape[1]

    # --------------------------------------------------------------------------
    def getIntersectionWithRay(self, u, o=None):
        """Get the distance to the origin of the ray if an intersection exists.

        Args:
            u (np.array(N)): the(unit) direction vector of the ray

        Returns:
            (float): the distance from the origin (0., 0., 0.) to the intersection 
                     point with the triangle using the Möller and Trumbore algorithm 
                     (np.nan if no intersection).
        """
        assert np.abs(np.linalg.norm(u)-1) <= EPS, \
               "Not a unit direction vector !"
        d = np.nan

        p1 = self.tri[0, :]
        p2 = self.tri[1, :]
        p3 = self.tri[2, :]
        if not o:
            o = np.zeros(self.dim)

        e1 = p2 - p1
        e2 = p3 - p1

        # compute determinant value
        pvec = np.cross(u, e2)
        det = np.dot(e1, pvec)

        if abs(det) >= EPS:
            invdet = 1./det
            tvec = o-p1
            cx = np.dot(tvec, pvec)*invdet

            qvec = np.cross(tvec, e1)
            cy = np.dot(u, qvec)*invdet

            if cx >= 0. and cy >= 0. and (cx+cy) <= 1.:
                d = np.dot(e2, qvec)*invdet
            # else:
            #     logger.debug("Ray outside the triangle")
        # else:
        #     logger.debug("Ray parallel to the triangle")
        return d


# ==============================================================================
class Zonotope():
    """Basic features to deal with zonotopes."""

    # --------------------------------------------------------------------------
    def __init__(self, Vmin, Vmax, M, coord_names=None):
        """Create the zonotope.

        Args:
            Vmin (np.array(nm): the vector of min values
            Vmax (np.array(nm): the vector of max values
            M (np.array(nc, nm)): the transformation matrix
            coord_name (list of strings): the names of the coordinates
                (optional input, if empty will be set to 'DOF_x', x in [1..nm])
        """
        # --- Initialisations ---
        assert Vmin.size == Vmax.size, \
                "Not the same number of actuators in upper and lower bounds"
        nm = Vmin.size
        assert M.shape[1] == nm, "Invalid number of actuators"
        nc = M.shape[0]
        if coord_names is None:
            coord_names = ["DOF %d" % (i+1) for i in range(nc)]
        assert nc == len(coord_names), "Invalid coordinate names"
        self.cNames = coord_names
        self.nc = nc
        self.nm = nm

        # --- Compute the H-representation of the Zonotope ---
        H, d = self.computeHrepr(Vmin, Vmax, M)

        # --- Compute the V-representation of the Zonotope (using the `polytope` python package)
        dd = np.array(d).reshape(len(d), 1)
        self.P = pt.Polytope(H.T, dd)
        logger.debug("Zonotope's H-repr:\n%s", self.P)
        ptVerts = pt.extreme(self.P) # caution: duplicate vertices are possible !
        rv = ptVerts.round(4).tolist()  # rounded vertice coordinates, possibly duplicated
        nodup = []
        for v in rv:    # remove duplicated vertices
            if v not in nodup:
                nodup.append(v)
        self.verts = np.array(nodup)
        logger.debug("Found %d vertices:\n%s",
                     len(self.verts), self.verts)

        # --- Compute the faces of the Zonotope (used for plotting)
        self.faces = self.getFaces()

        # --- Compute the faces of the Zonotope (used for plotting)
        self.edges = self.getEdges()

    # --------------------------------------------------------------------------
    def computeHrepr(self, Vmin, Vmax, M):
        """Compute the H-representation of the zonotope.
    
        Compute the H-representation of the zonotope through the hyperplane
        shifting method, based on the article by Gouttefarde et Krut (2010).

        This representation consists in a matrix H and a vector d which
        describe the hyperplanes forming the faces of the zonotope.

        The matrix of hyperplanes (dim 2p x nc) is built line-by-line.
        We use a loop on all possible combination of (n-1) column vectors of M 
        to find 'p' non colinear raw vectors noted c_ortho_{j in [0:p-1]}, and 
        the associated pairs {d_j_pos ; d_j_neg} of real values.

        A point X in IR^n (column vector of dimension 'n') belongs to the 
        zonotope if H.X < d, d column vector of dim '2p' with
        d[2k] = d_k_pos, and d[2k+1] = d_k_neg

        Args:
            Vmin (np.array(nm): the vector of min values
            Vmax (np.array(nm): the vector of max values
            M (np.array(nc, nm)): the transformation matrix
        """
        H = np.array([], float)
        d = []
        nc = self.nc
        nm = self.nm

        logger.info("Computing the H-representation from %d actuators and %d joints",
                    nm, nc)

        # Loop on each possible combination of (nc-1) column vectors of M
            
                                                      
        combin = list(itertools.combinations(range(nm), nc-1))
            # all combinations of (nc-1) indices in [0:nm[
        nbCombin = len(combin)
        logger.info("Running Hyperplane Shifting Method...")
        logger.debug("Loop on %d combinations", nbCombin)

        for j, subset in enumerate(combin):
        # for j in tqdm.tqdm(range(nbCombin),
                           # desc="Zonotope progression:"):
            # # logger.debug("combination %d, %s", j+1, combin[j])
            # subset = combin[j]  # i.e a list of (nc-1) values in range(nm)

            # reinit test (are there duplicated vectors in the set of column
            # vectors of indices in "combin" ?)
            c_dup = False

            WI0 = M[:, subset]
            U = np.linalg.svd(WI0)[0]   # because np.linalg.svd returns U, D and VT
                                          
            c_ortho = U[:, -1]  # the last column vector of U is orthogonal
                                # to the (nc-1) column vectors of W. Its norm is 1
                                # /!\ this is a row vector !

            # Check if the current 'c_ortho' is redundant with a previous one
            if j:        # j != 0 so H should not be empty...
                # check if the current c_ortho is colinear with any previous
                # one (stored in matrix H). Since there may be rounding error
                # in the SVD process, check if dot product are not ~ 1
                if np.max(np.dot(c_ortho.T, H)) > DEF_DOT_TOL:
                    # logger.debug("Iter %d, c_ortho is redundant with other "
                    #              "columns of H, H unchanged (%d column(s))",
                    #              j+1, H.shape[1])
                    c_dup = True

            # update constraints
            if not c_dup:
                # Compute offsets
                I = np.dot(c_ortho, M)

                # raw vector with 0 at indices where I<0
                vMax_pos = (I > 0)*Vmax
                # raw vector with 0 at indices where I>0
                vMax_neg = (I < 0)*Vmax
                vMin_pos = (I > 0)*Vmin
                vMin_neg = (I < 0)*Vmin

                # i_pos = np.where(I > 0)[0].tolist()
                # i_neg = np.where(I < 0)[0].tolist()

                dpos = np.dot(vMax_pos, I) + np.dot(vMin_neg, I)
                dneg = -  np.dot(vMax_neg, I) - np.dot(vMin_pos, I)

                # First loop ?
                if not j:
                    # Initialize H
                    H = c_ortho.reshape((nc, 1))
                else:
                    H = np.column_stack((H, c_ortho))

                # Append H
                H = np.column_stack((H, -c_ortho))
                d.append(dpos)
                d.append(dneg)
        return H,d

    # --------------------------------------------------------------------------
    def getEdges(self):
        """Return the list of edges (N-D segments) of the Zonotope.

        An edge is described as a pair of neighbour vertices indices (not 2D
        point coordinates).
        """
        logger.debug("Computing the edges of the Zonotope...")
        edges = []
        if self.nc == 2:
            edges = sp.ConvexHull(self.verts).simplices
        else:
            for f in self.faces:
                vList = f[0]    # the ordered list of vertices in the face
                nv = len(vList)
                vList.append(vList[0])  # to close the loop
                for k in range(nv):
                    s = set([vList[k], vList[k+1]]) # entered as a set to get rid off index order
                    if s not in edges:
                        edges.append(s)
            edges = [list(s) for s in edges]
        logger.debug("Found %d edges", len(edges))
        return edges

    # --------------------------------------------------------------------------
    def getFaces(self):
        """Return the list of all faces (hyperplans) of the Zonotope.
        
        Returns:
            (list) the list of the faces of the Zonotope. Each face is described
            as an ordered list of adjacent vertices indices and a the N-D vector
            normal to this face.
        """
        logger.debug("Computing the faces of the Zonotope...")

        if self.nc <= 2:
            logger.debug("No faces, Zonotope dimension is <3")
            logger.debug("Done!")
            return []

        logger.debug("Identifying all the faces of the Zonotope and their normals")
        pts = self.verts
        ch = sp.ConvexHull(pts) # used to computes simplices.
        faces = []

        normals = ch.equations

        #--- loop on simplices (triangles) to identify facets ---
        for k, simplex in enumerate(ch.simplices):
#            logger.debug("Processing simplex %d %s", k, simplex)
#            logger.debug("Normal %s", normals[k,:-1])
            facet = set(simplex)   # set of vertices that are in the same plan as this simplex

            # Test if this simplex is a subset of a previously identified facet
            new_facet = True    # initialize
            for f in faces:
                if facet.issubset(f[0]):    # only the vertices indices, not the normals
                    new_facet = False
#                    logger.debug("This face is already known")
                    break

            if new_facet:
                # select candidate vertices for this facet (all but those of the current simplex)
                verticesSubset = set(ch.vertices)
                for i in simplex:
                    verticesSubset.remove(i)

                # loop on vertices to identify those is the same plane as the current simplex
                for idx in verticesSubset:
#                    logger.debug("Candidate for face: vertex %d", idx)

                    # check if the vector formes by vertex 'idx' and the origin
                    #  of the simplex is colinear to the normal of the simplex
                    v = pts[idx] - pts[simplex[0]]
                    n = normals[k,:-1]
                    if abs(np.dot(n, v)) < EPS:
                        # Vertex is in the plan formed by the simplex
                        facet.add(idx)
#                   logger.debug("The plan formed by %s contains vertices %s", simplex, facet)
                
                # order the vertices of the facet and store it
                faces.append([ self.orderVerticesInFace(facet, n),
                               normals[k] ] )

        logger.debug("Found %d faces", len(faces))
        return faces

    # --------------------------------------------------------------------------
    def orderVerticesInFace(self, face, n):
        """Order vertices indices in a sequence of adjacent vertices.
        
        Args:
            face(set): the set of the vertices indices forming a face.
            n (np.array): the N-D vector normal to the facet

        Returns:
            (list) the ordered list of vertices indices
        """
#        logger.debug("Ordering vertices of face %s", face)
        indices = list(face)
        pts = self.verts

        ordered = orderPointsInPlane(pts[indices], n)

        return [indices[i] for i in ordered]

    # --------------------------------------------------------------------------
    def getNbCoord(self):
        """Return the number of coordinates defining the Zonotope."""
        return self.nc

    # --------------------------------------------------------------------------
    def getCoordId(self, coord):
        """Return the id of the coordinate in the Zonotope list of coords."""
        return self.cNames.index(coord)

    # --------------------------------------------------------------------------
    def plot2D(self, ax=None, clr='k', lbl=None):
        """Plot a 2D zonotope.

        Args:
            ax (plt.Axes): where to plot the Zonotope. If None, create the figure.
        """
        logger.debug("2D-plot of the Zonotope")

        # initialization
        if ax is None:
            f = plt.figure()
            f.suptitle("2D actuation zonotope")
            ax = f.add_subplot(111)
            ax.autoscale()
            ax.grid("on")
            ax.set_xlabel("Torque on '%s' (N.m)" % self.cNames[0])
            ax.set_ylabel("Torque on '%s' (N.m)" % self.cNames[1])

        # origin
        ax.plot(0, 0, 'ro')

        # vertices
        verts = self.verts
        for i,v in enumerate(verts):
            ax.text(*v, s=str(i))

        # edges
        for i1, i2 in self.edges:
            ax.plot([verts[i1][0], verts[i2][0]],
                    [verts[i1][1], verts[i2][1]],
                    color=clr,
                    ls='-',
                    marker='.')
        line = None
        if lbl:
            # redraw last line with a label for legend purpose
            line = ax.plot([verts[i1][0], verts[i2][0]],
                           [verts[i1][1], verts[i2][1]],
                           color=clr,
                           label=lbl)
        return [ax, line]

    # --------------------------------------------------------------------------
    def plot3D(self, ax=None, clr='b'):
        """Plot a 3D zonotope.

        Plot each face of the zonotope

        TODO: choose the coordinate to display if more than 3...

        Args:
            ax (plt.Axes): where to plot the Zonotope. If None, create the figure.

        Returns:
            handle to the figure's axis
        """
        logger.debug("3D-plot of the Zonotope")
        assert self.nc == 3, "Not a 3D zonotope"

        if ax is None:
            f = plt.figure()
            f.suptitle("3D actuation zonotope")
            ax = f.add_subplot(111, projection="3d")

            # initialisation of the plot
            ax.plot(0, 0, 0, 'ro')     # origin
            ax.set_xlabel("Torque on '%s' (N.m)" % self.cNames[0])
            ax.set_ylabel("Torque on '%s' (N.m)" % self.cNames[1])
            ax.set_zlabel("Torque on '%s' (N.m)" % self.cNames[2])
        else:
            assert isinstance(ax, plt.Axes), "Invalid axes"


        # Plot vertices and faces
        pts = self.verts

        # plot vertices
        for k,p in enumerate(pts):
            ax.text(*p, str(k))
            ax.plot(*p, 'ko')

        # plot faces
        faces = self.faces
        polys = []
        for f,n in faces:
            f.append(f[0])  # close the polygon by adding the starting point
            v = [pts[idx,:].tolist() for idx in f]
            polys.append(v)
        ax.add_collection3d(Poly3DCollection(polys, facecolor='0.8', edgecolor = 'k', alpha=0.65))

        return ax

    # --------------------------------------------------------------------------
    def getIntersectionPointWithRay(self, u, withPlot=False):
        r"""Return the intersection point with a ray, if any.

        The ray is defined by its unit vector 'u' (origin (0, 0, 0) ).

        As a zonotope is convex, there can be:
        - no intersection;
        - one intersection (the ray crosses a "corner vertex");
        - two intersections(opposites edges or faces). In this case, we only
        keep the intersection point with the same orientation as u (the dot
        product u-> . OP-> is positive).

        Args:
            u (np.array(nc)): a unit vector describing the direction of the ray
            withPlot (bool): plot the zonotope and the intersection point

        Returns:
            intersection point P_j \in IR**nc, None if no intersection.
        """
        logger.debug("Computing Zonotope intersection with ray %s", u)
        pt = None

        # 2D or >= 3D zonotope ?
        if self.nc == 2:
            simplex = N_Segment
        else:
            simplex = N_Triangle

        # loop on simplices
        vv = self.verts  # vertices of the zonotope
        ch = sp.ConvexHull(vv)
        
        for k, simp in enumerate(ch.simplices):  # loop on segments/triangles
            s = simplex(vv[simp, :])
            d = s.getIntersectionWithRay(u)
            if d >= 0:
                logger.debug("Found 1 intersection point with simplex #%d %s,"
                             " at distance d=%.2f from origin",
                             k, simp, d)
                pt = d*u

        if pt is None:
            logger.warning("No intersection found with ray %s", u)

        if withPlot:
            if self.nc == 2:
                ax = self.plot2D()[0]
                if pt is not None:
                    ax.plot([0., pt[0]], [0., pt[1]], 'ro:')
            if self.nc == 3:
                ax = self.plot3D()
                if pt is not None:
                    ax.plot([0., pt[0]], [0., pt[1]], [0., pt[2]], 'ro:')
        return pt

    # --------------------------------------------------------------------------
    def getIntersectionPointsWithPlane(self, nf, withPlot=False):
        r"""Return the list of intersection points with a plane.

        Args:
            nf (np.array(N+1)): the normal form of hyperplane, [n0, d] with n0
                                the unit N-D vector normal to the hyperplane
                                and d the shortest distance to the origin.
            withPlot (bool): plot the zonotope and the intersection points

        Returns:
            inter (list of N-D points) intersection(s) point(s),
                   ordered as sequence of consecutive points.
        """
        logger.debug("Looking for intersection with plan %s", nf)
        pts = self.verts
        edges = self.edges
        inter = []
        
        # loop on edges
        for idV1, idV2 in edges:
            seg = N_Segment([pts[idV1, :], pts[idV2, :]])
            lp = seg.getIntersectionWithPlane(nf)
            # update inter but avoid duplicated points
            for p in lp:
                # check if p was already identified
                # It's not easy to test with array objects
                # => use tuple and lists, and round to compensate for
                # computation approximations
                p_ = p.round(3).tolist()
                if not inter.count(p_):
                    inter.append(p_)
        logger.debug("Found %d intersection point(s)!", len(inter))

        # Reorder by adjacencies for plotting
        orderedInt = orderPointsInPlane(np.array(inter), nf[:-1])

        # Plots
        if withPlot:
            ax = self.plot3D()

            # plot intersection points
            for k, p in enumerate(inter):
                ax.plot(*p, 'r.')
                ax.text(*p, "I%0d"%k, c='r')

            lines = []
            loop = orderedInt + [orderedInt[0],]
            for i1, i2 in zip(loop[0:-1], loop[1:]):
                lines.append([inter[i1], inter[i2]])
            ax.add_collection3d(Line3DCollection(lines,color='r', ls='--'))

        return inter

    # --------------------------------------------------------------------------
    def __str__(self):
        """Return string representation of the Zonotope."""
        return str(self.P)


# =====================================
# MAIN CODE
# =====================================
def main():
    """Test the class ZNT_Processor."""
    # initialize the zonotope's characteristics
    # ------------------------------------------
    logger.info("Testing the class ZNT_Processor...")


    # Pseudo "biomechanical" 2D Zonotope:
    #     6 muscles, 3 joints, coupled joints and muscles
    # 2D-test...
    #-----------------------------------------
    M = DEF_M2D
    nc, nm = M.shape
    FMax = np.ones(nm, float)
    FMin = 0.15*FMax
    my2dZnt = Zonotope(FMin, FMax, M)
#    my2dZnt.plot2D()

    # Pseudo "biomechanical" 3D Zonotope:
    #       6 muscles, 3 joints, coupled joints and muscles
    # 3D-test...
    #-----------------------------------------
    M = DEF_M3D
#    M = DEF_CUB
    nc, nm = M.shape
    FMax = np.ones(nm, float)
    FMin = 0.15*FMax
    my3dZnt = Zonotope(FMin, FMax, M)
#    my3dZnt.plot3D()


    # Search and plot intersections with a ray
    #   (2D or 3D test) and a plane (3D test)
    # ----------------------------------------

    # 3D-zonotope / ray and plane intersection
    theta = DEF_THETA
    phi = DEF_PHI
    u = np.array([math.sin(theta) * math.cos(phi),
                  math.sin(theta) * math.sin(phi),
                  math.cos(theta)])

    hp = np.array([0., 0., 1., 0.25])
    my3dZnt.getIntersectionPointsWithPlane(hp, withPlot=True)
    my3dZnt.getIntersectionPointWithRay(u, withPlot=True)

    # 2D-zonotope / ray and plane intersection
    theta = DEF_THETA
    u = np.array([math.cos(theta), math.sin(theta)])
    my2dZnt.getIntersectionPointWithRay(u, withPlot=True)

    plt.show()


# ==============================================================================
if __name__ == "__main__":
    # logging stuff
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        '%(name)s - %(levelname)s - %(filename)s:%(lineno)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    main()
