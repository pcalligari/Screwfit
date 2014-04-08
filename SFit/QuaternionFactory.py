## Automatically adapted for numpy.oldnumeric Feb 13, 2009 by 

import numpy.oldnumeric as N
from MMTK import *
from Scientific.Geometry import Quaternion, Vector, Transformation
from Scientific.IO.ArrayIO import *
from sys import argv
from copy import deepcopy
from getopt import *
import string
import os



def findQuaternionMatrix(collection, point_ref, conf1, conf2 = None, matrix = True ):
	
        universe = collection.universe()
        if conf1.universe != universe:
                raise ValueError, "conformation is for a different universe"
        if conf2 is None:
                conf1, conf2 = conf2, conf1
        else:
                if conf2.universe != universe:
                        raise ValueError, "conformation is for a different universe"
        ref = conf1
        conf = conf2
        weights = universe.masses()
        weights = weights/collection.mass()
        ref_cms = point_ref.position().array
        pos = N.zeros((3,), N.Float)
        pos = point_ref.position(conf).array
        possq = 0.
        cross = N.zeros((3, 3), N.Float)
        for a in collection.atomList():
                r = a.position(conf).array - pos
                r_ref = a.position(ref).array-ref_cms
                w = weights[a]
                possq = possq + w*N.add.reduce(r*r) \
                                                + w*N.add.reduce(r_ref*r_ref)
                cross = cross + w*r[:, N.NewAxis]*r_ref[N.NewAxis, :]
        k = N.zeros((4, 4), N.Float)
        k[0, 0] = -cross[0, 0]-cross[1, 1]-cross[2, 2]
        k[0, 1] = cross[1, 2]-cross[2, 1]
        k[0, 2] = cross[2, 0]-cross[0, 2]
        k[0, 3] = cross[0, 1]-cross[1, 0]
        k[1, 1] = -cross[0, 0]+cross[1, 1]+cross[2, 2]
        k[1, 2] = -cross[0, 1]-cross[1, 0]
        k[1, 3] = -cross[0, 2]-cross[2, 0]
        k[2, 2] = cross[0, 0]-cross[1, 1]+cross[2, 2]
        k[2, 3] = -cross[1, 2]-cross[2, 1]
        k[3, 3] = cross[0, 0]+cross[1, 1]-cross[2, 2]
        for i in range(1, 4):
                for j in range(i):
                        k[i, j] = k[j, i]
        k = 2.*k
        for i in range(4):
                k[i, i] = k[i, i] + possq - N.add.reduce(pos*pos)
        import numpy.oldnumeric.linear_algebra as LinearAlgebra
        e, v = LinearAlgebra.eigenvectors(k)
        i = N.argmin(e)
        v = v[i]
        if v[0] < 0: v = -v
        if e[i] <= 0.:
                rms = 0.
        else:
                rms = N.sqrt(e[i])
	if matrix:
		emax = N.argmax(e)
		QuatMatrix = v
		return Quaternion.Quaternion(QuatMatrix),v, e, e[i],e[emax], rms
	else:
		return Quaternion.Quaternion(v), Vector(ref_cms), Vector(pos), rms


def findGenericTransformationAsQuaternion(collection, point_ref, conf1, conf2 = None):
        return findQuaternionMatrix(collection, point_ref, conf1, conf2 = None, matrix = False)

def findGenericTransformation(collection,point_ref,conf1, conf2 = None):
        q, cm1, cm2, rms = findGenericTransformationAsQuaternion(collection,point_ref,conf1,conf2)
        return Transformation.Translation(cm2),q.asRotation(),Transformation.Translation(-cm1),rms
