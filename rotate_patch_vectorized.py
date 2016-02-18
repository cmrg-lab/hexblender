# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 08:34:38 2011

@author: sturgeon
"""

import numpy

def rotate_patch_new(face, crossedEdge, uwdir, Bx, By, Bz,nei):

    # Rotate so that the paramater along the crossedEdge matches directions 
    # across the patches.  ie when the crossedEdge is 2 3 then the new patch 
    # should be rotated so edge 1 4 matches

    Bx = Bx.transpose([1,2,0])
    By = By.transpose([1,2,0])
    Bz = Bz.transpose([1,2,0])
    
    if (uwdir[0] == 1 and uwdir[1] == 2):
        fn0 = 0
        fn1 = 3
    elif (uwdir[0] == 0 and uwdir[1] == 3):
        fn0 = 1
        fn1 = 2
    elif (uwdir[0] == 3 and uwdir[1] == 2):
        fn0 = 0
        fn1 = 1
    elif (uwdir[0] == 0 and uwdir[1] == 1):
        fn0 = 3
        fn1 = 2
    
    rotate = 1
    needsRotating = - numpy.logical_and(face[:,fn0] == crossedEdge[:,0] ,face[:,fn1] == crossedEdge[:,1])
    
    while (any(needsRotating) and  rotate <=10):
        face[needsRotating] = numpy.vstack([face[needsRotating,1], face[needsRotating,2], face[needsRotating,3], face[needsRotating,0]]).transpose()
        nei[needsRotating] = numpy.vstack([nei[needsRotating,1], nei[needsRotating,2], nei[needsRotating,3], nei[needsRotating,0]]).transpose()
        
        rotate = rotate +1
        Bx[:,:,needsRotating] = rotate_parameterization(Bx[:,:,needsRotating])
        By[:,:,needsRotating] = rotate_parameterization(By[:,:,needsRotating])
        Bz[:,:,needsRotating] = rotate_parameterization(Bz[:,:,needsRotating])
        
        if rotate==5:
            face[needsRotating] = numpy.vstack([face[needsRotating,3], face[needsRotating,2], face[needsRotating,1], face[needsRotating,0]]).transpose()
            nei[needsRotating] = numpy.vstack([nei[needsRotating,2], nei[needsRotating,1], nei[needsRotating,0], nei[needsRotating,3]]).transpose()
            Bx[:,:,needsRotating] = flip_paramaterization(Bx[:,:,needsRotating])
            By[:,:,needsRotating] = flip_paramaterization(By[:,:,needsRotating])
            Bz[:,:,needsRotating] = flip_paramaterization(Bz[:,:,needsRotating])
        
        needsRotating = - numpy.logical_and(face[:,fn0] == crossedEdge[:,0] ,face[:,fn1] == crossedEdge[:,1])
        if rotate == 11:
            print("ERROR: Rotated 360 degrees")
    
    Bx = Bx.transpose([2,0,1])
    By = By.transpose([2,0,1])
    Bz = Bz.transpose([2,0,1])
    
    return face, Bx, By, Bz, nei

def rotate_parameterization(B):
    
    Bul = B[0:2,0:2,:]
    Bur = B[0:2,2:4,:]
    Bll = B[2:4,0:2,:]
    Blr = B[2:4,2:4,:]
    
    Bn = numpy.zeros((B.shape), 'float32')
    
    Bn[0:2,0:2,:] = numpy.rot90(Bul,3).copy()
    Bn[0:2,2:4,:] = -numpy.rot90(Bll,3).copy()
    Bn[2:4,0:2,:] = numpy.rot90(Bur,3).copy()
    Bn[2:4,2:4,:] = -numpy.rot90(Blr,3).copy()

    return Bn

def flip_paramaterization(B):
    Bn = numpy.zeros((B.shape), 'float32')
    for k in range(B.shape[2]):
        Bn[0:2,0:2,k] = numpy.fliplr(B[0:2,0:2,k])
        Bn[0:2,2:4,k] = -numpy.fliplr(B[0:2,2:4,k])
        Bn[2:4,0:2,k] = numpy.fliplr(B[2:4,0:2,k])
        Bn[2:4,2:4,k] = -numpy.fliplr(B[2:4,2:4,k])
    
    return Bn

