# -*- coding: utf-8 -*-
"""
Created on Thu Sep 01 16:03:53 2011

@author: sturgeon
"""
import numpy

def find_exterior_faces_with_consistent_cs(e,num_n):
    # Find exterior faces with a consistent CS in the element
    num_e = numpy.size(e,0)
    
    # Find the exterior faces
    sf = numpy.array([[0,1,2,3],[4,5,6,7],[0,1,5,4],[3,2,6,7],[0,3,7,4],[1,2,6,5]])

    fM = {}
    for k in range(num_e):
        for i in range(6):
            face = [e[k][sf[i][0]],e[k][sf[i][1]],e[k][sf[i][2]],e[k][sf[i][3]]]
            face.sort()
            fInd = tuple(face)
            if fInd in fM:
                fM[fInd] = numpy.concatenate((fM[fInd],numpy.array([[k,i]])))
            else:
                fM[fInd] = numpy.array([[k,i]])
    
    faceKeys = list(fM.keys())
    faceInd = list(fM.values())
    faceIsExtFace = numpy.zeros((len(faceKeys),1),'bool')
    nIsOnExt = [0]*num_n
    fExt = []
    fExtElNum = [] # element that the exterior face belongs to
    fExtFaceNum = [] # face number from sf
    for i in range(len(faceKeys)):
        if numpy.size(faceInd[i],0)==1:
            faceIsExtFace[i] = 1
            nIsOnExt[faceKeys[i][0]] = 1
            nIsOnExt[faceKeys[i][1]] = 1
            nIsOnExt[faceKeys[i][2]] = 1
            nIsOnExt[faceKeys[i][3]] = 1
            fExt.append(e[faceInd[i][0][0],sf[faceInd[i][0][1]]])
            fExtElNum.append(faceInd[i][0][0])
            fExtFaceNum.append(faceInd[i][0][1])
    
    
    fExt = numpy.array(fExt)
    fExtFaceNum = numpy.array(fExtFaceNum)
    fExtElNum = numpy.array(fExtElNum)
    nIsOnExt = numpy.array(nIsOnExt,'bool')
    
    return [fExt,nIsOnExt,fExtFaceNum,fExtElNum]

def find_exterior_faces(e,num_n):
    # Find exterior faces with consistent normals
    num_e = numpy.size(e,0)
    
    # Find the exterior faces
    sf = numpy.array([[0,1,2,3],[5,4,7,6],[4,5,1,0],[3,2,6,7],[4,0,3,7],[1,5,6,2]])
    fM = {}
    for k in range(num_e):
        for i in range(6):
            face = [e[k][sf[i][0]],e[k][sf[i][1]],e[k][sf[i][2]],e[k][sf[i][3]]]
            face.sort()
            fInd = tuple(face)
            if fInd in fM:
                fM[fInd] = numpy.concatenate((fM[fInd],numpy.array([[k,i]])))
            else:
                fM[fInd] = numpy.array([[k,i]])
    
    faceKeys = list(fM.keys())
    faceInd = list(fM.values())
    faceIsExtFace = numpy.zeros((len(faceKeys),1),'bool')
    nIsOnExt = [0]*num_n
    fExt = []
    fExtElNum = [] # element that the exterior face belongs to
    fExtFaceNum = [] # face number from sf
    for i in range(len(faceKeys)):
        if numpy.size(faceInd[i],0)==1:
            faceIsExtFace[i] = 1
            nIsOnExt[faceKeys[i][0]] = 1
            nIsOnExt[faceKeys[i][1]] = 1
            nIsOnExt[faceKeys[i][2]] = 1
            nIsOnExt[faceKeys[i][3]] = 1
            fExt.append(e[faceInd[i][0][0],sf[faceInd[i][0][1]]])
            fExtElNum.append(faceInd[i][0][0])
            fExtFaceNum.append(faceInd[i][0][1])
    
    
    fExt = numpy.array(fExt)
    fExtFaceNum = numpy.array(fExtFaceNum)
    fExtElNum = numpy.array(fExtElNum)
    nIsOnExt = numpy.array(nIsOnExt,'bool')
    
    return [fExt,nIsOnExt,fExtFaceNum,fExtElNum]

def find_faces_that_are_edge_neighbors_to_ext_faces(e,fExt):

    # Store a dictionary of all the edges used by the exterior faces (fExt) and 
    # the fExt number that uses the edge
    se = numpy.array([[0,1],[1,2],[2,3],[3,0]])
    fExtUseEdge = {}
    for k in range(fExt.shape[0]):
        for i in range(4):
            eInd = (min(fExt[k][se[i,0]],fExt[k][se[i,1]]),max(fExt[k][se[i,0]],fExt[k][se[i,1]]))
            if eInd in fExtUseEdge:
                fExtUseEdge[eInd].append(k)
            else:
                fExtUseEdge[eInd] = [k]
    
    f_Nei = numpy.ones(fExt.shape,'int')
    f_Nei = -f_Nei
    for k in range(fExt.shape[0]):
        for i in range(4):
            eInd = (min(fExt[k][se[i,0]],fExt[k][se[i,1]]),max(fExt[k][se[i,0]],fExt[k][se[i,1]]))
            if eInd in fExtUseEdge:
                neighbors = fExtUseEdge[eInd]
                if neighbors[0] == k:
                    f_Nei[k,i] = neighbors[1]
                else:
                    f_Nei[k,i] = neighbors[0]
    return f_Nei

def find_boundary_crvs(e,num_n):
    # Find the boundary curves
    eM = {}
    se = numpy.array([[0,1],[1,2],[2,3],[3,0],[4,5],[5,6],[6,7],[7,4],[0,4],[1,5],[2,6],[3,7]])
    for k in range(len(e)):
        for i in range(12):
            eInd = (min(e[k][se[i,0]],e[k][se[i,1]]),max(e[k][se[i,0]],e[k][se[i,1]]))
            if eInd in eM:
                eM[eInd].append(k)
            else:
                eM[eInd] = [k]
    
    # Identify boundary edges
    edges = list(eM.keys())
    edgeCount = list(eM.values())
    usedOnce = [0]*len(edges)
    for i in range(len(edges)):
        if numpy.size(edgeCount[i])==1:
            usedOnce[i] = 1
    
    usedOnce = numpy.array(usedOnce)
    edges = numpy.array(edges)
    bnd_edges = edges[numpy.nonzero(usedOnce),:][0]
    
    # Mark the boundary nodes and store their neighboring boundary nodes
    eB = numpy.ones((num_n,2),int)
    eB = -eB
    for k in range(bnd_edges.shape[0]):
        if eB[bnd_edges[k,0],0]==-1:
            eB[bnd_edges[k,0],0] = bnd_edges[k,1]
        else:
            eB[bnd_edges[k,0],1] = bnd_edges[k,1]
        
        if eB[bnd_edges[k,1],0]==-1:
            eB[bnd_edges[k,1],0] = bnd_edges[k,0]
        else:
            eB[bnd_edges[k,1],1] = bnd_edges[k,0]
    
    return(bnd_edges,eB) #eB stores the node before and after a node on a bnd_edge

