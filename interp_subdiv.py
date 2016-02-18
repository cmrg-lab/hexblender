# -*- coding: utf-8 -*-
"""
Perform Interpolatory subdivision surface on a quad mesh
Based on G.Li, W. Ma and H. Bao, A New Interpolatory subdivision for 
Quadrilateral Meshes, Computer Graphics forum, Vol 24 (2005) 1 pp3-16.
Similar to Kobbelt interpolatory subdivision scheme.

e = num_e x 8 array of the elements
n = num_n x 3 array of the vertex coordinates
It is assumed that the exterior faces are on e[:,0:4] and e[:,4:8] (this allows 
for exterior faces not on the endo or epi to be subdivided linearly).

Called by hexInterpSubDiv(e,n)

Returns:    fPn - dictionary storing the new face points
                the key is the sorted indices of the face
                fPn[face[k].sort()] = [x, y, z]
            ePn - dictionary storing the new edge points
                the key is the sorted indices of the edge
                ePn[(min(edge[k][0],edge[k][1]),max(edge[k][0],edge[k][1]))] = [x, y, z]

Revisions:
1.1
    Update to find the fEndoEpi and fExtOther.  The subdivision is run on both sets of
    faces and the new points are collected and returned together.
1.2
    Eliminate the need to find the fEndoEpi and fExtOther, which was based on the Endo
    and Epi being on eta3 surfaces.  Surfaces that are seperated by a 'boundary edge' are
    treated as seperate surfaces.  A 'boundary edge' is identified where an edge between
    exterior faces belongs to the same element.
1.3
    Identify areas on the exterior faces where an edge is not used twice (extraordinary 
    exterior edges) and computes a C0 subdivision across these edges.  The 'boundary edge'
    now includes exterior edges that are used by > 2 elements and < 2 elements.
    
            
@author: Greg Sturgeon
gregorymsturgeon@hotmail.com
October 31, 2011
"""
import numpy
from hexblender.interp_subdiv_functions import include_duplicate_point
from hexblender.interp_subdiv_functions import rotate_face_align_right
from hexblender.interp_subdiv_functions import vertex_weights_betas_gammas
from hexblender.interp_subdiv_functions import circ_shift

def interp_subdiv(e,n):
    if e.min()==1:
        e = e-1
    
    num_n = n.shape[0]
    num_e = e.shape[0]
    
    # Weighting values
    alpha0 = -1./8.
    alpha1 = 6./8.
    alpha2 = 3./8.
    
    fPn = {}
    ePn = {}
    
    # Find the exterior faces
    sf = numpy.array([[0,1,2,3],[4,5,6,7],[0,4,5,1],[3,2,6,7],[0,3,7,4],[1,5,6,2]])
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
    
    nIsOnExt = [0]*num_n
    fExt = []
    for i in range(len(faceKeys)):
        if numpy.size(faceInd[i],0)==1:
            fExt.append(e[faceInd[i][0][0],sf[faceInd[i][0][1]]])
            nIsOnExt[faceKeys[i][0]] = 1
            nIsOnExt[faceKeys[i][1]] = 1
            nIsOnExt[faceKeys[i][2]] = 1
            nIsOnExt[faceKeys[i][3]] = 1
    
    
    fExt = numpy.array(fExt)
    f = fExt
    num_f = f.shape[0]
    
    # Find the number of connected exterior faces for each EXTERIOR edge
    sef = numpy.array([[0,1],[1,2],[2,3],[3,0]])
    eeM = {} # Exterior edge face Matrix
    for k in range(numpy.size(fExt,0)):
        for i in range(4):
            eInd = (min(fExt[k,sef[i,0]],fExt[k,sef[i,1]]),max(fExt[k,sef[i,0]],fExt[k,sef[i,1]]))
            if eInd in eeM:
                eeM[eInd].append(k)
            else:
                eeM[eInd] = [k]
    
    
    # Find the boundary curves
    # Find the number of elements for each EXTERIOR edge
    eM = {} # Edge Matrix
    
    se = numpy.array([[0,1],[1,2],[2,3],[3,0],[4,5],[5,6],[6,7],[7,4],[0,4],[1,5],[2,6],[3,7]])
    for k in range(len(e)):
        for i in range(12):
            eInd = (min(e[k][se[i,0]],e[k][se[i,1]]),max(e[k][se[i,0]],e[k][se[i,1]]))
            if eInd in eeM:
                if eInd in eM:
                    eM[eInd].append(k)
                else:
                    eM[eInd] = [k]
    
    # Identify boundary edges or extraordinary ext edges
    edges = list(eM.keys())
    edgeCount = list(eM.values())
        
    notUsedTwice = [1]*len(edges)
    for i in range(len(edges)):
        if numpy.size(edgeCount[i])==2:
            notUsedTwice[i] = 0
    
    notUsedTwice = numpy.array(notUsedTwice)
    edges = numpy.array(edges)
    bnd_edges = edges[numpy.nonzero(notUsedTwice)[0],:]
    
    
    nNotOnBndCrv = numpy.ones((num_n,1),int)
    nNotOnBndCrv[bnd_edges[:,0]] = 0
    nNotOnBndCrv[bnd_edges[:,1]] = 0    
    
    eM = {}
    se = numpy.array([[0,1],[1,2],[2,3],[3,0]])
    eConnV = {}
    connFs = {}
    for k in range(num_f):
        for i in range(4):
            eInd = (min(f[k][se[i,0]],f[k][se[i,1]]),max(f[k][se[i,0]],f[k][se[i,1]]))
            if eInd in eM:
                eM[eInd].append(k)
            else:
                eM[eInd] = [k]
            if f[k,se[i,0]] in eConnV:
                eConnV[f[k,se[i,0]]].append(f[k,se[i,1]])
            else:
                eConnV[f[k,se[i,0]]] = [f[k,se[i,1]]]
            if f[k,se[i,1]] in eConnV:
                eConnV[f[k,se[i,1]]].append(f[k,se[i,0]])
            else:
                eConnV[f[k,se[i,1]]] = [f[k,se[i,0]]]
            if f[k,i] in connFs:
                connFs[f[k,i]].append(k)
            else:
                connFs[f[k,i]] = [k]
    
    # Store only the unique values in eConnV
    eConnV_keys = list(eConnV.keys())
    eConnV_values = list(eConnV.values())
    for k in range(len(eConnV)):
        eConnV_k = numpy.array(eConnV_values[k])
        [eConnV_k_val,eConnV_k_ind] = numpy.unique(eConnV_k,return_index=True)
        eConnV_k_ind.sort()
        eConnV[eConnV_keys[k]] = eConnV_k[eConnV_k_ind]
    
    # Store only the unique values in connFs
    connFs_keys = list(connFs.keys())
    connFs_values = list(connFs.values())
    for k in range(len(connFs)):
        connFs_k = numpy.array(connFs_values[k])
        [connFs_k_val,connFs_k_ind] = numpy.unique(connFs_k,return_index=True)
        connFs_k_ind.sort()
        connFs[connFs_keys[k]] = connFs_k[connFs_k_ind]
    
    
    
    qF = {}
    qE = {}
    # Perform the subdivision
    for k in range(num_n):
        if all(numpy.logical_and(nNotOnBndCrv[k],nIsOnExt[k])):
            fNext = connFs[k][0]
            fConn = f[connFs[k][0],:]
            # find if any point should be included in the 1 neighborhood twice
            pDup = include_duplicate_point(connFs[k],f,k)
            num_conn_edges = numpy.size(eConnV[k])
            # p stores the vertex indices in the 1 neighborhood stencil (see Figure 4)
            p = numpy.ones((2*num_conn_edges,1),int)
            p = -p
            ip = 0; # index into p
            p[0] = eConnV[k][0]
            fCCW = [] # Store the faces touching vertex k in a CCW direction
            # Create the 1 neighborhood stencil
            while numpy.any(p==-1):
                if not all(numpy.in1d(fNext,fCCW)):
                    fCCW.append(fNext)
                # Rotate fConn so that edge k p(ip) is the right edge
                fConn = rotate_face_align_right(fConn,[k,p[ip]])
                # Store the vertex indices in the 1 neighborhood stencil
                if not all(numpy.in1d(fConn[1],p)):
                    ip+=1
                    p[ip] = fConn[1]
                if numpy.logical_and(fConn[1]==pDup, not p[ip]==pDup):
                    ip+=1
                    p[ip] = fConn[1]
                if not all(numpy.in1d(fConn[0],p)):
                    ip+=1
                    p[ip] = fConn[0]
                if numpy.logical_and(fConn[0]==pDup, not p[ip]==pDup):
                    ip+=1
                    p[ip] = fConn[0]
                if numpy.logical_or(not all(numpy.in1d(fConn[3],p)),fConn[3]==pDup):
                    ip+=1
                    p[ip] = fConn[3]
                if numpy.logical_and(fConn[3]==pDup, not p[ip]==pDup):
                    ip+=1
                    p[ip] = fConn[3]
                # Find the next face that touches the current edge k p(ip)
                # which hasn't been stored in fCCW already
                # May need to check that eM has two values            
                indCur = (min(k,p[ip,0]),max(k,p[ip,0]))
                if indCur in eM:
                    if not all(numpy.in1d(numpy.array(eM[indCur][0]),numpy.array(fCCW))):
                        fNext = eM[indCur][0]
                        fConn = f[fNext,:]
                    elif len(eM[indCur])>1:
                        if not all(numpy.in1d(numpy.array(eM[indCur][1]),numpy.array(fCCW))):
                            fNext = eM[indCur][1]
                            fConn = f[fNext,:]
                        
                    elif any(p==-1):
                        print('ERROR: Check to verify that the mesh is manifold')
                    
                elif any(p==-1):
                    print('ERROR: Check to verify that the mesh is manifold')
            
            val=len(p)/2
            # Compute the Betas and Gammas
            [B,G] = vertex_weights_betas_gammas(val,alpha0,alpha1,alpha2)
            fCCW = numpy.array(fCCW) 
            # Rotate the order of p and fCCW to fill in all the qF's that touch the vertex k
            for j in range(len(fCCW)):
                pCirc = circ_shift(p.copy(),-2*j)
                fCCWcirc = circ_shift(fCCW.copy(),-j)
                # pIncCent includes the center vertex k in the 1 neighborhood
                c = numpy.ones((1,1),int)
                c[0] = k
                pIncCent = numpy.concatenate((c,pCirc),axis=0)
                if fCCWcirc[0] in qF:
                    qF[fCCWcirc[0]].append((B*n[pIncCent[:,0],:]).sum(0))
                else:
                    qF[fCCWcirc[0]] = [(B*n[pIncCent[:,0],:]).sum(0)]
                
                qEind = (min(k,pCirc[0,0]),max(k,pCirc[0,0]))
                if qEind in qE:
                    qE[qEind].append((G*n[pIncCent[:,0],:]).sum(0))
                else:
                    qE[qEind] = [(G*n[pIncCent[:,0],:]).sum(0)]
        
    qE_keys = list(qE.keys())
    qE_values = list(qE.values())
    
    qF_keys = list(qF.keys())
    qF_values = list(qF.values())
    
    for k in range(len(qF_values)):
        L = len(qF_values[k])
        face = f[qF_keys[k]].copy()
        face.sort()
        fInd = tuple(face)
        fPn[fInd] = (numpy.array(qF_values[k])/L).sum(0)
    
    
    for k in range(len(qE_values)):
        L = len(qE_values[k])
        ePn[qE_keys[k]] = (numpy.array(qE_values[k])/L).sum(0)
    
    
    # Work on the position of the vertices on extraordinary ext edges
    # Store the vertices on the boundary edges that touch the kth vertex
    eB = numpy.ones((num_n,3),int)
    eB = -eB
    for k in range(numpy.size(bnd_edges ,0)):
        if eB[bnd_edges[k,0],0]==-1:
            eB[bnd_edges[k,0],0] = bnd_edges[k,1]
        elif eB[bnd_edges[k,0],1]==-1:
            eB[bnd_edges[k,0],1] = bnd_edges[k,1]
        else:
            eB[bnd_edges[k,0],2] = 1
        
        if eB[bnd_edges[k,1],0]==-1:
            eB[bnd_edges[k,1],0] = bnd_edges[k,0]
        elif eB[bnd_edges[k,1],1]==-1:
            eB[bnd_edges[k,1],1] = bnd_edges[k,0]
        else:
            eB[bnd_edges[k,1],2] = 1
    
    
    isCorner = numpy.zeros((num_n,1),bool)
    for k in range(num_n):
        if (eB[k,2]==1):
            isCorner[k] = 1
    
    eB[numpy.nonzero(isCorner==1)[0],:] = -1
    
    # Weighting and averaging of the boundary vertices
    for k in range(numpy.size(bnd_edges,0)):
        p0 = bnd_edges[k,0]
        p1 = bnd_edges[k,1]
        if numpy.logical_or(eB[p0][0] ==-1, eB[p1][0] ==-1):
            eInd = (min(bnd_edges[k,0],bnd_edges[k,1]),max(bnd_edges[k,0],bnd_edges[k,1]))
            ePn[eInd] = .5 * (n[p0,:]+n[p1,:])
        else:
            if eB[p0][0] == p1:
                pminus1 = eB[p0][1]
            else:
                pminus1 = eB[p0][0]
            
            if eB[p1][0] == p0:
                p2 = eB[p1][1]
            else:
                p2 = eB[p1][0]
            
            qbE0 = alpha0*n[pminus1,:] + alpha1*n[p0,:] + alpha2*n[p1,:]
            qbE1 = alpha2*n[p0,:] + alpha1*n[p1,:] + alpha0*n[p2,:]
            eInd = (min(bnd_edges[k,0],bnd_edges[k,1]),max(bnd_edges[k,0],bnd_edges[k,1]))
            ePn[eInd] = .5 * (qbE0+qbE1)
    
    return [fPn,ePn]
