# -*- coding: utf-8 -*-
"""
regularizeElementsVect(e,n,eNN,nNN,itt)

e - elements of base mesh (num_e,8)
n - nodes of base mesh (num_n,3)
eNN - elements of twice subdivided mesh
nNN - nodes of twice subdivided mesh
itt - number of itterations to perform


Regularization of the elements where the motion is constrained to the surface 
through a pull-back and push-forward.

1) Hermite patches are fit to the exterior surfaces of (eNN,nNN) based on the 
structure of the base mesh. 
2) extraordinary exterior edges are found and marked so that a node on one of 
these edges will move along that edge. 
3) obatin the motion vectors - The HexDual (Vartziotis & Wipper) is used to 
determine the motion vectors.  
4) the motion is pulled-back to the parameter plane, where the node's u,v 
location on a patch is updated along with information about the patch the node 
lies on.  This motion is performed by the function MoveInSurface.  
5) Some motion is allowed in the normal direction with in a tolerance zone.  
6) invalid element handling where the nodal positions are not updated if they 
create poorly shaped elements.


Note: the trade-off between improving element skew and aspect ratio could be improved upon.

July 13, 2012
Greg Sturgeon
"""

import numpy
from hexblender.regularize_elements_pullback_functions import initialize_hermite_patch
from hexblender.hex_dual_non_iso_scale import hex_dual_non_iso
from hexblender.hex_dual import invalid_element_handling
from hexblender.move_in_surface import move_in_surface
from hexblender.misc_functions import find_faces_that_are_edge_neighbors_to_ext_faces
from hexblender.misc_functions import find_exterior_faces_with_consistent_cs
from hexblender.regularize_elements import extraordinary_ext_edges
from hexblender.regularize_elements import vert_normal
from hexblender.regularize_elements import normalize_v3

def regularize_elements_vect(e, n, eNN, nNN, itt, immobilizeRidges = True):

    # residual to the surface
    r = 1
    
    if e.min()==1:
        e = e-1
    
    if eNN.min()==1:
        eNN = eNN-1
    
    if n.dtype.name == 'float64':
        n = numpy.float32(n)
    
    num_e = numpy.size(e,0)
    num_n = numpy.size(n,0)
    
    [ptch_face,nIsOnExt,ptchFaceNum,ptchElNum] = find_exterior_faces_with_consistent_cs(e,num_n)
    
    # For the exterior faces fExt find the nodal derivatives from the subdivision
    [ptch_Bx,ptch_By,ptch_Bz,p_u,p_w,p_Bx,p_By,p_Bz,p_ptchNum,p_face] = initialize_hermite_patch(ptch_face,eNN,nNN,ptchFaceNum,ptchElNum)
    
    ptch_Nei = find_faces_that_are_edge_neighbors_to_ext_faces(e,ptch_face)
    
    ################################################################################
    # Now work on the subdivided mesh
    
    e = eNN
    n = nNN
    num_e = numpy.size(e,0)
    num_n = numpy.size(n,0)
    
    [fExt,nIsOnExt,eB,isCorner,eBval] = extraordinary_ext_edges(e,n)
    
    
    # Create a list of the exterior faces that use each node
    fExtUseV = [[-1]]*num_n
    for k in range(numpy.size(fExt,0)):
        for i in range(4):
            if fExtUseV[fExt[k,i]][0]==-1:
                fExtUseV[fExt[k,i]] = [k]
            else:
                fExtUseV[fExt[k,i]].append(k)
    
    # Turn fExtUseV into a matrix by filling it in with indices beyond 
    # the normal size of fExt
    num_ext_f = len(fExt)
    mxLen_fExtUseV = 0
    for k in range(num_n):
        if len(fExtUseV[k])>mxLen_fExtUseV:
            mxLen_fExtUseV = len(fExtUseV[k])
    
    for k in range(num_n):
        for i in range(len(fExtUseV[k]),mxLen_fExtUseV):
            fExtUseV[k].append(num_ext_f)
    
    fExtUseV = numpy.array(fExtUseV)
    
    # Create a list of the elements that use each node
    elUseV = [[-1]]*num_n
    for k in range(numpy.size(e,0)):
        for i in range(8):
            if elUseV[e[k,i]][0]==-1:
                elUseV[e[k,i]] = [k]
            else:
                elUseV[e[k,i]].append(k)
    
    # Turn elUseV into a matrix by filling it in with indices beyond 
    # the normal size of e
    mxLen_elUseV = 0
    for k in range(num_n):
        if len(elUseV[k])>mxLen_elUseV:
            mxLen_elUseV = len(elUseV[k])
    
    for k in range(num_n):
        for i in range(len(elUseV[k]),mxLen_elUseV):
            elUseV[k].append(num_e)
    
    elUseV = numpy.array(elUseV)
    
    
    # Extend the exterior faces, elements, and nodes
    fExt = numpy.concatenate((fExt,numpy.array([[num_n,num_n,num_n,num_n]])))
    e = numpy.concatenate((e,numpy.array([[num_n,num_n,num_n,num_n,num_n,num_n,num_n,num_n]])))
    n = numpy.concatenate((n,numpy.array([[0,0,0]],'float32')))
    nIsOnExt = numpy.concatenate((nIsOnExt,numpy.array([0],'bool')))
    
    
    meanVecUpdate = numpy.zeros((itt,1),float)
    
    nExt = n[numpy.nonzero(nIsOnExt)]
    
    isExtCorner = isCorner[numpy.nonzero(nIsOnExt[0:-1])]
    
    num_nExt = numpy.size(nExt,0)
    
    fExtUseV = fExtUseV[numpy.nonzero(nIsOnExt[0:-1])]
    
    eBExt = eB.copy()
    eBExt = eBExt[numpy.nonzero(nIsOnExt[0:-1])]
    
    p_Bx = p_Bx[nIsOnExt,:,:]
    p_By = p_By[nIsOnExt,:,:]
    p_Bz = p_Bz[nIsOnExt,:,:]

    p_u = p_u[nIsOnExt]
    p_w = p_w[nIsOnExt]
    p_ptchNum = p_ptchNum[nIsOnExt]

    p_face = p_face[nIsOnExt,:]
    p_Nei = ptch_Nei[p_ptchNum,:]
      
    p_isOnUedge = numpy.zeros((p_u.shape),'bool')
    p_isOnWedge = numpy.zeros((p_u.shape),'bool')
    for k in range(p_u.shape[0]):
        if eBExt[k,0]>=0:
            bnd1 = eBExt[k,:]
            bnd1 = bnd1.reshape((bnd1.size,1))
            bnd2 = eB[eBExt[k,:],:]
            bnd2 = bnd2.reshape((bnd2.size,1))
            bnd3 = eB[eB[eBExt[k,:],:],:]
            bnd3 = bnd3.reshape((bnd3.size,1))
            bnd = numpy.concatenate((bnd1,bnd2,bnd3))
            if len(numpy.setdiff1d(p_face[k,[0,1]],bnd))==0:
                p_isOnUedge[k] = 1
            if len(numpy.setdiff1d(p_face[k,[3,2]],bnd))==0:
                p_isOnUedge[k] = 1
            if len(numpy.setdiff1d(p_face[k,[0,3]],bnd))==0:
                p_isOnWedge[k] = 1
            if len(numpy.setdiff1d(p_face[k,[1,2]],bnd))==0:
                p_isOnWedge[k] = 1
        
    for k_itt in range(itt):
        # HexDUAL
        alpha = 1.0

        p_new = hex_dual_non_iso(e[0:-1,:],n[0:-1,:])
        
        mf = alpha*p_new[nIsOnExt[0:-1],:] + (1.0-alpha)*n[nIsOnExt[0:-1],:]
        
        vec = mf-nExt
        normal = vert_normal(fExt,n)
        normal = normal[numpy.nonzero(nIsOnExt)]
        normal[numpy.nonzero(numpy.isnan(normal))] = 0
        
        vecTangent = vec - ((vec*normal).sum(1)).reshape((num_nExt,1))*normal
        vecNormalMag = ((vec*normal).sum(1)).reshape((num_nExt,1))

        # Need to update vecTangent for the nodes on boundary curves
        indBndCrv = numpy.nonzero(numpy.any((-(eB==-1)),axis=1))
        p = n[indBndCrv]
        pplus = n[eB[indBndCrv,0]].squeeze()
        pminus = n[eB[indBndCrv,1]].squeeze()
        
        tminus = (.1/numpy.sqrt(numpy.sum((p-pminus)**2,axis=1))).reshape((p.shape[0],1))
        tplus =  (.1/numpy.sqrt(numpy.sum((p-pplus)**2,axis=1))).reshape((p.shape[0],1))
        
        tminus[numpy.nonzero(numpy.isnan(tminus))]=0
        tplus[numpy.nonzero(numpy.isnan(tplus))]=0
        
        crvTangent = tplus*pplus + (1-tplus)*p - ((1-tminus)*p+tminus*pminus)
        normalize_v3(crvTangent)
        indBndCrv = numpy.nonzero(numpy.any((-(eBExt==-1)),axis=1))
        vecTangent[indBndCrv,:] = ((vec[indBndCrv]*crvTangent).sum(1)).reshape((numpy.size(crvTangent,0),1))*crvTangent
        
        vecTangent[numpy.nonzero(isExtCorner)[0],:] = 0
        
        if immobilizeRidges:
            vecTangent[indBndCrv,:] = 0.0
        
        mf = nExt + vecTangent
        
        #######################################################################        
        [nExt,
         p_face_next,
         p_u_next,
         p_w_next,
         p_Bx_next,
         p_By_next,
         p_Bz_next,
         p_Nei_next] = move_in_surface(nExt.copy(),
                                       mf.copy(),
                                       p_isOnUedge.copy(),
                                       p_isOnWedge.copy(),
                                       p_face.copy(),
                                       p_u.copy(),
                                       p_w.copy(),
                                       p_Bx.copy(),
                                       p_By.copy(),
                                       p_Bz.copy(),
                                       p_Nei.copy(),
                                       ptch_face.copy(),
                                       ptch_Bx.copy(),
                                       ptch_By.copy(),
                                       ptch_Bz.copy(),
                                       ptch_Nei.copy())
        
        # Allow for motion in the normal direction within a tolerance zone
        for k in range(num_nExt):
            if abs(vecNormalMag[k]) > r:
                vecNormalMag[k] = numpy.sign(vecNormalMag[k])*r
        
        vecNormal = vecNormalMag*normal
        nExt = nExt + vecNormal
        
        meanVecUpdate[k_itt] = (numpy.sqrt(((n[nIsOnExt,:] - nExt)**2).sum(1))).mean()
        print("\nIteration: %s \t Mean Vector Update: %s" % (k_itt, meanVecUpdate[k_itt]))
        
        n_next = n.copy()
        n_next[nIsOnExt,:] = nExt
        # Update the nodal position of internal nodes
        n_next[numpy.nonzero(-nIsOnExt[0:-1])] = alpha*p_new[numpy.nonzero(-nIsOnExt[0:-1])] + (1.0-alpha)*n[numpy.nonzero(-nIsOnExt[0:-1])]
        
        # Perform invalid element handling
        n = invalid_element_handling(e[0:-1,:],n_next.copy(),n.copy())
        
        onValidElement = numpy.all(n[nIsOnExt] == n_next[nIsOnExt],axis=1)
        p_face[onValidElement] = p_face_next[onValidElement]
        p_u[onValidElement] = p_u_next[onValidElement]
        p_w[onValidElement] = p_w_next[onValidElement]
        p_Bx[onValidElement] = p_Bx_next[onValidElement]
        p_By[onValidElement] = p_By_next[onValidElement]
        p_Bz[onValidElement] = p_Bz_next[onValidElement]
        p_Nei[onValidElement] = p_Nei_next[onValidElement]
        
        nExt = n[nIsOnExt,:]
        
    n = n[0:num_n,:]
    return n

