"""
MoveInSurface(nExt,mf,p_isOnUedge,p_isOnWedge,p_face,p_u,p_w,p_Bx,p_By,p_Bz,p_Nei,ptch_face,ptch_Bx,ptch_By,ptch_Bz,ptch_Nei):

nExt - nodes on the Exterior surfaces
mf - motion vector of nExt
p_isOnUedge - boolean identifying if the node is on an extraordinary exterior edge in the u parameter direction
p_isOnWedge - boolean identifying if the node is on an extraordinary exterior edge in the w parameter direction
p_face - identifies the node ordering of the patch that the node is in (as rotated)
p_u - u parameter for position in current patch
p_w - w parameter for position in current patch
p_Bx,p_By,p_Bz - matrices of the Hermite derivatives for the current patch
p_Nei - identifies the neighboring patches to the current patch
ptch_face - stores the node order for each patch face
ptch_Bx, ptch_By, ptch_Bz - stores the Hermite derivatives for each patch
ptch_Nei - stores the neighboring patches of each patch

MoveInSurface 1) pulls the motion vector back from 3-space to the parameter plane 
to obtain vTu and vTw. 2) determines if the motion causes a node to cross from one
patch to another.  3) when a node crosses a patch the neighboring patch that it 
crosses into must be identified and then rotated to keep consistent 
parameterization directions (a step in u continues in u). 4) when a node is on 
a U or W edge is is constrained to stay on that edge.

January 17, 2012
Greg Sturgeon
"""
import numpy
from hexblender.regularize_elements_pullback_functions import hermite_points
from hexblender.rotate_patch_vectorized import rotate_patch_new
from hexblender.uw_component import uw_component


def move_in_surface(nExt, mf, p_isOnUedge, p_isOnWedge, p_face, p_u, p_w, p_Bx, 
                    p_By, p_Bz, p_Nei, ptch_face, ptch_Bx, ptch_By,
                    ptch_Bz, ptch_Nei):
    
    indPrev =  list(range(p_u.shape[0]))
    indPrev = numpy.array(indPrev,'int')
    
    nExtOut = nExt.copy()
    p_faceOut = p_face.copy()
    p_uOut = p_u.copy()
    p_wOut = p_w.copy()
    p_BxOut = p_Bx.copy()
    p_ByOut = p_By.copy()
    p_BzOut = p_Bz.copy()
    p_NeiOut = p_Nei.copy()
    
    numitt = 3
    for itt in range(numitt):
        [vTu,vTw] = uw_component(mf,nExt,p_Bx, p_By, p_Bz, p_u, p_w)
        
        # Check if on a u or w edge
        vTu[p_isOnWedge] = 0
        vTw[p_isOnUedge] = 0
        
        crossedEdge = numpy.zeros((vTu.shape[0],4),'bool')
        crossedEdge[:,0] = (p_w + vTw)<0
        crossedEdge[:,1] = (p_u + vTu)>1
        crossedEdge[:,2] = (p_w + vTw)>1
        crossedEdge[:,3] = (p_u + vTu)<0
        
        stayInPtch = -(crossedEdge.any(axis=1))
        
        # Apply the motion for the points that stay in the patch
        notCrossedUedge = numpy.logical_and(-crossedEdge[:,1],-crossedEdge[:,3])
        notCrossedWedge  = numpy.logical_and(-crossedEdge[:,0],-crossedEdge[:,2])
        p_u[notCrossedUedge ] = (p_u[notCrossedUedge] + vTu[notCrossedUedge])
        p_w[notCrossedWedge] = (p_w[notCrossedWedge] + vTw[notCrossedWedge])
        
        # Move to the edge of the patch before crossing it.
        p_w[crossedEdge[:,0]] = 0
        p_u[crossedEdge[:,1]] = 1
        p_w[crossedEdge[:,2]] = 1
        p_u[crossedEdge[:,3]] = 0
        
        # Check to find the patch it moves into
        # If there is no neighboring patch then it is on a boundary (set crossedEdge accordingly)
        iEN=[1,3,2,0]
        uwdir = numpy.array(([[1,2],[0,3],[3,2],[0,1]]),'int')
        for j in range(4):
            edgeN = -numpy.ones((vTu.shape[0]),'int')
            edgeN[crossedEdge[:,iEN[j]]] = p_Nei[crossedEdge[:,iEN[j]],iEN[j]]
            crossedEdge[edgeN==-1,iEN[j]] = 0
            iCE = numpy.nonzero(crossedEdge[:,iEN[j]])[0]
            edgeCrossed = numpy.array([p_face[iCE,uwdir[j,0]],p_face[iCE,uwdir[j,1]]])
            [p_face[iCE,:],
             p_Bx[iCE,:,:],
             p_By[iCE,:,:],
             p_Bz[iCE,:,:],
             p_Nei[iCE,:]] = rotate_patch_new(
                                ptch_face[edgeN[iCE]].copy(),
                                edgeCrossed.transpose(),
                                uwdir[j,:],
                                ptch_Bx[edgeN[iCE],:,:].copy(),
                                ptch_By[edgeN[iCE],:,:].copy(),
                                ptch_Bz[edgeN[iCE],:,:].copy(),
                                ptch_Nei[edgeN[iCE],:].copy())
            
            if iEN[j]==1:
                p_u[iCE] = 0
            if iEN[j]==3:
                p_u[iCE] = 1
            if iEN[j]==2:
                p_w[iCE] = 0
            if iEN[j]==0:
                p_w[iCE] = 1
        
        nExt = hermite_points(p_Bx, p_By, p_Bz, p_u, p_w)
        mf[stayInPtch] = nExt[stayInPtch]
        
        
        # Fill the points that stay in the patch into the output
        nExtOut[indPrev[stayInPtch],:] = nExt[stayInPtch,:] .copy()
        p_faceOut[indPrev[stayInPtch],:] = p_face[stayInPtch,:] .copy()
        p_uOut[indPrev[stayInPtch]] = p_u[stayInPtch].copy()
        p_wOut[indPrev[stayInPtch]] = p_w[stayInPtch].copy()
        p_BxOut[indPrev[stayInPtch],:,:] = p_Bx[stayInPtch,:,:].copy()
        p_ByOut[indPrev[stayInPtch],:,:] = p_By[stayInPtch,:,:].copy()
        p_BzOut[indPrev[stayInPtch],:,:] = p_Bz[stayInPtch,:,:].copy()
        p_NeiOut[indPrev[stayInPtch],:] = p_Nei[stayInPtch,:].copy()
        
        indPrev = indPrev[-stayInPtch]
        
        nExt = nExt[-stayInPtch,:]
        mf = mf[-stayInPtch,:]
        
        p_isOnUedge = p_isOnUedge[-stayInPtch]
        p_isOnWedge = p_isOnWedge[-stayInPtch]
        p_face = p_face[-stayInPtch,:]
        p_u = p_u[-stayInPtch]
        p_w = p_w[-stayInPtch]
        p_Bx = p_Bx[-stayInPtch,:,:]
        p_By = p_By[-stayInPtch,:,:]
        p_Bz = p_Bz[-stayInPtch,:,:]
        p_Nei = p_Nei[-stayInPtch,:]
        
        if stayInPtch.all():
            break
        
        if itt == numitt:
            nExtOut[indPrev,:] = nExt.copy()
            p_faceOut[indPrev,:] = p_face.copy()
            p_uOut[indPrev] = p_u.copy()
            p_wOut[indPrev] = p_w.copy()
            p_BxOut[indPrev,:,:] = p_Bx.copy()
            p_ByOut[indPrev,:,:] = p_By.copy()
            p_BzOut[indPrev,:,:] = p_Bz.copy()
            p_NeiOut[indPrev,:] = p_Nei.copy()
    
    return (nExtOut,p_faceOut,p_uOut,p_wOut,p_BxOut,p_ByOut,p_BzOut,p_NeiOut)
