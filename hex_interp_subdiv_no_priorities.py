# -*- coding: utf-8 -*-
"""
Performs linear subdivision on the interior of hexahedral elements and 
interpolatory subdivision on the exterior faces.

e = num_e x 8 array of the elements
n = num_n x 3 array of the vertex coordinates

Note that as a result of the linear subdivision on the interior and the C1 
subdivision on the exterior some of the elements may be inverted.  This will be 
fixed in the regularization.

Revisions:
1.1
    Additional call to adjust the nodes on extraordinary exterior edges based on
    and approximating scheme.  This moves the nodes on the extraordinary exterior 
    edges to approximate C1 (fillets these edges).  
    
1.2
    Add priorities and MatList

1.3
    Perform a spatial mapping from a linear hex subdivision to fit the exterior 
    points determined from an interpolatory subdivision surface.  A thin-plate-
    splines mapping is used to update the intermediate nodes.


1.4
    Fix bug which allowed the parent nodes to move in the daughter during the 
    thinplatesplinesSMA()

1.5
    Apply the thin plate splines element by element
    
@author: Greg Sturgeon
gregorymsturgeon@hotmail.com
December 2, 2013
"""
import numpy
from numpy import zeros
from numpy import linalg
from hexblender.interp_subdiv_on_bnd_crv import interp_subdiv
from hexblender.adjust_extraordinary_ext_edge_nodes import adjust_extraordinary_ext_edge_nodes

#from tpsSMA import tpsSMA

def hex_interp_subdiv_no_priorities(e,
                                    n, 
                                    MatList=None,
                                    priorities=None,
                                    allowRidgeMotion=False,
                                    thinPlateMapping=False,
                                    thinPlateRegions=None,
                                    LinearSubdivOnly=False):
    
    print("!! Using hex_interp_subdiv_no_priorities")
    if e.min()==1:
        e = e-1
    
    num_n = n.shape[0]
    num_e = e.shape[0]
    se = numpy.array([[0,1],[1,2],[2,3],[3,0],[4,5],[5,6],[6,7],[7,4],[0,4],[1,5],[2,6],[3,7]])
    sf = numpy.array([[0,1,2,3],[4,5,6,7],[0,1,4,5],[2,3,6,7],[0,3,4,7],[1,2,5,6]])
    eM = {}
    fM = {}
    eP = {}
    fP = {}
    cP = zeros((num_e,3))
    for k in range(num_e):
        for i in range(12):
            eInd = (min(e[k][se[i][0]],e[k][se[i][1]]),max(e[k][se[i][0]],e[k][se[i][1]]))
            eP[eInd] = .5*(n[e[k][se[i][0]]] + n[e[k][se[i][1]]])
        for i in range(6):
            face = [e[k][sf[i][0]],e[k][sf[i][1]],e[k][sf[i][2]],e[k][sf[i][3]]]
            face.sort()
            fInd = tuple(face)
            fP[fInd] = .25*(n[face[0]]+n[face[1]]+n[face[2]]+n[face[3]])
        cP[k] = (1.0/8.0)*n[e[k]].sum(0)
    
    edges = list(eP.keys())
    edgePoints = numpy.array(list(eP.values()))
    num_edges = len(edges)
    
    faces = list(fP.keys())
    facePoints = numpy.array(list(fP.values()))
    num_faces = len(faces)
    
    nN = numpy.concatenate((n,edgePoints,facePoints,cP),axis=0)
    
    cPInd = num_n+num_edges+num_faces
    for k in range(num_edges):
        eM[edges[k]] = num_n+k
    
    for k in range(num_faces):
        fM[faces[k]] = num_n+num_edges+k
    
    for k in range(num_e):
        eInd = list([0]*12)
        for i in range(12):
            eInd[i] = (min(e[k][se[i][0]],e[k][se[i][1]]),max(e[k][se[i][0]],e[k][se[i][1]]))
        fInd = list([0]*6)
        for i in range(6):
            face = [e[k][sf[i][0]],e[k][sf[i][1]],e[k][sf[i][2]],e[k][sf[i][3]]]
            face.sort()
            fInd[i] = tuple(face)
        eN1 = numpy.array([e[k,0],eM[eInd[0]],fM[fInd[0]],eM[eInd[3]],eM[eInd[8]],fM[fInd[2]],k+cPInd,fM[fInd[4]]])
        eN2 = numpy.array([eM[eInd[0]],e[k,1],eM[eInd[1]],fM[fInd[0]],fM[fInd[2]],eM[eInd[9]],fM[fInd[5]],k+cPInd])
        eN3 = numpy.array([fM[fInd[0]],eM[eInd[1]],e[k,2],eM[eInd[2]],k+cPInd,fM[fInd[5]],eM[eInd[10]],fM[fInd[3]]])
        eN4 = numpy.array([eM[eInd[3]],fM[fInd[0]],eM[eInd[2]],e[k,3],fM[fInd[4]],k+cPInd,fM[fInd[3]],eM[eInd[11]]])
        eN5 = numpy.array([eM[eInd[8]],fM[fInd[2]],k+cPInd,fM[fInd[4]],e[k,4],eM[eInd[4]],fM[fInd[1]],eM[eInd[7]]])
        eN6 = numpy.array([fM[fInd[2]],eM[eInd[9]],fM[fInd[5]],k+cPInd,eM[eInd[4]],e[k,5],eM[eInd[5]],fM[fInd[1]]])
        eN7 = numpy.array([k+cPInd,fM[fInd[5]],eM[eInd[10]],fM[fInd[3]],fM[fInd[1]],eM[eInd[5]],e[k,6],eM[eInd[6]]])
        eN8 = numpy.array([fM[fInd[4]],k+cPInd,fM[fInd[3]],eM[eInd[11]],eM[eInd[7]],fM[fInd[1]],eM[eInd[6]],e[k,7]])
        eN1.resize(1,8)
        eN2.resize(1,8)
        eN3.resize(1,8)
        eN4.resize(1,8)
        eN5.resize(1,8)
        eN6.resize(1,8)
        eN7.resize(1,8)
        eN8.resize(1,8)
        if k==0:
            eN = numpy.concatenate((eN1,eN2,eN3,eN4,eN5,eN6,eN7,eN8),axis=0)
        else:
            eN = numpy.concatenate((eN,eN1,eN2,eN3,eN4,eN5,eN6,eN7,eN8),axis=0)
    
    if LinearSubdivOnly:
        return list([eN,nN])
    
    nN_linear = nN.copy()
    
    # InterpSubDiv
    #[fPn,ePn] = interpSubDiv(e,n)
    [fPn,ePn,nOnBndCrv] = interp_subdiv(e,n,returnnOnBndCrv=True)
    fPn_keys = list(fPn.keys())
    fPn_values = list(fPn.values())
    
    interpP = []
    for k in range(len(fPn_keys)):
        nN[fM[fPn_keys[k]]] = fPn_values[k]
        interpP.append(fM[fPn_keys[k]])
    
    ePn_keys = list(ePn.keys())
    ePn_values = list(ePn.values())
    for k in range(len(ePn_keys)):
        nN[eM[ePn_keys[k]]] = ePn_values[k]
        interpP.append(eM[ePn_keys[k]])
    
    interpP = numpy.array(interpP)
    
    if thinPlateMapping:
        nN = applyTPSbyElement(eN,nN,nN_linear)
#        interpP = numpy.setdiff1d(interpP,numpy.arange(0,numpy.size(n,0)))
#        ctrlpoints = nN[interpP]
#        points = nN_linear[interpP]
#        
#        nN_linearWarped = tpsSMA(points,ctrlpoints,nN_linear)
#        
#        eP = numpy.arange(num_n,num_n+num_edges)
#        
#        eP_linear = numpy.setdiff1d(eP,interpP)
#        
#        b = numpy.setdiff1d(numpy.arange(0,numpy.size(nN,0)),interpP)
#        b = numpy.setdiff1d(b,numpy.arange(0,num_n))
#        nN[b,:] = nN_linearWarped[b,:]
#        nN[eP_linear,:] = nN_linear[eP_linear,:]
    
    if allowRidgeMotion:
        nN = adjust_extraordinary_ext_edge_nodes(eN,nN)
    
    return list([eN,nN])



def tpsSMA(points,ctrlpoints,toWarp):
    npnts = points.shape[0]
    a = numpy.reshape(numpy.sum((points.copy()**2),1),[1,npnts])
    k = numpy.dot(numpy.ones([npnts,1]),a).T + numpy.dot(numpy.ones([npnts,1]), a) - 2.0*numpy.dot(points,points.T)
    
    k[k<1e-320] = 1e-320
    
    k = numpy.sqrt(k)
    # Calculate P matrix
    p = numpy.concatenate((numpy.ones((npnts,1)),points.copy()),1)
    
    # Calculate L matrix
    l = numpy.concatenate((numpy.concatenate((k,p),1),numpy.concatenate((p.T,numpy.zeros((4,4))),1)),0)
    
    param = numpy.dot( linalg.pinv(l) , numpy.concatenate((ctrlpoints,numpy.zeros((4,3))),0) )
    
    
    # Calculate new coordinates (x',y',z') for each points 
    
    pntsNum = toWarp.shape[0]
     
    k = numpy.zeros((pntsNum,npnts))
    
    gx = toWarp[:,0]
    gy = toWarp[:,1]
    gz = toWarp[:,2]
    
    for nn in range(npnts):
        k[:,nn] = (gx - points[nn,0])**2 + (gy - points[nn,1])**2 + (gz - points[nn,2])**2
    
    k[k<1e-320] = 1e-320
    
    k = numpy.sqrt(k)
    
    p = numpy.concatenate((numpy.ones((pntsNum,1)),toWarp.copy()),1)
    
    l = numpy.concatenate((k,p),1)
    
    warped = numpy.dot( l, param )
    return warped


def applyTPSbyElement(eN,nN,nN_linear):
    nN_TPS = numpy.zeros(numpy.shape(nN))
    
    nN_TPS_count = numpy.zeros(numpy.shape(nN)[0])
    
    
    for k in range(0,int(numpy.shape(eN)[0]/8)):
        curNodes = numpy.unique(numpy.array([eN[8*k,:], eN[8*k+1,:], eN[8*k+2,:], eN[8*k+3,:], eN[8*k+4,:], eN[8*k+5,:], eN[8*k+6,:], eN[8*k+7,:]]))
        
        eCurLin = nN_linear[curNodes,:]
        eCurInterp = nN[curNodes,:]
        
        iInterp = curNodes[(numpy.any((eCurLin!=eCurInterp),axis=1))]
        
        # Also need to include the corner nodes from one level coarser
        iCorner = numpy.array([eN[8*k,0],eN[8*k+1,1], eN[8*k+2,2], eN[8*k+3,3], eN[8*k+4,4], eN[8*k+5,5], eN[8*k+6,6], eN[8*k+7,7]])
        
        iCP = numpy.concatenate([iInterp,iCorner],0)
        points = nN_linear[iCP,:]
        ctrlpoints = nN[iCP,:]
        
        iLin = numpy.setdiff1d(curNodes,iCP)
        obj = nN_linear[iLin,:]
        
        wobj = tpsSMA(points, ctrlpoints,obj)
        
        nN_TPS[iLin,:] += wobj
        nN_TPS_count[iLin] += 1
    
    
    
    nN_TPS = numpy.vstack([nN_TPS[:,0]/nN_TPS_count, nN_TPS[:,1]/nN_TPS_count, nN_TPS[:,2]/nN_TPS_count]).T
    nN_TPS[nN_TPS_count==0,:] = nN[nN_TPS_count==0,:]
    
    return nN_TPS


