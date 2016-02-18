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
    ThinPlateSplines()
    
1.4.1
    Changed linalg package to use numpy instead of scipy since the only function 
    using it is the pseudoinverse pinv. Hopefully it should be the same

@author: Greg Sturgeon
gregorymsturgeon@hotmail.com
August 20, 2012
Changed : Adarsh Krishnamurthy
April 8, 2013

"""

import time
import numpy
import collections
from hexblender.interp_subdiv import interp_subdiv
from numpy import linalg
from numpy.matlib import repmat, repeat
from hexblender.adjust_extraordinary_ext_edge_nodes import adjust_extraordinary_ext_edge_nodes

def hex_interp_subdiv(e, 
                      n, 
                      MatList=None, 
                      priorities=None,
                      allowRidgeMotion=False,
                      thinPlateMapping=False,
                      thinPlateRegions=None,
                      LinearSubdivOnly=False):
    
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
    cP = numpy.zeros((num_e,3))
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
    if priorities:
        interpP = []
        priorities.reverse() #a method for lists. Go backward through to execute on lower priority regions first at the 'back' of the list
        for matGroup in priorities:
            currElemsInds = numpy.array([ind for ind,val in enumerate(MatList) if val in matGroup])
            allVertInds = numpy.unique(e[currElemsInds,:]).flatten()
            allVertIndsSet = set(allVertInds)
            isVertIncluded = numpy.array([vert in allVertIndsSet for vert in range(numpy.shape(n)[0])])
            vertSelectionMap = numpy.cumsum(isVertIncluded)-1
            
            eReduced = vertSelectionMap[e[currElemsInds,:]]
            nReduced = numpy.zeros([numpy.shape(allVertInds)[0],3])
            for ind,val in enumerate(nReduced):
                nReduced[ind,0] = n[allVertInds[ind],0]
                nReduced[ind,1] = n[allVertInds[ind],1]
                nReduced[ind,2] = n[allVertInds[ind],2]
                
            [fPn,ePn] = interp_subdiv(eReduced,nReduced)
            fPn_keys = list(fPn.keys())
            fPn_values = list(fPn.values())
            
            #need to convert the "real" vertex numbers to the "reduced" index numbers in the keys
            fM_keys = list(fM.keys())
            fM_values = list(fM.values())
            fM_keys_reduced = [newList for newList in fM_keys if all(myNewVal in allVertIndsSet for myNewVal in newList)]
            indicesForFMvalues = [ind for ind, newList in enumerate(fM_keys) if all(myNewVal in allVertIndsSet for myNewVal in newList)]
            fM_values_reduced = [fM_values[val] for val in indicesForFMvalues]
            fM_keys_reduced_mapped = [[vertSelectionMap[fullNodeNum] for fullNodeNum in myList] for myList in fM_keys_reduced]
            
            fM_reduced = {}
            for ind, key in enumerate(fM_keys_reduced_mapped):
                fM_reduced[tuple(key)] = fM_values_reduced[ind]
            
            eM_keys = list(eM.keys())
            eM_values = list(eM.values())
            eM_keys_reduced = [newList for newList in eM_keys if all(myNewVal in allVertIndsSet for myNewVal in newList)]
            indicesForEMvalues= [ind for ind, newList in enumerate(eM_keys) if all(myNewVal in allVertIndsSet for myNewVal in newList)]
            eM_values_reduced = [eM_values[val] for val in indicesForEMvalues]
            eM_keys_reduced_mapped = [[vertSelectionMap[fullNodeNum] for fullNodeNum in myList] for myList in eM_keys_reduced]
            eM_reduced = {}
            for ind, key in enumerate(eM_keys_reduced_mapped):
                eM_reduced[tuple(key)] = eM_values_reduced[ind]
            
            for k in range(len(fPn_keys)):
                nN[fM_reduced[fPn_keys[k]]] = fPn_values[k]
                interpP.append(fM_reduced[fPn_keys[k]])
            ePn_keys = list(ePn.keys())
            ePn_values = list(ePn.values())
            for k in range(len(ePn_keys)):
                nN[eM_reduced[ePn_keys[k]]] = ePn_values[k]
                interpP.append(eM_reduced[ePn_keys[k]])

    else:
        [fPn,ePn] = interp_subdiv(e,n)
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
    interpP = numpy.concatenate((interpP,numpy.array(list(range(num_n)))),0)
    
    if thinPlateMapping:
        if MatList and thinPlateRegions: #only execute by region if MatList and thinPlateRegions are present
            
            thinPlateRegions.reverse() #a method for lists. Go backward through to execute on lower priority regions first at the 'back' of the list
            MatListOnceSubdiv = sum([8*[val] for val in MatList],[]) #sum (X,[]) flattens the list of lists made by the list comprehension
            for matGroup in thinPlateRegions:
                
                origElemsInds = numpy.array([ind for ind,val in enumerate(MatList) if val in matGroup])
                currElemsInds = numpy.array([ind for ind,val in enumerate(MatListOnceSubdiv) if val in matGroup])
                allVertInds = numpy.unique(eN[currElemsInds,:]).flatten()
                allVertIndsSet = set(allVertInds)
                isVertIncluded = numpy.array([vert in allVertIndsSet for vert in range(numpy.shape(nN)[0])])
                vertSelectionMap = numpy.cumsum(isVertIncluded)-1
                
                #interpP here plays the role of "allVertInds". Need this step and next to be consistent with vertex numberings determined in the lines above
                interpPvertsIncluded = numpy.array([vertNum for vertNum in interpP if vertNum in allVertIndsSet])
                interpPReduced = vertSelectionMap[interpPvertsIncluded]
                
                eReduced = vertSelectionMap[e[origElemsInds,:]]
                eNReduced = vertSelectionMap[eN[currElemsInds,:]]
                nNReduced = numpy.zeros([numpy.shape(allVertInds)[0],3])
                nN_linReduced = numpy.zeros([numpy.shape(allVertInds)[0],3])

                for ind,val in enumerate(nNReduced):
                    nNReduced[ind,0] = nN[allVertInds[ind],0]
                    nNReduced[ind,1] = nN[allVertInds[ind],1]
                    nNReduced[ind,2] = nN[allVertInds[ind],2]
                for ind,val in enumerate(nN_linReduced):
                    nN_linReduced[ind,0] = nN_linear[allVertInds[ind],0]
                    nN_linReduced[ind,1] = nN_linear[allVertInds[ind],1]
                    nN_linReduced[ind,2] = nN_linear[allVertInds[ind],2]
            
                nN_transformed = ThinPlateSplines(eReduced,nNReduced,nN_linReduced,interpPReduced)
                
                for ind,vert_num in enumerate(allVertInds):
                    nN[vert_num,0] = nN_transformed[ind,0]
                    nN[vert_num,1] = nN_transformed[ind,1]
                    nN[vert_num,2] = nN_transformed[ind,2]
            
        else: #execute everything by default if there is no MatList
            nN = ThinPlateSplines(e,nN,nN_linear,interpP)
    
    if allowRidgeMotion:
        nN = adjust_extraordinary_ext_edge_nodes(eN, nN)
    
    return list([eN,nN])



def ThinPlateSplines(e,nN,nN_linear,interpP):
    # Adjust linear by thin plate splines
    now = time.time()
    num_e = e.shape[0]
    
    ctrlpoints = nN[numpy.concatenate((numpy.array(list(range(num_e))),interpP),0),:]
    points = nN_linear[numpy.concatenate((numpy.array(list(range(num_e))),interpP),0),:]
    
    npnts = points.shape[0]
    
    distMat = numpy.sum((repmat(points, npnts, 1) - repeat(points, npnts, axis=0))**2, axis=1)
    k = distMat.reshape((npnts, npnts))

    k[k<1e-320] = 1e-320
    
    k = numpy.sqrt(k)
    # Calculate P matrix
    p = numpy.concatenate((numpy.ones((npnts,1)),points.copy()),1)
    
    # Calculate L matrix
    l = numpy.concatenate((numpy.concatenate((k,p),1),numpy.concatenate((p.T,numpy.zeros((4,4))),1)),0)
    
    param = numpy.dot( linalg.pinv(l) , numpy.concatenate((ctrlpoints,numpy.zeros((4,3))),0) )
    # Calculate new coordinates (x',y',z') for each points 
    
    pntsNum = nN_linear.shape[0]
     
    k = numpy.zeros((pntsNum,npnts))
    
    gx = nN_linear[:,0]
    gy = nN_linear[:,1]
    gz = nN_linear[:,2]
    
    for nn in range(npnts):
        k[:,nn] = (gx - points[nn,0])**2 + (gy - points[nn,1])**2 + (gz - points[nn,2])**2
    
    k[k<1e-320] = 1e-320
    
    k = numpy.sqrt(k)
    
    p = numpy.concatenate((numpy.ones((pntsNum,1)),nN_linear.copy()),1)
    
    l = numpy.concatenate((k,p),1)
    
    wnN_linear = numpy.dot( l, param )
    
    wnN_linear[numpy.concatenate((numpy.array(list(range(num_e))),interpP),0),:] = ctrlpoints
    print("ThinPlateSplines took: %f" % (time.time() - now))
    
    return wnN_linear
    
