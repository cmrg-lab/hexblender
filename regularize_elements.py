    # -*- coding: utf-8 -*-
"""
Regularize hexahedral elements to improve the element skew, uniformity of 
element size, and aspect ratios.
  
Inputs:
e = num_e x 8 array of the elements
n = num_n x 3 array of the vertex coordinates
itt = number of itterations to perform

Based on steps 2 & 3 of ‘Surface Smoothing and Quality Improvement of 
Quadrilateral / Hexahedral Meshes with Geometric Flow’ by Zhang et.al

Nodes on the surface of the mesh are moved in the tangent plane.  Interior nodes
are repositioned to a weighted average of their neighboring elements centroids.

Since the motion is in the tangent plane this algorithm will have minimal effect
on the objects shape if the mesh well approximates smooth objects.

Rev 1.0
Correct for division by zeros warnings

Rev 1.1
Improve performance and memory usage by only updating the interior nodes every 
fifth itteration (initially) or every itteration for the final ten itterations.  
float32 was also used instead of float64

Rev 1.2
Correct for a bug where vertex 0 may remain stationary durring the regularization.
vecTangent[numpy.nonzero(isExtCorner),:] = 0   ->  vecTangent[numpy.nonzero(isExtCorner)[0],:] = 0
n[numpy.nonzero(nIsOnExt)] = nExt  ->  n[numpy.nonzero(nIsOnExt)[0]] = nExt
n[numpy.nonzero(-nIsOnExt)] = nInt  ->  n[numpy.nonzero(-nIsOnExt)[0]] = nInt

Rev 1.3
Include a function to determine the valance of extraordinary exterior edges.  The edges with valance ~= 2
are included in the eB array which identifies the nodes on an extraordinary exterior edges and provides its 
neighboring nodes.  In the event of a Y junction formed by edges with different valance the neighbors the 
neighbors with the highest valance are stored in eB.  This allows for corners at valance 1 and valance 3 edges
to move tangent to the valance 3 edge. 

Rev 1.4
Update the order of fExt to maintain outward facing normals.
    sf = numpy.array([[0,3,2,1],[5,6,7,4],[0,1,5,4],[3,7,6,2],[0,4,7,3],[1,2,6,5]])
was
    sf = numpy.array([[0,1,2,3],[4,5,6,7],[0,4,5,1],[3,2,6,7],[0,3,7,4],[1,5,6,2]])

@author: Greg Sturgeon
gregorymsturgeon@hotmail.com
November 16, 2011
"""

import numpy as np

def regularize_elements(e,
                        n,
                        itt, 
                        preserveRidges = True,
                        immobExtraord = False,
                        tangentMotion = True,
                        normalMotion = True,
                        internalMotion = True):
    
    if e.min()==1:
        e = e-1
    
    if n.dtype.name == 'float64':
        n = np.float32(n)
    
    num_e = np.size(e,0)
    num_n = np.size(n,0)
    
    [fExt,nIsOnExt,eB,isCorner,eBval] = extraordinary_ext_edges(e,n)
    
    
    # Create a list of the exterior faces that use each node
    fExtUseV = [[-1]]*num_n
    for k in range(np.size(fExt,0)):
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
    
    fExtUseV = np.array(fExtUseV)
    maxValence = mxLen_fExtUseV
    
    # Create a list of the elements that use each node
    elUseV = [[-1]]*num_n
    for k in range(np.size(e,0)):
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
    
    elUseV = np.array(elUseV)
    
    #orderedOneNeigh[num_nodes, 2*max_valence]: stores NODE NUMBERS
    #orderedOneNeigh ;= waterfall algorithm orders the one-neighborhood 
    # nodes in the order they are needed for vecNormal calculation
    orderedOneNeigh = [[]]*num_n
    orderedOneNeighFaces = [[] for val in range(num_n)]
    rotatedFaceVertMap = [[] for val in range(num_n)]
    ExtValence = np.zeros([num_n],np.int32)

    # ind := current node number, num_n_list := list of face numbers 
    # for that node
    for ind,num_n_list in enumerate(fExtUseV):  

        #short circuit if val == -1, i.e not an exterior node
        if num_n_list[0] == -1: 
            continue
            
        curr_face_list = num_n_list[np.nonzero(num_n_list < np.size(fExt,0) )[0]].copy()
        curr_node_list = np.array(fExt[curr_face_list]) #ind 0: face instance of node "ind"; ind 1: node within face instance 1-->4
        curr_face_list = curr_face_list.tolist()
        ExtValence[ind] = len(curr_face_list)
        # turned into array because will be square and np.roll() should 
        # be faster than other circshifting strategies
        #~ print("curr face list = %s" % curr_face_list)
        #~ print("curr node list = %s\n" % curr_node_list)
        
        # rotate faces locally until the local node is as the first position
        for ind2,face_instance in enumerate(curr_node_list):
        
            if face_instance[0] == ind:
                continue_while_loop = False
            else:
                continue_while_loop = True
            while continue_while_loop and any(face_instance == ind):
                face_instance = np.roll(face_instance,1)
                if face_instance[0] == ind:
                    continue_while_loop = False
                    curr_node_list[ind2] = face_instance

        # now, curr_node_list, for each node, [face instance, 1:4] is now
        # circshifted with curr_node("ind") placed at position 0 in index 1
        
        #~ print("new node list = %s\n" % curr_node_list)
        
        # store a copy of curr_node_list, with its properly-rotated faces,
        curr_node_list_copy = curr_node_list.copy().tolist()
        rotatedFaceVertMap[ind].append(curr_node_list_copy[0])
        
        curr_unique_node_list = list(np.unique(curr_node_list.flatten()))
        orderedOneNeigh[ind] = list(curr_node_list[0][1:])

        # pluck out first face from the list and then delete it
        orderedOneNeighFaces[ind].append(curr_face_list[0]) 

        # remove first value which was just stored
        curr_face_list.pop(0) 
        curr_unique_node_list = [val for val in curr_unique_node_list if val not in orderedOneNeigh[ind] and val != ind]
        #~ print("orderedOneNeigh = %s" % orderedOneNeigh[ind])
        #~ print("curr unique node list = %s" % curr_unique_node_list)
        # delete first row / column of curr_node_list which have been stored already / are unneeded, respectively
        curr_node_list = np.delete(curr_node_list,0,0)
        curr_node_list = np.delete(curr_node_list,0,1)
        stayInLoop = True
        while stayInLoop:
            hit = np.nonzero(curr_node_list ==orderedOneNeigh[ind][-1])
            if hit[1][0] == 2: #last vertex position in curr_node_list
                row_index = hit[0][0]
                addList = list(curr_node_list[row_index][-2::-1])
                addList2 = list(curr_node_list[row_index][-1::-1])
                orderedOneNeigh[ind].extend(addList)
                orderedOneNeighFaces[ind].append(curr_face_list[row_index])
            elif hit[1][0] == 0: #first vertex position in curr_node_list
                row_index = hit[0][0]
                addList = list(curr_node_list[row_index][1:])
                addList2 = list(curr_node_list[row_index])
                orderedOneNeigh[ind].extend(addList)
                orderedOneNeighFaces[ind].append(curr_face_list[row_index])
            
            # delete row and face no longer needed
            rotatedFaceVertMap[ind].append([ind] + addList2)
            curr_node_list = np.delete(curr_node_list,row_index,0)
            curr_face_list.pop(row_index)
            curr_unique_node_list = [val for val in curr_unique_node_list if val not in addList]
            if not curr_unique_node_list: #i.e, is empty
                stayInLoop = False
        
        # first vertex is double-counted, so pop it out
        orderedOneNeigh[ind].pop(-1)
        #~ print("rotatedFaceVertMap = %s" % rotatedFaceVertMap[ind])
        #~ print("orderedOneNeigh = %s" % orderedOneNeigh[ind])
    
    # turn orderedOneNeigh into a matrix in the same way as is done above
    mxLen_orderedOneNeigh = 0
    for k in range(num_n):
        if len(orderedOneNeigh[k])>mxLen_orderedOneNeigh:
            mxLen_orderedOneNeigh = len(orderedOneNeigh[k])
    
    # append vertex number so that the displacement (p_i - p_k) = 0 for each vertex for k > n
    for k in range(num_n):
        for i in range(len(orderedOneNeigh[k]),mxLen_orderedOneNeigh):
            orderedOneNeigh[k].append(k) 

    orderedOneNeigh = np.array(orderedOneNeigh)
    orderedOneNeigh = orderedOneNeigh[np.nonzero(nIsOnExt)]
    
    # "cumsum - 1" is a way to re-index the array orderedOneNeigh into 
    # the node number in nExt arrays
    inverse_nIsOnExt = np.cumsum(nIsOnExt)-1
    inverse_orderedOneNeigh = inverse_nIsOnExt[orderedOneNeigh]
    
    # turn orderedOneNeighFaces into a matrix with num_n appended at end
    mxLen_orderedOneNeighFaces = 0
    for k in range(num_n):
        if len(orderedOneNeighFaces[k])>mxLen_orderedOneNeighFaces:
            mxLen_orderedOneNeighFaces = len(orderedOneNeighFaces[k])
    
    # append vertex number so that the displacement (p_i - p_k) = 0 
    # for each vertex for k > n
    for k in range(num_n):
        for i in range(len(orderedOneNeighFaces[k]),mxLen_orderedOneNeighFaces):
            orderedOneNeighFaces[k].append(-1)
            
    orderedOneNeighFaces = np.array(orderedOneNeighFaces)
    orderedOneNeighFaces = orderedOneNeighFaces[np.nonzero(nIsOnExt)]
    
    #same for rotatedFaceVertMap
    mxLen_rotatedFaceVertMap = 0
    for k in range(num_n):
        if len(rotatedFaceVertMap[k])>mxLen_rotatedFaceVertMap:
            mxLen_rotatedFaceVertMap = len(rotatedFaceVertMap[k])

    for k in range(num_n):
        for i in range(len(rotatedFaceVertMap[k]), mxLen_rotatedFaceVertMap):
                rotatedFaceVertMap[k].append([num_n,num_n,num_n,num_n])
            
    rotatedFaceVertMap = np.array(rotatedFaceVertMap)
    rotatedFaceVertMap = rotatedFaceVertMap[np.nonzero(nIsOnExt)]
    
    # Extend the exterior faces, elements, and nodes
    fExt = np.concatenate((fExt,np.array([[num_n,num_n,num_n,num_n]])))
    e = np.concatenate((e,np.array([[num_n,num_n,num_n,num_n,num_n,num_n,num_n,num_n]])))
    n = np.concatenate((n,np.array([[0,0,0]],'float32')))
    nIsOnExt = np.concatenate((nIsOnExt,np.array([0],'bool')))
     
    meanVecTangent = np.zeros((itt,1),float)
    
    nExt = n[np.nonzero(nIsOnExt)]
    nInt = n[np.nonzero(-nIsOnExt)]
    num_nExt = np.size(nExt,0)
    
    fExtUseV = fExtUseV[np.nonzero(nIsOnExt[0:-1])]
    elUseVInt = elUseV[np.nonzero(-nIsOnExt[0:-1])]
    
    eBExt = eB.copy()
    eBExt = eBExt[np.nonzero(nIsOnExt[0:-1])]
    
    isExtCorner = isCorner[np.nonzero(nIsOnExt[0:-1])]
    
    indBndCrvArray = np.nonzero(np.any((-(eB==-1)),axis=1))[0]
    isExtCornerArray = np.nonzero(isExtCorner)[0]
    
    # find array locations in orderedOneNeigh that correspond to corners or ridges
    inverse_orderedOneNeigh_normals = inverse_orderedOneNeigh.copy()
    row_inds = np.array([],np.int32)
    col_inds = np.array([],np.int32)
    for ind1, row in enumerate(orderedOneNeigh):
        for ind2, val in enumerate(row):
            if (np.any(val == isExtCornerArray) or np.any(val == indBndCrvArray)):
                row_inds = np.concatenate((row_inds, np.array([ind1])))
                col_inds = np.concatenate((col_inds, np.array([ind2])))
                inverse_orderedOneNeigh_normals[ind1,ind2] = ind1
    
    CornerOrEB1neigh = (row_inds, col_inds)
    
    # look for extraordinary nodes and their one-neighborhood and demand 
    # that they have no tangent motion without an implicit integration 
    # scheme, these areas can be disproportionately unstable since the
    # algorithm does not directly penalize skewed elements
    isValGT5 = (ExtValence != 4)
    isValGT5Ext = isValGT5[np.nonzero(nIsOnExt[0:-1])]
    # concatenate these nodes' one-neighborhoods
    isValGT5Ext = np.unique(np.concatenate((np.nonzero(isValGT5Ext)[0], inverse_orderedOneNeigh[isValGT5Ext].flatten())))
    
    print("\nmeanVecTangent")
    print("---------------") 

    for k_itt in range(itt):
        # Face Points for the exterior faces
        fP = (n[fExt].sum(1))/4.
        
        fAreas = face_area(fExt,n)
        
        mf = np.zeros((np.size(fExtUseV,0),3),'float32')
        w_areas= fAreas[fExtUseV]
        totAreas = (fAreas[fExtUseV].sum(1)).reshape((np.size(fExtUseV,0),1))
        totAreas[np.nonzero(totAreas==0)] = 1.0
        mf[:,0] = (fP[fExtUseV,0]*w_areas).sum(1)
        mf[:,1] = (fP[fExtUseV,1]*w_areas).sum(1)
        mf[:,2] = (fP[fExtUseV,2]*w_areas).sum(1)
        mf /= totAreas
        
        vec = mf-nExt
        normal = vert_normal(fExt,n)
        
        normal = normal[np.nonzero(nIsOnExt)]
        normal[np.nonzero(np.isnan(normal))] = 0.0
        
        vecTangent = vec - ((vec*normal).sum(1)).reshape((num_nExt,1))*normal
        vecNormal = calc_vec_normal(n,nExt,normal,totAreas,inverse_orderedOneNeigh,maxValence,rotatedFaceVertMap,inverse_orderedOneNeigh_normals,preserveRidges)
        
        #get indices of the boundary contours/curves ("ridges")
        indBndCrv = np.nonzero(np.any((-(eB==-1)),axis=1))
        indBndCrv_Ext = inverse_nIsOnExt[indBndCrv[0]]

        #zero out vecNormal motions on nodes that we don't want to move
        vecNormal[np.nonzero(isExtCorner)[0],:] = 0.0
        vecNormal[indBndCrv_Ext,:] = 0.0
        vecNormal[isValGT5Ext,:] = 0.0
        
        # Need to update vecTangent for the nodes on boundary curves
        p = n[indBndCrv]
        pplus = n[eB[indBndCrv,0]].squeeze()
        pminus = n[eB[indBndCrv,1]].squeeze()
        
        tminus = (.1/np.sqrt(np.sum((p-pminus)**2,axis=1))).reshape((p.shape[0],1))
        tplus =  (.1/np.sqrt(np.sum((p-pplus)**2,axis=1))).reshape((p.shape[0],1))
        
        tminus[np.nonzero(np.isnan(tminus))]=0.0
        tplus[np.nonzero(np.isnan(tplus))]=0.0
        
        crvTangent = tplus*pplus + (1.0-tplus)*p - ((1.0-tminus)*p+tminus*pminus)
        normalize_v3(crvTangent)
        indBndCrv = np.nonzero(np.any((-(eBExt==-1)),axis=1))
        vecTangent[indBndCrv,:] = ((vec[indBndCrv]*crvTangent).sum(1)).reshape((np.size(crvTangent,0),1))*crvTangent
        
        vecTangent[np.nonzero(isExtCorner)[0],:] = 0.0

        
        meanVecTangent[k_itt] = (np.sqrt( vecTangent[:,0]**2 + vecTangent[:,1]**2 + vecTangent[:,2]**2 )).mean()

        print("iteration %s: \t %.5f" % (str(k_itt+1), meanVecTangent[k_itt][0]))
        
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # dampen the motion for the first 10 or so iterations, which will 
        # help prevent oscillations ( 5 is arbitrary damping factor chosen)
        # done this way, there will be no motion on the first iteration, 
        # which is a cheap way to only update interior nodes if function 
        # called with itt = 1
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        vecTangent *= (1.0 - np.exp((k_itt) / -5.0))
        vecNormal   *=  (1.0 - np.exp((k_itt) / -5.0))
        
        # dampen the motion if the magnitude of the displacement is too large
        
        # find the length of the sides around each exterior node by 
        # calculating a "displacement"
        disp = np.zeros([np.size(nExt,0), 3])
        disp -= nExt
        disp = disp[:,np.newaxis,:]*np.ones([2*maxValence])[:,np.newaxis]
        disp += nExt[inverse_orderedOneNeigh]
        
        #calculate the characteristic length - sum(-1) is over x,y,z 
        # coordinates, then find max coordinate
        # replace zeros with nan and then use nanmin - BEWARE OF DIFFERENT 
        # NUMPYS POSSIBLY FAILING TO FIND 0.0 ALLOCATED FROM NP.ZEROS()
        char_length = np.sqrt((disp*disp).sum(-1))
        char_length[char_length == 0.0] = np.nan
        char_length = np.nanmin(char_length, axis = 1)
        
        # find the indices of vector magnitudes which are too large:
        # use 2.5% of the characteristic length as an arbitrary parameter
        vecTangTooBig = np.sqrt((vecTangent*vecTangent).sum(-1)) > char_length*0.025 
        vecNormTooBig = np.sqrt((vecNormal*vecNormal).sum(-1)) > char_length*0.025
        
        if np.any(vecTangTooBig) or np.any(vecNormTooBig):
            #find values that are too large, and then shorten the vectors 
            # to have length 2.5% of the characteristic length
            vecTangTooBig = np.nonzero(vecTangTooBig)[0]
            vecNormTooBig = np.nonzero(vecNormTooBig)[0]
            normedVecTang = vecTangent / np.sqrt((vecTangent*vecTangent).sum(-1))[:,np.newaxis]*np.ones([3])[np.newaxis,:]
            normedVecNorm = vecNormal / np.sqrt((vecNormal*vecNormal).sum(-1))[:,np.newaxis]*np.ones([3])[np.newaxis,:]
            vecTangent[vecTangTooBig] = 0.025*char_length[vecTangTooBig][:,np.newaxis]*np.ones([3])[np.newaxis,:] * normedVecTang[vecTangTooBig]
            vecNormal[vecNormTooBig] = 0.025*char_length[vecNormTooBig][:,np.newaxis]*np.ones([3])[np.newaxis,:] * normedVecNorm[vecNormTooBig]

        # keyword argument to disable tangent motion for extraordinary 
        # nodes in the tangent direction, sloppy since they are 
        # calculated and zeroed out
        if immobExtraord:
            vecTangent[isValGT5Ext,:] = 0.0
                
        # Update the exterior nodes
        if tangentMotion:
            nExt += vecTangent
        if normalMotion:
            nExt -= vecNormal
        
        # Store the exterior nodes
        n[np.nonzero(nIsOnExt)[0]] = nExt
        
        if internalMotion:
            cP = (n[e].sum(axis=1))/8.0
            eVols = hex_volume(e,n)
            
            mv = np.zeros((np.size(elUseVInt,0),3),'float32')   
            w_vols = eVols[elUseVInt]
            totVols = (eVols[elUseVInt].sum(1)).reshape((np.size(elUseVInt,0),1))
            totVols[np.nonzero(totVols==0)] = 1
            mv[:,0] = (cP[elUseVInt,0]*w_vols).sum(1)
            mv[:,1] = (cP[elUseVInt,1]*w_vols).sum(1) 
            mv[:,2] = (cP[elUseVInt,2]*w_vols).sum(1) 
            mv /= totVols
            
            nInt[0:-1] = mv
            n[np.nonzero(-nIsOnExt)[0]] = nInt
    
    eVols = hex_volume(e,n)
    if any(eVols < 0.0):
        print("some element volumes were negative")
        if len(np.nonzero(eVols < 0.0)[0]) < 500:
            print("negative elements:\n", np.nonzero(eVols < 0.0)[0].tolist())
        else:
            print("length of element volumes > 500, suppressed!")
            
    n = n[0:num_n,:]
    return n

def normalize_v3(arr):
    ''' Normalize a np array of 3 component vectors shape=(n,3) '''
    lens = np.sqrt( arr[:,0]**2 + arr[:,1]**2 + arr[:,2]**2 )
    # Correct for division by zero    
    lens[np.nonzero(lens==0)]=1
    arr[:,0] /= lens
    arr[:,1] /= lens
    arr[:,2] /= lens                
    return arr

def vert_normal(faces,vertices):
    # Need to create trifaces if faces are quads here two sets of trifaces are created
    # This will help decrease the influence of the choice of triangularion
    if np.size(faces,1)==4:
        faces = np.concatenate((faces[:,[0,1,2]],faces[:,[1,2,3]],faces[:,[2,3,0]],faces[:,[3,0,1]]))
    
    norm = np.zeros( vertices.shape, dtype=vertices.dtype )
    
    tris = vertices[faces]
    n = np.cross( tris[::,1 ] - tris[::,0]  , tris[::,2 ] - tris[::,0] )
    normalize_v3(n)
    
    norm[ faces[:,0] ] += n
    norm[ faces[:,1] ] += n
    norm[ faces[:,2] ] += n
    normalize_v3(norm)
    return norm

def hex_volume(e,n):
    # Compute the volume of hexahedrons by creating 24 tet elements using the center point
    cP = (n[e].sum(axis=1))/8.0
    vol = np.zeros((np.size(e,0)),'float32')
    
    faces = np.array([[0,3,2,1],[5,6,7,4],[0,1,5,4],[3,7,6,2],[0,4,7,3],[1,2,6,5]])
    for k in range(6):
        curface = e[:,faces[k]]
        fP = (n[curface].sum(axis=1))/4.0
        for j in range(4):
            a = n[curface[:,j]]
            b = n[curface[:,np.mod(j+1,4)]]
            c = fP
            d = cP
            vol += ((a-d)*np.cross(b-d,c-d)).sum(1)/6.0
    return vol

def face_area(f,n):
    # Compute the area of a quad face by subdividing it into 4 triangles
    fP = (n[f].sum(1))/4.
    a = np.zeros((np.size(f,0)),'float32')
    for i in range(4):
        a += 0.5*np.sqrt((np.cross(fP-n[f[:,i]],fP-n[f[:,np.mod(i+1,4)]])**2).sum(1))
    return a

def extraordinary_ext_edges(e,n):
    if e.min()==1:
        e = e-1
    
    num_e = np.size(e,0)
    num_n = np.size(n,0)
    
    # Find the exterior faces
    sf = np.array([[0,3,2,1],[5,6,7,4],[0,1,5,4],[3,7,6,2],[0,4,7,3],[1,2,6,5]])
    fM = {}
    for k in range(num_e):
        for i in range(6):
            face = [e[k][sf[i][0]],e[k][sf[i][1]],e[k][sf[i][2]],e[k][sf[i][3]]]
            face.sort()
            fInd = tuple(face)
            if fInd in fM:
                fM[fInd] = np.concatenate((fM[fInd],np.array([[k,i]])))
            else:
                fM[fInd] = np.array([[k,i]])
    
    faceKeys = list(fM.keys())
    faceInd = list(fM.values())
    nIsOnExt = [0]*num_n
    
    fExt = []
    fEndoEpi = []
    for i in range(len(faceKeys)):
        if np.size(faceInd[i],0)==1:
            nIsOnExt[faceKeys[i][0]] = 1
            nIsOnExt[faceKeys[i][1]] = 1
            nIsOnExt[faceKeys[i][2]] = 1
            nIsOnExt[faceKeys[i][3]] = 1
            fExt.append(e[faceInd[i][0][0],sf[faceInd[i][0][1]]])
            if faceInd[i][0][1] == 0 or faceInd[i][0][1] == 1:
                fEndoEpi.append(e[faceInd[i][0][0],sf[faceInd[i][0][1]]])
    
    fExt = np.array(fExt)

    # Find the number of connected exterior faces for each EXTERIOR edge
    sef = np.array([[0,1],[1,2],[2,3],[3,0]])
    eeM = {} # Exterior edge face Matrix
    for k in range(np.size(fExt,0)):
        for i in range(4):
            eInd = (min(fExt[k,sef[i,0]],fExt[k,sef[i,1]]),max(fExt[k,sef[i,0]],fExt[k,sef[i,1]]))
            if eInd in eeM:
                eeM[eInd].append(k)
            else:
                eeM[eInd] = [k]
    
    
    # Find the boundary curves
    # Find the number of elements for each EXTERIOR edge
    eM = {} # Edge Matrix
    
    se = np.array([[0,1],[1,2],[2,3],[3,0],[4,5],[5,6],[6,7],[7,4],[0,4],[1,5],[2,6],[3,7]])
    for k in range(len(e)):
        for i in range(12):
            eInd = (min(e[k][se[i,0]],e[k][se[i,1]]),max(e[k][se[i,0]],e[k][se[i,1]]))
            if eInd in eeM:
                if eInd in eM:
                    eM[eInd].append(k)
                else:
                    eM[eInd] = [k]
    
    
    edges = list(eM.keys())
    elUse_edge = list(eM.values())
    notUsedTwice = [1]*len(edges)
    edgeElCount = [-1]*len(edges)
    for i in range(len(edges)):
        edgeElCount[i] = np.size(elUse_edge[i])
        if edgeElCount[i]==2:
            notUsedTwice[i] = 0
    
    edgeElCount = np.array(edgeElCount)
    
    edges = np.array(edges)
    notUsedTwice = np.array(notUsedTwice)
    notUsedTwice = np.nonzero(notUsedTwice)[0]
    
    irr_ext_edge = edges[notUsedTwice,:]
    irr_ext_edge_valance = edgeElCount[notUsedTwice]
    
    
    # Check the number of connected exterior faces to the ext edge?
    # eeM would have k repeated - Not done since unlikely / bad practice to have
    # 'double-diamond' shape
    
    
    # Mark the boundary nodes and store their neighboring boundary nodes
    eB = np.ones((num_n,3,irr_ext_edge_valance.max()),int)
    eB = -eB
    for j in range(irr_ext_edge_valance.max()):
        for k in range(np.size(irr_ext_edge ,0)):
            if irr_ext_edge_valance[k] == j+1:
                if eB[irr_ext_edge[k,0],0,j]==-1:
                    eB[irr_ext_edge[k,0],0,j] = irr_ext_edge[k,1]
                elif eB[irr_ext_edge[k,0],1,j]==-1:
                    eB[irr_ext_edge[k,0],1,j] = irr_ext_edge[k,1]
                else:
                    eB[irr_ext_edge[k,0],2,j] = 1
                
                if eB[irr_ext_edge[k,1],0,j]==-1:
                    eB[irr_ext_edge[k,1],0,j] = irr_ext_edge[k,0]
                elif eB[irr_ext_edge[k,1],1,j]==-1:
                    eB[irr_ext_edge[k,1],1,j] = irr_ext_edge[k,0]
                else:
                    eB[irr_ext_edge[k,1],2,j] = 1
    
    # Find Corners
    isCorner = np.zeros((num_n,1),bool)
    for j in range(irr_ext_edge_valance.max()):
        for k in range(num_n):
            if (eB[k,2,j]==1):
                isCorner[k] = 1
    
    
    eBNew = np.ones((num_n,2),int)
    eBNew = -eBNew
    for j in range(irr_ext_edge_valance.max()):
        for k in range(num_n):
            if (eB[k,0,j]>-1):
                eBNew[k,:]  = eB[k,0:2,j]
    
    
     # Use the exterior faces
    nIsOnExt = np.array(nIsOnExt,'bool')
    
    return [fExt,nIsOnExt,eBNew,isCorner,eB]
    
def calc_vec_normal(n,
                    nExt,
                    normal,
                    totAreas,
                    orderedOneNeigh,
                    maxValence,
                    rotFaceVertMap,
                    inverse_orderedOneNeigh_normals, 
                    preserveRidges=True):
    
    np.seterr(divide = 'ignore', invalid = 'ignore')
    
    gauss_pts_u = np.array([0.5-np.sqrt(3.0)/6.0, 0.5+np.sqrt(3.0)/6.0, 0.5-np.sqrt(3.0)/6.0, 0.5+np.sqrt(3.0)/6.0])
    gauss_pts_v = np.array([0.5-np.sqrt(3.0)/6.0, 0.5-np.sqrt(3.0)/6.0, 0.5+np.sqrt(3.0)/6.0, 0.5+np.sqrt(3.0)/6.0])
    
    #making two copies, (u,v) and (u2,v2) with different dimensions to improve code readability
        
    # (u,v) : [valence, gauss points]
    #(u2,v2) : [valence, gauss points, coords x-y-z]
    u = np.ones([maxValence])[:,np.newaxis] * gauss_pts_u[np.newaxis,:]
    v = np.ones([maxValence])[:,np.newaxis] * gauss_pts_v[np.newaxis,:]
    u2 = u[:,:,np.newaxis] * np.ones([3])[np.newaxis,:]
    v2 = v[:,:,np.newaxis] * np.ones([3])[np.newaxis,:]

    num_ext_n = np.size(nExt,0)
    weights = np.zeros([num_ext_n, 2*maxValence])
    displacements = np.zeros([num_ext_n, 3])
    vecNormal = np.zeros([num_ext_n,3])
    valences = np.zeros([num_ext_n],np.int32)
    
    for node in range(num_ext_n):
        
        # curr_coords := [valence, 4, 3] --> stores coordinates.  For each node, find ordered faces' ordered vertices' xyz components
        curr_coords = n[rotFaceVertMap[node,:,:]]
    
        #ordering is different than the paper for the second term because we order our vertices (p1,p2,p3,p4) i.e. CCW/CW vs. (p1,p2,p4,p3)
        tang_u = (1.0 - v2 ) * ((curr_coords[:,1,:]-curr_coords[:,0,:])[:,np.newaxis,:]*np.ones([4])[:,np.newaxis]) + \
                      (       v2 ) * ((curr_coords[:,2,:]-curr_coords[:,3,:])[:,np.newaxis,:]*np.ones([4])[:,np.newaxis])
        tang_v = (1.0 - u2 ) * ((curr_coords[:,3,:]-curr_coords[:,0,:])[:,np.newaxis,:]*np.ones([4])[:,np.newaxis]) + \
                      (       u2 ) * ((curr_coords[:,2,:]-curr_coords[:,1,:])[:,np.newaxis,:]*np.ones([4])[:,np.newaxis])
                      
        denom = np.sqrt( (tang_u * tang_u).sum(-1) * (tang_v * tang_v).sum(-1) - ( (tang_u * tang_v).sum(-1) * (tang_u * tang_v).sum(-1) ) )
            
        alpha21 = ( ( (1.0 - v ) * ( tang_v * ( (v2 - 1.0)*tang_v - (u2 - 1.0)*tang_u ) ).sum(-1) ) / denom ).sum(-1) / 4.0
        alpha43 = ( ( ( v )        * ( tang_v * ( (v2 - 1.0)*tang_v - (u2 - 1.0)*tang_u ) ).sum(-1) ) / denom ).sum(-1) / 4.0
        alpha31 = ( ( (1.0 - u ) * ( tang_u * ( (u2 - 1.0)*tang_u - (v2 - 1.0)*tang_v ) ).sum(-1) ) / denom ).sum(-1) / 4.0
        alpha42 = ( ( ( u )        * ( tang_u * ( (u2 - 1.0)*tang_u - (v2 - 1.0)*tang_v) ).sum(-1) ) / denom ).sum(-1) / 4.0
        
        alpha = alpha42 - alpha21
        beta = alpha31 - alpha43
        gamma = alpha43 + alpha42
        
        alpha = alpha[-np.isnan(alpha)]
        beta = beta[-np.isnan(beta)]
        gamma = gamma[-np.isnan(gamma)]
        
        #alpha/beta/gamma have length equal to Valence at this point
        
        #find the last beta (or alpha, gamma...) that is not a nan... this is the current node's valence
        curr_valence = max(np.nonzero(-np.isnan(beta))[0])+1
        valences[node] = curr_valence
        
        #collate the alpha, beta, and gamma into the weights - note that
        #in the paper there is an indexing ambiguity that is resolved by assuming out-of-bounds indices cycle back to beginning/end

        weights[node, 1:2*curr_valence : 2] = 2.0 * gamma / totAreas[node]
        weights[node, 0 ]                          = 2.0 * (alpha[0] + beta[-1] ) / totAreas[node]
        weights[node, 2:2*curr_valence : 2] = 2.0 * (alpha[1:] + beta[:-1]) / totAreas[node]
    
    #now evaluate equation 14
    displacements -= nExt
    displacements = displacements[:,np.newaxis,:]*np.ones([2*maxValence])[:,np.newaxis]
    displacements += nExt[orderedOneNeigh]
    
    char_length = np.sqrt((displacements*displacements).sum(-1))
    char_length = char_length.sum(-1) / (2*valences)
    char_length = char_length[:,np.newaxis]*np.ones([3])[np.newaxis,:]
    
    weights = weights[:,:,np.newaxis]*np.ones([3])[np.newaxis,:]
    meanCrvNrml = ( weights*displacements ).sum(1)

    #now evaluate equation 16
    for node in range(num_ext_n):
        
        #evaluate "term 1" := n(p_i) * n(p_k).T * H(p_k)
        if preserveRidges:
            n_pk = normal[inverse_orderedOneNeigh_normals[node,:]]
            h_pk = meanCrvNrml[inverse_orderedOneNeigh_normals[node,:]]
        else:
            n_pk = normal[orderedOneNeigh[node,:]]
            h_pk = meanCrvNrml[orderedOneNeigh[node,:]]
        
        n_pi = normal[node,np.newaxis,:]*np.ones([2*maxValence])[:,np.newaxis]
        h_pi = meanCrvNrml[node,np.newaxis,:]*np.ones([2*maxValence])[:,np.newaxis]
        
        n_outer = n_pi[:,:,np.newaxis] * n_pk[:,np.newaxis,:]

        vecNormal[node,:] = ( weights[node,:,:] * ( np.einsum('ijk,ik->ij', n_outer, h_pk) - h_pi ) ).sum(0)
        vecNormal[node,:] *= 2.0
        
    #maybe i should vectorize in the node dimension when i figure out how much memory it would require for big meshes...
    #~ vecNormal = (weights* \
        #~ (np.einsum('ijkl,ijl->ijk',expanded_normals[:,:,:,np.newaxis] * normal[orderedOneNeigh][:,:,np.newaxis,:], meanCrvNrml[orderedOneNeigh]) - \
                    #~ meanCrvNrml[:,np.newaxis,:]*np.ones([2*maxValence])[:,np.newaxis]) \
                      #~ ).sum(1)
    #~ vecNormal *= 2
    
    #step size determined by characteristic length
    vecNormal *= (0.25*char_length*char_length*char_length)
    
    return vecNormal
