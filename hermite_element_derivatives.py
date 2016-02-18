# -*- coding: utf-8 -*-
"""
HermiteElementDerivatives(eNN,nNN)

Computes the nodal derivatives for tricubic Hermites by element from a twice 
subdivided hexahedral mesh.

Input:
    eNN - (num_base_elements*64,8) Elements of a twice subdivided mesh as 
        produced from:
            [eN,nN] = hexInterpSubDiv(e,n)
            [eNN,nNN] = hexInterpSubDiv(eN,nN)     
    nNN - (num_nNN,3) Nodes of the twice subdivided mesh

Note that ideally the mesh should be regularized prior to computing the nodal 
derivatives.  This is due to the linear subdivision on the interior of the 
elements in the hexInterpSubDiv() function.
    nNN = regularizeElements(eNN,nNN,itt)

Output bxe, bye, bze:
    (8,8,num_base_elements) - nodal derivitives for each element

    ['000', '000u', '000w', '000t', '000uw', '000ut', '000wt', '000uwt']
    ['100', '100u', '100w', '100t', '100uw', '100ut', '100wt', '100uwt']
    ['110', '110u', '110w', '110t', '110uw', '110ut', '110wt', '110uwt']
    ['010', '010u', '010w', '010t', '010uw', '010ut', '010wt', '010uwt']
    ['001', '001u', '001w', '001t', '001uw', '001ut', '001wt', '001uwt']
    ['101', '101u', '101w', '101t', '101uw', '101ut', '101wt', '101uwt']
    ['111', '111u', '111w', '111t', '111uw', '111ut', '111wt', '111uwt']
    ['011', '011u', '011w', '011t', '011uw', '011ut', '011wt', '011uwt']


@author: Greg Sturgeon
gregorymsturgeon@hotmail.com
May 17, 2011
"""

import numpy

def hermite_element_derivatives(eNN,nNN):

    num_e = int(numpy.size(eNN,0)/64)
    
    # array of element node order such that
    # [eno0[0,0,4],eno1[0,0,4]] is the index into eNN for the node at eta3=0,eta1=0,eta2=1
    eno0 = numpy.array([[[0,1,8,9,9],[3,2,11,10,10],[24,25,16,17,17],[27,26,19,18,18],[27,26,19,18,18]],
    [[4,5,12,13,13],[7,6,15,14,14],[28,29,20,21,21],[31,30,23,22,22],[31,30,23,22,22]],
    [[32,33,40,41,41],[35,34,43,42,42],[56,57,48,49,49],[59,58,51,50,50],[59,58,51,50,50]],
    [[36,37,44,45,45],[39,38,47,46,46],[60,61,52,53,53],[63,62,55,54,54],[63,62,55,54,54]],
    [[36,37,44,45,45],[39,38,47,46,46],[60,61,52,53,53],[63,62,55,54,54],[63,62,55,54,54]]])
    
    eno1 = numpy.array([[[0,0,0,0,1],[0,0,0,0,1],[0,0,0,0,1],[0,0,0,0,1],[3,3,3,3,2]],
    [[0,0,0,0,1],[0,0,0,0,1],[0,0,0,0,1],[0,0,0,0,1],[3,3,3,3,2]],
    [[0,0,0,0,1],[0,0,0,0,1],[0,0,0,0,1],[0,0,0,0,1],[3,3,3,3,2]],
    [[0,0,0,0,1],[0,0,0,0,1],[0,0,0,0,1],[0,0,0,0,1],[3,3,3,3,2]],
    [[4,4,4,4,5],[4,4,4,4,5],[4,4,4,4,5],[4,4,4,4,5],[7,7,7,7,6]]])
    
    eno0ind = eno0.copy()
    eno1ind = eno1.copy()
    for k in range(5):
        eno0ind[:,:,k] = eno0[k,:,:].copy()
        eno1ind[:,:,k] = eno1[k,:,:].copy()
    
    # ind contains the indicies to the nodes in each element ordered:
    # ind[elementnumber,eta1,eta2,eta3] 
    ind125 = numpy.zeros((num_e,5,5,5),int)
    uwt125 = numpy.zeros((num_e,5,5,5,3),float)
    ind = numpy.zeros((num_e,4,4,4),int)
    uwt = numpy.zeros((num_e,4,4,4,3),float)
    
    # c will contain a row of the 64 u,w,t products for each of the 64 points
    c = numpy.zeros((64,64),float)
    xc = numpy.zeros((64),float)
    yc = numpy.zeros((64),float)
    zc = numpy.zeros((64),float)
    
    # array to rearrange the 64 element column vectors into 8x8 array ordered as described in the comments 
    re = numpy.array([[0,32,8,2,40,34,10,42],[16,48,24,18,56,50,26,58],[20,52,28,22,60,54,30,62],[4,36,12,6,44,38,14,46],[1,33,9,3,41,35,11,43],[17,49,25,19,57,51,27,59],[21,53,29,23,61,55,31,63],[5,37,13,7,45,39,15,47]])
    bxe = numpy.zeros((8,8,num_e),float)
    bye = numpy.zeros((8,8,num_e),float)
    bze = numpy.zeros((8,8,num_e),float)
    
    # create an array to map from the 4x4x4 array to a 64 element column vector
    i64 = numpy.array(list(range(0,64)))
    i64 = i64.reshape(4,4,4)
    
    for k_e in range(num_e):
        ind125[k_e] = eNN[k_e*64+eno0ind,eno1ind]
        # Determine the lengths
        lu = numpy.sqrt(((nNN[ind125[k_e,:,0:4,:]] - nNN[ind125[k_e,:,1:5,:]])**2).sum(3))
        lw = numpy.sqrt(((nNN[ind125[k_e,0:4,:,:]] - nNN[ind125[k_e,1:5,:,:]])**2).sum(3))
        lt = numpy.sqrt(((nNN[ind125[k_e,:,:,0:4]] - nNN[ind125[k_e,:,:,1:5]])**2).sum(3))
        
        lu = numpy.cumsum(lu,axis=1)
        lw = numpy.cumsum(lw,axis=0)
        lt = numpy.cumsum(lt,axis=2)
        
        u = numpy.zeros((5,5,5,1),float)
        w = numpy.zeros((5,5,5,1),float)
        t = numpy.zeros((5,5,5,1),float)
        
        for k in range(4):
            u[:,k+1,:,0] = lu[:,k,:]/lu[:,3,:]
            w[k+1,:,:,0] = lw[k,:,:]/lw[3,:,:]
            t[:,:,k+1,0] = lt[:,:,k]/lt[:,:,3]
        
        uwt125[k_e,:,:,:,:] = numpy.concatenate((u,w,t),axis=3)
        
        # select 64 points and paramater values from the 125
        ind[k_e,:,:,:] = ind125[numpy.ix_([k_e],[0,1,3,4],[0,1,3,4],[0,1,3,4])].copy()
        uwt[k_e,:,:,:] = uwt125[numpy.ix_([k_e],[0,1,3,4],[0,1,3,4],[0,1,3,4],[0,1,2])].copy()
        x = nNN[ind[k_e,:,:,:],0]
        y = nNN[ind[k_e,:,:,:],1]
        z = nNN[ind[k_e,:,:,:],2]
        
        c = numpy.zeros((64,64),float)
        # create the u,w,t products for each of the 64 points.
        for i1 in range(0,4):
            for i2 in range(0,4):
                for i3 in range(0,4):
                    for j1 in range(1,5):
                        for j2 in range(1,5):
                            for j3 in range(1,5):
                                c[i64[i1,i2,i3],i64[j1-1,j2-1,j3-1]] = fH(j1,uwt[k_e,i1,i2,i3,0])*fH(j2,uwt[k_e,i1,i2,i3,1])*fH(j3,uwt[k_e,i1,i2,i3,2])
                                xc[i64[i1,i2,i3]] = x[i1,i2,i3]
                                yc[i64[i1,i2,i3]] = y[i1,i2,i3]
                                zc[i64[i1,i2,i3]] = z[i1,i2,i3]
        
        # Find the geometric coefficients bx,by,bz
        cinv = numpy.linalg.inv(c)
        bx = numpy.dot(cinv,xc)
        by = numpy.dot(cinv,yc)
        bz = numpy.dot(cinv,zc)
        
        for i in range(8):
            for j in range(8):
                bxe[i,j,k_e] = bx[re[i,j]]
                bye[i,j,k_e] = by[re[i,j]]
                bze[i,j,k_e ]= bz[re[i,j]]

    return bxe,bye,bze

def fH(j,u):
    if j==1:
        f = 2*u**3 - 3*u**2 + 1
    elif j==2:
        f = -2*u**3 + 3*u**2
    elif j==3:
        f = u**3 - 2*u**2 + u
    elif j==4:
        f = u**3 - u**2
    return f
