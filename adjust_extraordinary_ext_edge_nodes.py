# -*- coding: utf-8 -*-
"""
Adjust nodes on extraordinary exterior edges based on an approximating subdivision
scheme.  Only the positions of nodes on the extraordinary exterior edges are updated.

e = num_e x 8 array of the elements
n = num_n x 3 array of the vertex coordinates

@author: Greg Sturgeon
gregorymsturgeon@hotmail.com
October 31, 2011
"""

import numpy
from hexblender.regularize_elements import extraordinary_ext_edges

def adjust_extraordinary_ext_edge_nodes(e, n):
    
    if e.min()==1:
            e = e-1
    
    num_n = numpy.size(n,0)
    [fExt,nIsOnExt,eB,isCorner,eBval] = ExtraordinaryExtEdges(e,n)
    
    
    # Find nodes on val3 or higher ext edges
    isOnIrrExtEdge = numpy.zeros((num_n,1),bool)
    for k in range(2,numpy.size(eBval,2)):
        for i in range(num_n):
            if eBval[i,0,k]>-1:
                isOnIrrExtEdge[i] = 1
    
    
    # Find edge connected and face connected for fExt
    sef = numpy.array([[0,1],[1,2],[2,3],[3,0]])
    eConnV = [[-1]]*num_n
    for k in range(numpy.size(fExt,0)):
        for i in range(4):
            if isOnIrrExtEdge[fExt[k,sef[i,0]]]:
                if eConnV[fExt[k,sef[i,0]]][0] ==-1:
                    eConnV[fExt[k,sef[i,0]]] = [fExt[k,sef[i,1]]]
                else:
                    eConnV[fExt[k,sef[i,0]]].append(fExt[k,sef[i,1]])
            if isOnIrrExtEdge[fExt[k,sef[i,1]]]:
                if eConnV[fExt[k,sef[i,1]]][0] ==-1:
                    eConnV[fExt[k,sef[i,1]]] = [fExt[k,sef[i,0]]]
                else:
                    eConnV[fExt[k,sef[i,1]]].append(fExt[k,sef[i,0]])
    
    for k in range(num_n):
        if eConnV[k][0]>-1:
            eConnV[k] = numpy.unique(numpy.array(eConnV[k]))
    
    # fExt  connected to each node
    fExtUseV = [[-1]]*num_n
    for k in range(numpy.size(fExt,0)):
        for i in range(4):
            if fExtUseV[fExt[k,i]][0]==-1:
                fExtUseV[fExt[k,i]] = [k]
            else:
                fExtUseV[fExt[k,i]].append(k)
                
    # Face Points connected to each exterior node
    fP = (n[fExt].sum(1))/4
    
    
    # Calculate the adjusted position of nodes on the val3 and higher ext edges
    # based on the face points of adjacent faces and midpoints of adjacent edges
    nNew = n.copy()
    m = numpy.zeros((num_n,1),int)
    for k in range(num_n):
        if eConnV[k][0]>-1:
            m[k] = numpy.size(fExtUseV[k])
            nNew[k,:] = ((n[eConnV[k],:].sum(0))/(2*numpy.size(eConnV[k])) + (2*(fP[fExtUseV[k],:].sum(0)))/m[k] + (m[k]-3)*n[k,:] + (n[k,:]/2))/m[k]
    
    return nNew
