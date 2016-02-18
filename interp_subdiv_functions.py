"""
Collection of functions called by InterpSubdiv
@author: Greg Sturgeon
gregorymsturgeon@hotmail.com
April 26, 2011
"""

import numpy

def rotate_face_align_right(face,edge):
    rotate = 1
    while numpy.logical_and(-numpy.logical_and(face[2]==edge[0],face[1]==edge[1]),rotate<=8):
        face = [face[3], face[0], face[1], face[2]]
        rotate += 1
        if rotate == 5:
            face = [face[1], face[0], face[3], face[2]]
    return face

def include_duplicate_point(connFsTmp,f,k):
    # Find if any point should be included twice in the 1 neighborhood stencil
    # by finding if any edges that don't include k are shared by two faces in connFs
    eDup = []
    pDup = -1
    edgesTmp = []
    se = numpy.array([[0,1],[1,2],[2,3],[3,0]])
    for j in range(len(connFsTmp)):
        for i in range(4):
            eInd = [min(f[connFsTmp[j],se[i,0]],f[connFsTmp[j],se[i,1]]),max(f[connFsTmp[j],se[i,0]],f[connFsTmp[j],se[i,1]])]
            edgesTmp.append(eInd)
    
    edgesTmp = numpy.array(edgesTmp)    
    # Find edges which don't include k    
    maskNotk = numpy.in1d(numpy.array(edgesTmp),numpy.array(k))
    maskNotk.resize(len(connFsTmp)*4,2)    
    edgesNotk = edgesTmp[numpy.nonzero(-maskNotk.any(1))[0],:]
    
    # Find the duplicate edges that don't include k
    eMd = {}    
    for i in range(numpy.size(edgesNotk,0)):
        eInd = (min(edgesNotk[i,0],edgesNotk[i,1]),max(edgesNotk[i,0],edgesNotk[i,1]))
        if eInd in eMd:
            eMd[eInd] = [eMd[eInd],i]
        else:
            eMd[eInd] = i
    
    edgesNotk = list(eMd.keys())
    edgeCount = list(eMd.values())
    usedMult = [0]*len(edgesNotk)
    for i in range(len(edgesNotk)):
        if numpy.size(edgeCount[i])>1:
            usedMult[i] = 1
    
    usedMult = numpy.array(usedMult)
    edgesNotk = numpy.array(edgesNotk)
    eDup = edgesNotk[numpy.nonzero(usedMult),:][0]
    
    if numpy.size(eDup,0) > 0:
        print("###\n %s %s %s %s" % k, eDup[0], shape(eDup), edgesTmp)
        print("###\n %s" % eDup[1])
        print("###\n %s" % minimum(k,eDup[0]))
        # print "###\n", array([min(k,eDup[0]), max(k,eDup[0])]), edgesTmp
        print("###\n %s %s" % (any(numpy.in1d(numpy.array([minimum(k,eDup[0]),maximum(k,eDup[0])]),edgesTmp))))
        if not any(numpy.in1d(numpy.array([minimum(k,eDup[0]),maximum(k,eDup[0])]),edgesTmp)):
            pDup = eDup[0]
        elif not any(numpy.in1d(numpy.array([minimum(k,eDup[1]),maximum(k,eDup[1])]),edgesTmp)):
            pDup = eDup[1]
    
    return pDup

def vertex_weights_betas_gammas(n,alpha0,alpha1,alpha2):
    theta = 2.0*numpy.pi/n
    # Compute the Betas
    B = numpy.zeros((2*n+1,1),float)
    Bb = numpy.zeros((2*n+1,1),float)
    Bb[0] = alpha1**2
    Bb[1] = (4.0/n)*alpha1*alpha2
    Bb[2] = (alpha2*alpha2) - (2*(1-(4.0/n))*alpha0*alpha2) + ((1-(4.0/n))*(alpha0*alpha0))
    
    for j in range(1, int(n)):
        Bb[(2*j)+1] = (2.0/n) * alpha1 * ((1+numpy.cos(j*theta) + numpy.sin(j*theta))*alpha2 + (1-numpy.cos(j*theta) -numpy.sin(j*theta))*alpha0)
        Bb[(2*j)+2] = (4.0/n) * alpha0*(alpha2+(alpha2-alpha0)*numpy.cos(j*theta))
    
    B[0] = Bb[0]
    B[1] = .5*(Bb[1]+Bb[3])
    B[2] = Bb[2]
    B[3] = B[1]
    
    for j in range(4,2*int(n)+1):
        B[j] = .5 * (Bb[j] + Bb[2*n+4-j])
    
    # Compute the Gammas
    G = numpy.zeros((2*n+1,1),float)
    G[0] = .75
    G[1] = .25 + (1.0/(2*n))
    for j in range(1,int(n)):
        G[2*j] = 0
        G[(2*j)+1] = (1.0/(2*n))*numpy.cos(j*theta)
    
    return[B,G]

def circ_shift(a,n):
    if n<0:
        b = numpy.concatenate((a[abs(n):len(a)],a[0:abs(n)]),axis=0)
    else:
        b = numpy.concatenate((a[len(a)-n:len(a)],a[0:len(a)-n]),axis=0)
    
    return b
