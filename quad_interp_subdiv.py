import numpy
from math import pi, cos, sin
from hexblender.set_xor_2d import set_xor_2d
from functools import cmp_to_key

def quad_interp_subdiv(f,v):
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % Perform Interpolatory subdivision surface on a quad mesh
# % Based on G.Li, W. Ma and H. Bao, A New Interpolatory subdivision for
# % Quadrilateral Meshes, Computer Graphics forum, Vol 24 (2005) 1 pp3-16.
# % Similar to Kobbelt interpolatory subdivision scheme.
# %
# % [fN vN] = interpSubDiv(f,v)
# % [fN vN] = interpSubDiv(f,v,'pre-cleaned')
# %
# % f = num_f x 4 numpy.array of the faces
# % v = num_v x 3 numpy.array of the vertex coordinates
# % optional argument 'pre-cleaned' specifies that the replaceSharedVertices
# % and deleteSharedFace functions have already been applied to the mesh
# %
# % fN faces of subdivided surface
# % vN vertices of sudvivided surface
# %
# %
# % Revision 1.0
# % Initial Release June 2010
# %
# % Revision 1.1
# % Clean mesh to remove shared vertices and duplicate faces
# % Account for thin strips of quad meshes (with all vertices on the mesh
# % boundary)
# % Remove call to external function 'edgesQuadFace'
# % Check if vertex ismember of stencil prior to adding it
# %
# % September 2010
# % Greg Sturgeon gregorymsturgeon@hotmai.com
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    # Weighting values
    alpha0 = -1/8.
    alpha1 = 6/8.
    alpha2 = 3/8.

    # Number of vertices and faces in Base Mesh
    num_v = numpy.shape(v)[0]
    num_f = numpy.shape(f)[0]

    # Mark an edge (vertex pair) in edgeMatrix1 the first time an edge is used
    # and edgeMatrix2 the second time.  The edgeMatrix contains the face # that
    # uses the edge.  These edges are ordered by node number.

    # eM1 = numpy.zeros((num_v,num_v))
    # eM2 = numpy.zeros((num_v,num_v))
    eM1 = Matrix(num_v,num_v)
    eM2 = Matrix(num_v,num_v)

    # eConnV = cell(num_v,1);
    eConnV = list(range(num_v))
    for i in range(num_v):
        eConnV[i] = []
    # connFs = cell(num_v,1);
    connFs = list(range(num_v))
    for i in range(num_v):
        connFs[i] = []



    for i in range(num_f):

        if not eM1[numpy.minimum(f[i,0],f[i,1]),numpy.maximum(f[i,0],f[i,1])]:
            eM1[numpy.minimum(f[i,0],f[i,1]),numpy.maximum(f[i,0],f[i,1])]=i+1
        else:
            eM2[numpy.minimum(f[i,0],f[i,1]),numpy.maximum(f[i,0],f[i,1])]=i+1

        if not eM1[numpy.minimum(f[i,1],f[i,2]),numpy.maximum(f[i,1],f[i,2])]:
            eM1[numpy.minimum(f[i,1],f[i,2]),numpy.maximum(f[i,1],f[i,2])]=i+1
        else:
            eM2[numpy.minimum(f[i,1],f[i,2]),numpy.maximum(f[i,1],f[i,2])]=i+1

        if not eM1[numpy.minimum(f[i,2],f[i,3]),numpy.maximum(f[i,2],f[i,3])]:
            eM1[numpy.minimum(f[i,2],f[i,3]),numpy.maximum(f[i,2],f[i,3])]=i+1
        else:
            eM2[numpy.minimum(f[i,2],f[i,3]),numpy.maximum(f[i,2],f[i,3])]=i+1

        if not eM1[numpy.minimum(f[i,3],f[i,0]),numpy.maximum(f[i,3],f[i,0])]:
            eM1[numpy.minimum(f[i,3],f[i,0]),numpy.maximum(f[i,3],f[i,0])]=i+1
        else:
            eM2[numpy.minimum(f[i,3],f[i,0]),numpy.maximum(f[i,3],f[i,0])]=i+1

        # eConnV is a cell numpy.array storing all the vertices connected to a
        # given vertex.  eConnV{k} has an numpy.array of all the vertices that are
        # connected to vertex k.

        if f[i,3] not in eConnV[f[i,0]]:
            eConnV[f[i,0]].append(f[i,3])

        if f[i,1] not in eConnV[f[i,0]]:
            eConnV[f[i,0]].append(f[i,1])

        if f[i,0] not in eConnV[f[i,1]]:
            eConnV[f[i,1]].append(f[i,0])

        if f[i,2] not in eConnV[f[i,1]]:
            eConnV[f[i,1]].append(f[i,2])

        if f[i,1] not in eConnV[f[i,2]]:
            eConnV[f[i,2]].append(f[i,1])

        if f[i,3] not in eConnV[f[i,2]]:
            eConnV[f[i,2]].append(f[i,3])

        if f[i,2] not in eConnV[f[i,3]]:
            eConnV[f[i,3]].append(f[i,2])

        if f[i,0] not in eConnV[f[i,3]]:
            eConnV[f[i,3]].append(f[i,0])

        # connFs is a cell numpy.array storing all the faces connected to a given
        # vertex
        if i not in connFs[f[i,0]]:
           connFs[f[i,0]].append(i)

        if i not in connFs[f[i,1]]:
           connFs[f[i,1]].append(i)

        if i not in connFs[f[i,2]]:
           connFs[f[i,2]].append(i)

        if i not in connFs[f[i,3]]:
           connFs[f[i,3]].append(i)

    # Find the boundary edges those used by only one face
    inds = list(eM1.keys())
    inds.sort()
    i1 = numpy.fromiter((ind[0] for ind in inds), dtype='d')
    j1 = numpy.fromiter((ind[1] for ind in inds), dtype='d')
    inds = list(eM2.keys())
    inds.sort()
    i2 = numpy.fromiter((ind[0] for ind in inds), dtype='d')
    j2 = numpy.fromiter((ind[1] for ind in inds), dtype='d')

    a = numpy.vstack((i1,j1))
    b = numpy.vstack((i2,j2))

    a = a.transpose()
    b = b.transpose()
    e_boundary = set_xor_2d(a, b)
    e_boundary = numpy.array(e_boundary)

    b = numpy.zeros((num_v,1))
    b[list(e_boundary[:,0])]=numpy.nan
    b[list(e_boundary[:,1])]=numpy.nan

    # Allocate memory for qE and qF
    qE1_x = numpy.zeros((num_v,num_v))
    qE1_y = numpy.zeros((num_v,num_v))
    qE1_z = numpy.zeros((num_v,num_v))

    qE2_x = numpy.zeros((num_v,num_v))
    qE2_y = numpy.zeros((num_v,num_v))
    qE2_z = numpy.zeros((num_v,num_v))

    #qF = cell(num_f,1); in MatLab
    qF = list(range(num_f))
    for i in range(num_f):
        qF[i] = []

    eConnV = numpy.array(eConnV)
    connFs = numpy.array(connFs)

    for k in range(num_v):
        # Only operating on non-bounday vertices
        if not all(numpy.isnan(b[k])):

            # fConn is a face [v1 v2 v3 v4] connected to vertex k
            fNext = connFs[k][0]
            fConn = f[connFs[k][0],:]

            # number of connected edges
            n=len(eConnV[k])    # n=numpy.shape(eConnV[k])[1]

            # p stores the vertex indices in the 1 neighborhood stencil (see Figure 4).
            p=numpy.zeros((2*n,1))
            ip=0 # index into p
            p[0] = (eConnV[k][0])+1

            fCCW = [] # Store the faces touching vertex k in a CCW direction.
            fCCW = numpy.array(fCCW)

            while any([not val for val in p]):

                if fNext not in fCCW:
                    fCCW = numpy.hstack((fCCW,fNext))

                # Rotate fConn so that edge k p(ip) is the Right edge
                edge = (p[ip][0])
                edge = int(edge)-1
                fConn = rotate_face_align_right(fConn,[k,edge])

                # Store the vertex indices in the 1 neighborhood stencil.
                if fConn[1]+1 not in p:
                    ip=ip+1
                    p[ip] = fConn[1]+1

                if fConn[0]+1 not in p:
                    ip=ip+1
                    p[ip] = (fConn[0])+1

                if fConn[3]+1 not in p:
                    ip=ip+1
                    p[ip] = fConn[3]+1

                if  (eM1[numpy.minimum(int(k),int(p[ip]-1)),numpy.maximum(int(k),int(p[ip]-1))]-1) not in fCCW and eM1[numpy.minimum(int(k),int(p[ip]-1)), numpy.maximum(int(k),int(p[ip]-1))]:
                    fNext = eM1[numpy.minimum(int(k),int(p[ip]-1)), numpy.maximum(int(k),int(p[ip]-1))]-1
                    fConn = f[fNext,:]

                elif (eM2[numpy.minimum(int(k),int(p[ip]-1)), numpy.maximum(int(k),int(p[ip]-1))]-1) not in fCCW and eM2[numpy.minimum(int(k),int(p[ip]-1)), numpy.maximum(int(k),int(p[ip]-1))]:
                    fNext = eM2[numpy.minimum(int(k),int(p[ip]-1)), numpy.maximum(int(k),int(p[ip]-1))]-1
                    fConn = f[fNext,:]

                elif any([not val for val in p]):
                    print("Check to verify that the mesh is a manifold surface")
                    return

            for i in range(len(p)):
                p[i] = p[i]-1

            B,G = vertex_weights_betas_gammas(n,alpha0,alpha1,alpha2)

            fCCW = fCCW.transpose().copy()

            # Rotate the order of p and fCCW to fill in all the qF's that touch
            # the vertex k.
            for j in range(numpy.shape(fCCW)[0]):
                # pCirc = circshift(p,-2*(j-1)); in MatLab
                pCirc = numpy.roll(p, (-2*(j)), axis = 0)                                                 # Stand 05/20/2011 Fehler auch in def rotate...()
                # fCCWCirc = circshift(fCCW,-(j-1)); in MatLab
                fCCWCirc = numpy.roll(fCCW,(-(j)), axis = 0)

                # pIncCent includes the center vertex k in the 1-neighborhood
                pIncCent = numpy.vstack((numpy.array(k),pCirc))
                pIncCent = list(pIncCent)
                val1 = list(numpy.sum(B[:]*v[pIncCent,0],axis=0))
                val1 = val1[0]
                val2 = list(numpy.sum(B[:]*v[pIncCent,1],axis=0))
                val2 = val2[0]
                val3 = list(numpy.sum(B[:]*v[pIncCent,2],axis=0))
                val3 = val3[0]
                coord = numpy.array([val1, val2, val3])

                if qF[int(fCCWCirc[1-1])] == []:
                    qF[int(fCCWCirc[1-1])] = coord
                else:
                    qF[int(fCCWCirc[1-1])] = numpy.vstack((qF[int(fCCWCirc[1-1])],coord))

                if not qE1_x[int(numpy.minimum(k,pCirc[0])),int(numpy.maximum(k,pCirc[0]))]:
                    qE1_x[int(numpy.minimum(k,pCirc[0])),int(numpy.maximum(k,pCirc[0]))] = numpy.sum(G[:]*v[pIncCent,0],axis=0)[0]
                    qE1_y[int(numpy.minimum(k,pCirc[0])),int(numpy.maximum(k,pCirc[0]))] = numpy.sum(G[:]*v[pIncCent,1],axis=0)[0]
                    qE1_z[int(numpy.minimum(k,pCirc[0])),int(numpy.maximum(k,pCirc[0]))] = numpy.sum(G[:]*v[pIncCent,2],axis=0)[0]
                else:
                    qE2_x[int(numpy.minimum(k,pCirc[0])),int(numpy.maximum(k,pCirc[0]))] = numpy.sum(G[:]*v[pIncCent,0],axis=0)[0]
                    qE2_y[int(numpy.minimum(k,pCirc[0])),int(numpy.maximum(k,pCirc[0]))] = numpy.sum(G[:]*v[pIncCent,1],axis=0)[0]
                    qE2_z[int(numpy.minimum(k,pCirc[0])),int(numpy.maximum(k,pCirc[0]))] = numpy.sum(G[:]*v[pIncCent,2],axis=0)[0]

    # Average the qFs for the faces that have at least one interior vertex.
    fPn = numpy.zeros((num_f,3))
    for k in range(num_f):
        if len(qF[k]) == 0:
            L = 0
        elif len(numpy.shape(qF[k])) == 1:
            L = 1
        else:
            L = numpy.shape(qF[k])[0]
        if L != 0:
            if len(numpy.shape(qF[k])) == 1:
                fPn[k,:] = [numpy.sum(qF[k][0],axis=0)/L, numpy.sum(qF[k][1],axis=0)/L, numpy.sum(qF[k][2],axis=0)/L]
            elif len(numpy.shape(qF[k])) == 2:
                fPn[k,:] = [numpy.sum(qF[k][:,0],axis=0)/L, numpy.sum(qF[k][:,1],axis=0)/L, numpy.sum(qF[k][:,2],axis=0)/L]

        else:
            # Face has only exterior / boundary vertices.
            fPn[k,:] = [.25*(numpy.sum(v[f[k,:],0],axis=0)), .25*(numpy.sum(v[f[k,:],1],axis=0)), .25*(numpy.sum(v[f[k,:],2],axis=0))]

    inds = list(eM1.keys())

    def compare(x, y):
        if x[1] < y[1]:
            return -1
        elif (x[1] == y[1]) and (x[0] < y[0]):
            return -1
        else:
            return 0

    inds = sorted(inds, key = cmp_to_key(compare))
    # old way, not valid in Python 3
    #inds.sort(key=lambda x, y: -1 if x[1]<y[1] else ( -1 if (x[1] == y[1]) and (x[0] < y[0]) else 1))
    e1 = numpy.fromiter((ind[0] for ind in inds), dtype='d')
    e2 = numpy.fromiter((ind[1] for ind in inds), dtype='d')

    edgeArray = [numpy.minimum(e1,e2), numpy.maximum(e1,e2)]
    edgeArray = numpy.array(edgeArray)
    edgeArray = edgeArray.transpose()

    ismember = numpy.zeros(numpy.shape(edgeArray))
    for i in range(numpy.shape(edgeArray)[0]):
        for j in range(numpy.shape(edgeArray)[1]):
            if edgeArray[i,j] in e_boundary:
                ismember[i,j] = numpy.array([1.])

    numIntV = numpy.array([2 - numpy.sum(ismember,axis=1)])
    numIntV = numIntV.transpose()

    eWithNoIntV = []
    for i,Q in enumerate(numIntV):
        if Q == 0:
            eWithNoIntV.extend(numpy.array([edgeArray[i,:]]))
    eWithIntV = []
    L = []
    for i,Q in enumerate(numIntV):
        if Q > 0:
            eWithIntV.extend(numpy.array([edgeArray[i,:]]))
            L.extend(Q)

    eWithNoIntV = numpy.array(eWithNoIntV)
    eWithIntV = numpy.array(eWithIntV)

    intEwithTwoBndV = numpy.setdiff1d(eWithNoIntV, e_boundary)
    intEwithTwoBndV = numpy.array([intEwithTwoBndV])

    # Average the qEs for the edges that have at least one interior vertex
    eWithIntV = eWithIntV.astype(int)
    
    qE1_x = numpy.array(qE1_x[eWithIntV[:,0],eWithIntV[:,1]]).transpose()
    qE1_y = numpy.array(qE1_y[eWithIntV[:,0],eWithIntV[:,1]]).transpose()
    qE1_z = numpy.array(qE1_z[eWithIntV[:,0],eWithIntV[:,1]]).transpose()
    qE2_x = numpy.array(qE2_x[eWithIntV[:,0],eWithIntV[:,1]]).transpose()
    qE2_y = numpy.array(qE2_y[eWithIntV[:,0],eWithIntV[:,1]]).transpose()
    qE2_z = numpy.array(qE2_z[eWithIntV[:,0],eWithIntV[:,1]]).transpose()

    # inds = qE1_x.keys()
    # inds.sort(lambda x, y: -1 if x[1]<y[1] else ( -1 if (x[1] == y[1]) and (x[0] < y[0]) else 1))
    # val_qE1_x = numpy.array([numpy.fromiter((qE1_x[ind] for ind in inds), dtype='d')]).transpose()
        
    QE_x = (qE1_x+qE2_x)
    QE_y = qE1_y+qE2_y
    QE_z = qE1_z+qE2_z
    L = numpy.array(L)
    ePn = numpy.array([QE_x/L, QE_y/L, QE_z/L]).transpose()

    # New edge points for boundary edges.
    eB = numpy.zeros((num_v,2))
    # Store the boundary vertices that touch the kth vertex
    for k in range(numpy.shape(e_boundary)[0]):
        if  not eB[e_boundary[k,0],0]:
            eB[e_boundary[k,0],0] = e_boundary[k,1]
        else:
            eB[e_boundary[k,0],1] = e_boundary[k,1]

        if not eB[e_boundary[k,1],0]:
            eB[e_boundary[k,1],0] = e_boundary[k,0]
        else:
            eB[e_boundary[k,1],1] = e_boundary[k,0]


    ebPn = numpy.zeros((numpy.shape(e_boundary)[0],3))
    # Weighting and averaging of the boundary vertices
    for k in range(numpy.shape(e_boundary)[0]):
        p0 = e_boundary[k,0]
        p1 = e_boundary[k,1]
        pminus1 = list(numpy.setdiff1d(eB[p0,:],p1))
        p2 = list(numpy.setdiff1d(eB[p1,:],p0))

        qbE0 = alpha0*v[pminus1,:] + alpha1*v[p0,:] + alpha2*v[p1,:]
        qbE1 = alpha2*v[p0,:] + alpha1*v[p1,:] + alpha0*v[p2,:]

        ebPn[k,:] =  .5 *  (qbE0 +qbE1)

    # New edge points for intEwithTwoBndV
    if numpy.shape(intEwithTwoBndV)[1] >= 2:
        eintPn = .5*(v[list(intEwithTwoBndV[:,0]),:]+v[list(intEwithTwoBndV[:,1]),:])
    else:
        eintPn = numpy.array([])

    # Construct the topology for the new faces
    # Store where the edge point for each edge is in vN
    eM = numpy.zeros((num_v,num_v))

    for i in range(numpy.shape(eWithIntV)[0]):
        eM[eWithIntV[i,0],eWithIntV[i,1]] = i + num_v + num_f

    for i in range(numpy.shape(e_boundary)[0]):
        eM[e_boundary[i,0],e_boundary[i,1]] = i + num_v + num_f + numpy.shape(eWithIntV)[0]

    for i in range(numpy.shape(eintPn)[0]):
        eM[intEwithTwoBndV[i,0],intEwithTwoBndV[i,1]] = i + num_v + num_f + numpy.shape(eWithIntV)[0] + numpy.shape(e_boundary)[0]


    vN = v
    vN = numpy.vstack((vN,fPn))
    vN = numpy.vstack((vN,ePn))
    vN = numpy.vstack((vN,ebPn))
    if numpy.shape(eintPn) != (0,):
        vN = numpy.vstack((vN,eintPn))

    # Create the new face array
    f1 = numpy.zeros((num_f,4))
    f2 = numpy.zeros((num_f,4))
    f3 = numpy.zeros((num_f,4))
    f4 = numpy.zeros((num_f,4))
    fN = numpy.zeros((4*num_f,4))
    for k in range(num_f):
        f1[k,:] = [f[k,0], eM[numpy.minimum(f[k,0],f[k,1]),numpy.maximum(f[k,0],f[k,1])], num_v+k, eM[numpy.minimum(f[k,0],f[k,3]),numpy.maximum(f[k,0],f[k,3])]]
        f2[k,:] = [eM[numpy.minimum(f[k,0],f[k,1]),numpy.maximum(f[k,0],f[k,1])], f[k,1], eM[numpy.minimum(f[k,1],f[k,2]),numpy.maximum(f[k,1],f[k,2])], num_v+k]
        f3[k,:] = [num_v+k, eM[numpy.minimum(f[k,1],f[k,2]),numpy.maximum(f[k,1],f[k,2])], f[k,2], eM[numpy.minimum(f[k,2],f[k,3]),numpy.maximum(f[k,2],f[k,3])]]
        f4[k,:] = [eM[numpy.minimum(f[k,0],f[k,3]),numpy.maximum(f[k,0],f[k,3])], num_v+k, eM[numpy.minimum(f[k,2],f[k,3]),numpy.maximum(f[k,2],f[k,3])], f[k,3]]
        fN[4*(k),:] = f1[k,:]
        fN[4*(k)+1,:] = f2[k,:]
        fN[4*(k)+2,:] = f3[k,:]
        fN[4*(k)+3,:] = f4[k,:]

    fNlist = []
    for Array in fN:
        fNlist.append(list(Array))
    for i in range(len(fNlist)):
        for j in range(len(fNlist[i])):
            fNlist[i][j] = int(fNlist[i][j])

    return fNlist, vN

def rotate_face_align_right(face,edge):

    Qa = (not(face[2] == edge[0] and face[1] == edge[1]))
    Qb = face[2] == edge[0]
    Qc = face[1] == edge[1]


    rotate = 1
    while not(face[2] == edge[0] and face[1] == edge[1]) and rotate <=8:
        Qd = (not(face[2] == edge[0] and face[1] == edge[1]) and rotate <=8)
        face = [face[4-1], face[1-1], face[2-1], face[3-1]]
        rotate = rotate +1
        if rotate == 5:
            face = [face[1], face[0], face[3], face[2]]
        if rotate == 9:
            print("Could not align Faces.  Check that all the normals are outward facing.")

    return face

def vertex_weights_betas_gammas(n,alpha0,alpha1,alpha2):

    theta = 2*pi/n

    # Compute the Betas
    B = numpy.zeros((2*n+1,1))
    Bb = numpy.zeros((2*n+1,1))
    Bb[0] = alpha1**2
    Bb[1] = (4./n) *alpha1*alpha2
    Bb[2] = alpha2**2 - 2*(1-(4./n))*alpha0*alpha2 + (1-(4./n))*alpha0**2

    for j in range(1,n):
        Bb[(2*j+1)] = (2./n) * alpha1 * ((1+cos(j*theta) + sin(j*theta))*alpha2 + (1 - cos(j*theta) - sin(j*theta))*alpha0)
        Bb[(2*j+2)] = (4./n) * alpha0*(alpha2+(alpha2-alpha0)*cos(j*theta))

    B[0] = Bb[0]
    B[1] = (1/2.)*(Bb[1] + Bb[3])
    B[2] = Bb[2]
    B[3] = B[1]

    for j in range(3,(2*n)):
        B[(j+1)] = (1/2.) * (Bb[(j+1)]+Bb[(2*n+3-j)])

    # Compute the Gammas
    G = numpy.zeros((2*n+1,1))
    G[0] = .75
    G[1] = .25 + (1./(2.*n))
    for j in range(1,n):
       G[(2*j)] = 0
       G[(2*j+1)] = (1./(2.*n))*cos(j*theta)

    return B, G

class Matrix(dict):
    def __init__(self, m, n):
        self._m = m
        self._n = n
    def __getitem__(self, key):
        if not isinstance(key, tuple):
            raise TypeError("expected tuple")
        if len(key)!=2:
            raise TypeError("expected 2 indices")
        if key[0] >= self._m:
            raise IndexError
        if key[1] >= self._n:
            raise IndexError
        if key in self:
            return dict.__getitem__(self, key)
        return 0
