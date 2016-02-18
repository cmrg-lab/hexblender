"""

Created on Mon Jul 18 10:26:36 2011

@author: sturgeon
"""
import numpy

def initialize_hermite_patch(ptches,e,n,ptchFaceNum,ptchElNum):
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
    
    Mfinv = numpy.matrix([[0., 0., 0., 1.],[ 1., 1., 1., 1.],[0., 0., 1., 0.],[3., 2., 1., 0.]])
    
    Bx = numpy.zeros((ptches.shape[0],4,4),'float')
    By = numpy.zeros((ptches.shape[0],4,4),'float')
    Bz = numpy.zeros((ptches.shape[0],4,4),'float')
    
    p_u = numpy.zeros((n.shape[0]),'float')
    p_w = numpy.zeros((n.shape[0]),'float')
    p_ptchNum = numpy.zeros((n.shape[0]),'int')
    p_face = numpy.zeros((n.shape[0],4),'int')
    p_Bx = numpy.zeros((n.shape[0],4,4),'float')
    p_By = numpy.zeros((n.shape[0],4,4),'float')
    p_Bz = numpy.zeros((n.shape[0],4,4),'float')
    
    
    for k in range(ptches.shape[0]):
        fNum = ptchFaceNum[k]
        elNum = ptchElNum[k]
        if fNum == 0:
            fno0 = eno0[0,:,:]
            fno1 = eno1[0,:,:]
        elif fNum == 1:
            fno0 = eno0[4,:,:]
            fno1 = eno1[4,:,:]
        elif fNum == 4:
            fno0 = eno0[:,:,0]
            fno1 = eno1[:,:,0]
        elif fNum == 5:
            fno0 = eno0[:,:,4]
            fno1 = eno1[:,:,4]
        elif fNum == 2:
            fno0 = eno0[:,0,:]
            fno1 = eno1[:,0,:]
        elif fNum == 3:
            fno0 = eno0[:,4,:]
            fno1 = eno1[:,4,:]
        
        k_e = elNum
        ind25 = e[k_e*64+fno0,fno1]
        
        lu = numpy.sqrt(((n[ind25[:,0:4]] - n[ind25[:,1:5]])**2).sum(2))
        lw = numpy.sqrt(((n[ind25[0:4,:]] - n[ind25[1:5,:]])**2).sum(2))
        lu = numpy.cumsum(lu,axis=1)
        lw = numpy.cumsum(lw,axis=0) 
        
        u = numpy.zeros((lu.shape[0],5))
        w = numpy.zeros((5,lw.shape[1]))
        
        
        u[:,1:5] = lu/(numpy.array([lu[:,-1], lu[:,-1], lu[:,-1], lu[:,-1]])).transpose()
        w[1:5,:] = lw/(numpy.array([lw[-1,:], lw[-1,:], lw[-1,:], lw[-1,:]]))
        
        U = numpy.hstack((numpy.vstack((u[0:2,0:2],u[3:5,0:2])),numpy.vstack((u[0:2,3:5],u[3:5,3:5]))))
        W = numpy.hstack((numpy.vstack((w[0:2,0:2],w[3:5,0:2])),numpy.vstack((w[0:2,3:5],w[3:5,3:5]))))
        
        U = numpy.reshape(U.transpose(),16)
        W = numpy.reshape(W.transpose(),16)
        
        U = numpy.array([U**3, U**2, U, numpy.ones(U.shape,dtype=float)]).transpose()
        W = numpy.array([W**3, W**2, W, numpy.ones(W.shape,dtype=float)]).transpose()
        
        E = numpy.zeros((16,16))
        for i in range(16):
            UtW = numpy.reshape(U[i,:], (4,1))[:] * numpy.reshape(W[i,:], (1,4))[:]
            E_tmp = numpy.reshape(UtW.transpose(),(1,16))
            E[i,:] = E_tmp
        
        Einv = numpy.linalg.pinv(E)
        
        px = numpy.reshape(numpy.hstack((numpy.vstack((n[ind25[0:2, 0:2],0],n[ind25[3:5, 0:2],0])),numpy.vstack((n[ind25[0:2, 3:5],0],n[ind25[3:5, 3:5],0])))),(4,4))
        py = numpy.reshape(numpy.hstack((numpy.vstack((n[ind25[0:2, 0:2],1],n[ind25[3:5, 0:2],1])),numpy.vstack((n[ind25[0:2, 3:5],1],n[ind25[3:5, 3:5],1])))),(4,4))
        pz = numpy.reshape(numpy.hstack((numpy.vstack((n[ind25[0:2, 0:2],2],n[ind25[3:5, 0:2],2])),numpy.vstack((n[ind25[0:2, 3:5],2],n[ind25[3:5, 3:5],2])))),(4,4))
        
        px = numpy.reshape(px.transpose(), (16,1))
        py = numpy.reshape(py.transpose(), (16,1))
        pz = numpy.reshape(pz.transpose(), (16,1))
        
        ax = numpy.array(numpy.matrix(Einv)*px[:])
        ay = numpy.array(numpy.matrix(Einv)*py[:])
        az = numpy.array(numpy.matrix(Einv)*pz[:])
        
        # Algebraic Coefficients
        Ax = numpy.reshape(ax,(4,4)).transpose()
        Ay = numpy.reshape(ay,(4,4)).transpose()
        Az = numpy.reshape(az,(4,4)).transpose()
        
        # Geometric Coefficients
        Bx[k,:,:] = numpy.array(Mfinv*numpy.matrix(Ax)*Mfinv.transpose())
        By[k,:,:] = numpy.array(Mfinv*numpy.matrix(Ay)*Mfinv.transpose())
        Bz[k,:,:] = numpy.array(Mfinv*numpy.matrix(Az)*Mfinv.transpose())
        
        # Adjust for rounding errors in the matrix manipulation.
        
        Bx[k,0,0] = n[ind25[0,0],0]
        Bx[k,0,1] = n[ind25[4,0],0]
        Bx[k,1,0] = n[ind25[0,4],0]
        Bx[k,1,1] = n[ind25[4,4],0]
        
        By[k,0,0] = n[ind25[0,0],1]
        By[k,0,1] = n[ind25[4,0],1]
        By[k,1,0] = n[ind25[0,4],1]
        By[k,1,1] = n[ind25[4,4],1]
        
        Bz[k,0,0] = n[ind25[0,0],2]
        Bz[k,0,1] = n[ind25[4,0],2]
        Bz[k,1,0] = n[ind25[0,4],2]
        Bz[k,1,1] = n[ind25[4,4],2]
        
        # Store p.u , p.w, p.ptchnum, p.face, p.Bx, p.By, p.Bz
        for ii in range(5):
            for jj in range(5):
                p_u[ind25[ii,jj]] = u[ii,jj]
                p_w[ind25[ii,jj]] = w[ii,jj]
                p_Bx[ind25[ii,jj],:,:] = Bx[k,:,:]
                p_By[ind25[ii,jj],:,:] = By[k,:,:]
                p_Bz[ind25[ii,jj],:,:] = Bz[k,:,:]
                p_ptchNum[ind25[ii,jj]] = k
                p_face[ind25[ii,jj]] = ptches[k]
        
    return Bx,By,Bz,p_u,p_w,p_Bx,p_By,p_Bz,p_ptchNum,p_face

def hermite_points(p_bx, p_by, p_bz, p_u, p_w):
    mf = numpy.matrix(([[2, -2, 1, 1,], [-3, 3, -2, -1], [0, 0, 1, 0], [1, 0, 0, 0]]),'float')
    
    u = numpy.vstack([p_u**3,p_u**2,p_u,numpy.ones(p_u.shape)]).transpose()
    w = numpy.vstack([p_w**3,p_w**2,p_w,numpy.ones(p_w.shape)])
    
    a = u*mf
    b = mf.transpose()*w
    
    a = numpy.array(a)
    b = numpy.array(b)
    
    p = numpy.zeros((p_u.shape[0],3),'float')
    for i in range(4):
        for j in range(4):
            p[:,0] += a[:,i]*b[j,:]*p_bx[:,i,j]
            p[:,1] += a[:,i]*b[j,:]*p_by[:,i,j]
            p[:,2] += a[:,i]*b[j,:]*p_bz[:,i,j]
    
    return p

def hermite_derivative(p_bx, p_by, p_bz, p_u, p_w):
    mf = numpy.matrix(([[2, -2, 1, 1,], [-3, 3, -2, -1], [0, 0, 1, 0], [1, 0, 0, 0]]),'float')
    mfp = numpy.matrix(([[0, 0, 0, 0,], [6, -6, 3, 3], [-6, 6, -4, -2], [0, 0, 1, 0]]),'float')
    
    u = numpy.vstack([p_u**3,p_u**2,p_u,numpy.ones(p_u.shape)]).transpose()
    w = numpy.vstack([p_w**3,p_w**2,p_w,numpy.ones(p_w.shape)])
    
    # pu
    a = u*mfp
    b = mf.transpose()*w
    a = numpy.array(a)
    b = numpy.array(b)
    pu = numpy.zeros((p_u.shape[0],3),'float')
    for i in range(4):
        for j in range(4):
            pu[:,0] += a[:,i]*b[j,:]*p_bx[:,i,j]
            pu[:,1] += a[:,i]*b[j,:]*p_by[:,i,j]
            pu[:,2] += a[:,i]*b[j,:]*p_bz[:,i,j]
    
    # pw
    a = u*mf
    b = mfp.transpose()*w
    a = numpy.array(a)
    b = numpy.array(b)
    pw = numpy.zeros((p_u.shape[0],3),'float')
    for i in range(4):
        for j in range(4):
            pw[:,0] += a[:,i]*b[j,:]*p_bx[:,i,j]
            pw[:,1] += a[:,i]*b[j,:]*p_by[:,i,j]
            pw[:,2] += a[:,i]*b[j,:]*p_bz[:,i,j]
    
    return pu,pw

def hermite_second_derivative(p_bx, p_by, p_bz, p_u, p_w):
    mf = numpy.matrix(([[2, -2, 1, 1,], [-3, 3, -2, -1], [0, 0, 1, 0], [1, 0, 0, 0]]),'float')
    mfpp = numpy.matrix(([[0, 0, 0, 0,], [0, 0, 0, 0,], [12, -12, 6, 6], [-6, 6, -4, -2]]),'float')
    
    u = numpy.vstack([p_u**3,p_u**2,p_u,numpy.ones(p_u.shape)]).transpose()
    w = numpy.vstack([p_w**3,p_w**2,p_w,numpy.ones(p_w.shape)])
    
    # pu
    a = u*mfpp
    b = mf.transpose()*w
    a = numpy.array(a)
    b = numpy.array(b)
    puu = numpy.zeros((p_u.shape[0],3),'float')
    for i in range(4):
        for j in range(4):
            puu[:,0] += a[:,i]*b[j,:]*p_bx[:,i,j]
            puu[:,1] += a[:,i]*b[j,:]*p_by[:,i,j]
            puu[:,2] += a[:,i]*b[j,:]*p_bz[:,i,j]
    
    # pw
    a = u*mf
    b = mfpp.transpose()*w
    a = numpy.array(a)
    b = numpy.array(b)
    pww = numpy.zeros((p_u.shape[0],3),'float')
    for i in range(4):
        for j in range(4):
            pww[:,0] += a[:,i]*b[j,:]*p_bx[:,i,j]
            pww[:,1] += a[:,i]*b[j,:]*p_by[:,i,j]
            pww[:,2] += a[:,i]*b[j,:]*p_bz[:,i,j]
    
    return puu,pww

def normalize_v3(arrin):
    ''' Normalize a numpy array of 3 component vectors shape=(n,3) '''
    arr = arrin.copy()    
    lens = numpy.sqrt( arr[:,0]**2 + arr[:,1]**2 + arr[:,2]**2 )
    # Correct for division by zero    
    lens[numpy.nonzero(lens==0)]=1
    arr[:,0] /= lens
    arr[:,1] /= lens
    arr[:,2] /= lens                
    return arr

def magnitude_v3(arr):
    lens = numpy.sqrt( arr[:,0]**2 + arr[:,1]**2 + arr[:,2]**2 )
    lens = lens.reshape(lens.shape[0],1)
    return lens

def dot_v3(a,b):
    d = a[:,0]*b[:,0] + a[:,1]*b[:,1] + a[:,2]*b[:,2]
    d = d.reshape(a.shape[0],1)
    return d
    
def hex_volume(e,n):
    # Compute the volume of hexahedrons by creating 24 tet elements using the center point
    cP = (n[e].sum(axis=1))/8.0
    vol = numpy.zeros((numpy.size(e,0)),'float32')
    
    faces = numpy.array([[0,1,2,3],[4,5,6,7],[0,4,5,1],[3,7,6,2],[0,4,7,3],[1,5,6,2]])
    for k in range(6):
        curface = e[:,faces[k]]
        fP = (n[curface].sum(axis=1))/4.0
        for j in range(4):
            a = n[curface[:,j]]
            b = n[curface[:,numpy.mod(j+1,4)]]
            c = fP
            d = cP
            vol += abs(((a-d)*numpy.cross(b-d,c-d)).sum(1))/6
    return vol

def face_area(f,n):
    # Compute the area of a quad face by subdividing it into 4 triangles
    fP = (n[f].sum(1))/4.
    a = numpy.zeros((numpy.size(f,0)),'float32')
    for i in range(4):
        a += 0.5*numpy.sqrt((numpy.cross(fP-n[f[:,i]],fP-n[f[:,numpy.mod(i+1,4)]])**2).sum(1))
    return a


def rotate_patch(face, crossedEdge, uwdir, Bx, By, Bz,nei):
    # Rotate so that the paramater along the crossedEdge matches directions 
    # across the patches.  ie when the crossedEdge is 2 3 then the new patch 
    # should be rotated so edge 1 4 matches
    if (uwdir[0] == 1 and uwdir[1] == 2):
        fn0 = 0
        fn1 = 3
    elif (uwdir[0] == 0 and uwdir[1] == 3):
        fn0 = 1
        fn1 = 2
    elif (uwdir[0] == 3 and uwdir[1] == 2):
        fn0 = 0
        fn1 = 1
    elif (uwdir[0] == 0 and uwdir[1] == 1):
        fn0 = 3
        fn1 = 2
    rotate = 1
    while (not(face[fn0] == crossedEdge[0] and face[fn1] == crossedEdge[1]) and  rotate <=9):
        face = numpy.array([face[1], face[2], face[3], face[0]])
        nei = numpy.array([nei[1], nei[2], nei[3], nei[0]])
        rotate = rotate +1
        Bx = RotateParameterization(Bx)
        By = RotateParameterization(By)
        Bz = RotateParameterization(Bz)
        if rotate==5:
            face = numpy.array([face[3], face[2], face[1], face[0]])
#            nei = numpy.array([nei[3], nei[2], nei[1], nei[0]])
            nei = numpy.array([nei[2], nei[1], nei[0], nei[3]]) # Adjust how the neighbor is found if the patch is flipped
            Bx = FlipParamaterization(Bx)
            By = FlipParamaterization(By)
            Bz = FlipParamaterization(Bz)
    if rotate == 10:
        print("ERROR: Tryed to cross an edge boundary")
        print('face')
        print(face)
        print('crossedEdge')
        print(crossedEdge)
        print('uwdir')
        print(uwdir)
#        raw_input()
        
        face = []
        Bx = []
        By = []
        Bz = []
        nei = []
         
    return face, Bx, By, Bz, nei

def rotate_parameterization(B):
    
    Bul = B[0:2,0:2]
    Bur = B[0:2,2:4]
    Bll = B[2:4,0:2]
    Blr = B[2:4,2:4]
    
    Bn = numpy.zeros([4,4], numpy.float32)
    
    Bn[0:2,0:2] = numpy.rot90(Bul,3).copy()
    Bn[0:2,2:4] = -numpy.rot90(Bll,3).copy()
    Bn[2:4,0:2] = numpy.rot90(Bur,3).copy()
    Bn[2:4,2:4] = -numpy.rot90(Blr,3).copy()

    return Bn

    Bn = numpy.zeros([4,4], numpy.float32)
    Bn[0:2,0:2] = numpy.fliplr(B[0:2,0:2])
    Bn[0:2,2:4] = -numpy.fliplr(B[0:2,2:4])
    Bn[2:4,0:2] = numpy.fliplr(B[2:4,0:2])
    Bn[2:4,2:4] = -numpy.fliplr(B[2:4,2:4])
    return Bn
