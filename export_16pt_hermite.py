from numpy import *
from hexblender.quad_interp_subdiv import quad_interp_subdiv

def subdivide_quad_16pt_hermite(f,v):

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % function [Bx By Bz ] = subdivideQuad16PtHermite(f,v)
    # % Perform Interpolatory subdivision or Catmull-Clark subdivision surface
    # % on a quad mesh and interpolate the corners and 12 other points to obtain
    # % the Geometric Coefficients of the Hermite Surface.
    # %
    # % f = num_f x 4 array of the faces
    # % v = num_v x 3 array of the vertex coordinates
    # %
    # % Bx, By, Bz are cell arrays (one for each face) of the Geometric
    # % Coefficients defining the Hermite interpolation.
    # %
    # % Note if approximating subdivisions are used this interpolates the corners
    # % of the subdivided surface not the original surface.
    # %
    # %
    # % Revision 1.0
    # % Initial Release May 2010
    # %
    # % Revision 1.1
    # % Clean mesh to remove shared vertices and duplicate faces
    # % September 2010
    # %
    # % Greg Sturgeon
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # % Clean the mesh to remove shared vertices and shared faces.
    # [v, f] = replaceSharedVertices(v,f);
    # [f] = deleteSharedFace(f);


    num_f = len(f)

    # % Perform Subdivision twice.
    fN,vN = quad_interp_subdiv(f,v)
    fN = array(fN)
    fN,vN = quad_interp_subdiv(fN,vN)
    fN = array(fN)


    # % Uncomment for Catmull-Clark (approximating) Subdivision
    # % [fN,vN] = subdivideQuad(f,v);
    # % [fN,vN] = subdivideQuad(fN,vN);


    # Bx = cell(num_f,1);
    Bx = list(range(num_f))
    for i in range(num_f):
        Bx[i] = array([[0.,0.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.]])
    Bx = array(Bx)

    # By = cell(num_f,1);
    By = list(range(num_f))
    for i in range(num_f):
        By[i] = array([[0.,0.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.]])
    By = array(By)
    
    # Bz = cell(num_f,1);
    Bz = list(range(num_f))
    for i in range(num_f):
        Bz[i] = array([[0.,0.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.]])
    Bz = array(Bz)

    Mfinv = matrix([[0., 0., 0., 1.],[ 1., 1., 1., 1.],[0., 0., 1., 0.],[3., 2., 1., 0.]])

    for k in range(num_f):
        n=(k*16)-1
        
        # % Store the locations of the 25 points from subdividing each face.
        p25x = array([[vN[fN[1+n,0],0],  vN[fN[2+n,0],0],  vN[fN[5+n,0],0],  vN[fN[6+n,0],0],   vN[fN[6+n,1],0]],
                [vN[fN[4+n,0],0],  vN[fN[3+n,0],0],  vN[fN[8+n,0],0],  vN[fN[7+n,0],0],   vN[fN[7+n,1],0]],
                [vN[fN[13+n,0],0], vN[fN[14+n,0],0], vN[fN[9+n,0],0],  vN[fN[10+n,0],0],  vN[fN[10+n,1],0]],
                [vN[fN[16+n,0],0], vN[fN[15+n,0],0], vN[fN[12+n,0],0], vN[fN[11+n,0],0],  vN[fN[11+n,1],0]],
                [vN[fN[16+n,3],0], vN[fN[15+n,3],0], vN[fN[12+n,3],0], vN[fN[11+n,3],0],  vN[fN[11+n,2],0]]])
                
        p25y = array([[vN[fN[1+n,0],1],  vN[fN[2+n,0],1],  vN[fN[5+n,0],1],  vN[fN[6+n,0],1],   vN[fN[6+n,1],1]],
                [vN[fN[4+n,0],1],  vN[fN[3+n,0],1],  vN[fN[8+n,0],1],  vN[fN[7+n,0],1],   vN[fN[7+n,1],1]],
                [vN[fN[13+n,0],1], vN[fN[14+n,0],1], vN[fN[9+n,0],1],  vN[fN[10+n,0],1],  vN[fN[10+n,1],1]],
                [vN[fN[16+n,0],1], vN[fN[15+n,0],1], vN[fN[12+n,0],1], vN[fN[11+n,0],1],  vN[fN[11+n,1],1]],
                [vN[fN[16+n,3],1], vN[fN[15+n,3],1], vN[fN[12+n,3],1], vN[fN[11+n,3],1],  vN[fN[11+n,2],1]]])

        p25z = array([[vN[fN[1+n,0],2],  vN[fN[2+n,0],2],  vN[fN[5+n,0],2],  vN[fN[6+n,0],2],   vN[fN[6+n,1],2]],
                [vN[fN[4+n,0],2],  vN[fN[3+n,0],2],  vN[fN[8+n,0],2],  vN[fN[7+n,0],2],   vN[fN[7+n,1],2]],
                [vN[fN[13+n,0],2], vN[fN[14+n,0],2], vN[fN[9+n,0],2],  vN[fN[10+n,0],2],  vN[fN[10+n,1],2]],
                [vN[fN[16+n,0],2], vN[fN[15+n,0],2], vN[fN[12+n,0],2], vN[fN[11+n,0],2],  vN[fN[11+n,1],2]],
                [vN[fN[16+n,3],2], vN[fN[15+n,3],2], vN[fN[12+n,3],2], vN[fN[11+n,3],2],  vN[fN[11+n,2],2]]])

        # % Determine the cumulative length along the curves
        lu = sqrt((p25x[:,0:4]-p25x[:,1:5])**2 + (p25y[:,0:4]-p25y[:,1:5])**2 +(p25z[:,0:4]-p25z[:,1:5])**2)
        lu = cumsum(lu,axis=1)
        
        lw = sqrt((p25x[0:4,:]-p25x[1:5,:])**2 + (p25y[0:4,:]-p25y[1:5,:])**2 +(p25z[0:4,:]-p25z[1:5,:])**2)
        lw = cumsum(lw,axis=0)        

        u = zeros((lu.shape[0],5))
        w = zeros((5,lw.shape[1]))
        
        u[:,1:5] = lu/(array([lu[:,3], lu[:,3], lu[:,3], lu[:,3]])).transpose()
        w[1:5,:] = lw/(array([lw[3,:], lw[3,:], lw[3,:], lw[3,:]]))

        # U = u([1:2,4:5],[1:2,4:5]);
        U = hstack((vstack((u[0:2,0:2],u[3:5,0:2])),vstack((u[0:2,3:5],u[3:5,3:5]))))
        W = hstack((vstack((w[0:2,0:2],w[3:5,0:2])),vstack((w[0:2,3:5],w[3:5,3:5]))))

        U = reshape(U.transpose(),16)
        W = reshape(W.transpose(),16)
        
        U = array([U**3, U**2, U, ones(U.shape,dtype=float)]).transpose()
        W = array([W**3, W**2, W, ones(W.shape,dtype=float)]).transpose()
        
        E = zeros((16,16))
        for i in range(16):
            UtW = reshape(U[i,:], (4,1))[:] * reshape(W[i,:], (1,4))[:]
            E_tmp = reshape(UtW.transpose(),(1,16))
            E[i,:] = E_tmp

        Einv = linalg.pinv(E)

        px = reshape(hstack((vstack((p25x[0:2, 0:2],p25x[3:5, 0:2])),vstack((p25x[0:2, 3:5],p25x[3:5, 3:5])))),(4,4))
        py = reshape(hstack((vstack((p25y[0:2, 0:2],p25y[3:5, 0:2])),vstack((p25y[0:2, 3:5],p25y[3:5, 3:5])))),(4,4))
        pz = reshape(hstack((vstack((p25z[0:2, 0:2],p25z[3:5, 0:2])),vstack((p25z[0:2, 3:5],p25z[3:5, 3:5])))),(4,4))

        px = reshape(px.transpose(), (16,1))
        py = reshape(py.transpose(), (16,1))
        pz = reshape(pz.transpose(), (16,1))

        ax = array(matrix(Einv)*px[:])
        ay = array(matrix(Einv)*py[:])
        az = array(matrix(Einv)*pz[:])
        
        
        # % Algebraic Coefficients
        Ax = reshape(ax,(4,4)).transpose()
        Ay = reshape(ay,(4,4)).transpose()
        Az = reshape(az,(4,4)).transpose()

        
        # % Geometric Coefficients
        Bx[k,:,:] = array(Mfinv*matrix(Ax)*Mfinv.transpose())
        By[k,:,:] = array(Mfinv*matrix(Ay)*Mfinv.transpose())
        Bz[k,:,:] = array(Mfinv*matrix(Az)*Mfinv.transpose())
        # numpy.set_printoptions(precision=4)
        # numpy.set_printoptions(suppress=True)

        # % Adjust for rounding errors in the matrix manipulation.

        Bx[k,0,0] = p25x[0,0]
        Bx[k,1,0] = p25x[0,4]
        Bx[k,0,1] = p25x[4,0]
        Bx[k,1,1] = p25x[4,4]

        By[k,0,0] = p25y[0,0]
        By[k,1,0] = p25y[0,4]
        By[k,0,1] = p25y[4,0]
        By[k,1,1] = p25y[4,4]

        Bz[k,0,0] = p25z[0,0]
        Bz[k,1,0] = p25z[0,4]
        Bz[k,0,1] = p25z[4,0]
        Bz[k,1,1] = p25z[4,4]
        
    # print "Bx\n", Bx
    # print "By\n", By
    # print "Bz\n", Bz

    return Bx, By, Bz