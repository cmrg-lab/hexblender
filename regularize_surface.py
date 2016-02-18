import numpy as np
import scipy.linalg
from hexblender.quad_interp_subdiv import quad_interp_subdiv
import random

def regularize_surface(f, v, itt):

#regularizeSurface 
#Moves the vertices (v), on a smooth bicubic Hermite surface fit through 
#the quad mesh, in order to acheive more uniform edge lengths and improved
#skew of the elments.

#Note: the surface must be a manifold with consistent normal directions.  
#Open surfaces are permitted.

#Input
#f - (num_faces, 4) array of the faces of the mesh
#v - (num_vertices,3) array of the vertex location
#itt - number of iterations.
 
# Output
# f - faces of the mesh
# v - updated vertex locations
# p - cell array of structures storing the location of the vertices
#     relative to the Hermite patches.
# ptch - cell array of structures storing the Hermite surface description
#       of each face.
# gDistReg - geodesic distance for each edge in the regularized mesh
# gDistReg - geodesic distance for each edge in the original mesh

# Only function using scipy for eigh (generalized eigenvalue function)
# Could be replaced with numpy if substitute found.

# Greg Sturgeon
# December 2010

    p, ptch = move_node_to_random_connected_face(f,v)
    
    immobile_verts = np.array([ind for ind,vert in enumerate(p) if vert.isOnBnd == 1])
    
    #~ numpy.set_printoptions(precision = 4)
    #~ print "len p = ", len(p)
    #~ for i in range(len(p)):
        #~ print "i = ", i
        #~ print "xyz = ", p[i].xyz
        #~ print "u = ", p[i].u
        #~ print "w = ", p[i].w
        #~ print "Bx = "
        #~ print p[i].Bx
        #~ print "By = "
        #~ print p[i].By
        #~ print "Bz = "
        #~ print p[i].Bz
        #~ print "f = ", p[i].f
        #~ print "ptch num = ", p[i].ptch_num
        #~ print "isOnBnd = ", p[i].isOnBnd
        #~ print "isOnTwoBnd = ", p[i].isOnTwoBnd
        #~ print  "bnd edge dir = ", p[i].bnd_edge_dir
        #~ print " "

    #~ for j in range(len(ptch)):
        #~ print "j = ", j
        #~ print "patch neighbor = ", ptch[j].neighbor
        #~ print "patch boundary = ", ptch[j].bnd_edge
        #~ print "Bx = "
        #~ print ptch[j].Bx
        #~ print "By = "
        #~ print ptch[j].By
        #~ print "Bz = "
        #~ print ptch[j].Bz
        #~ print "patch f = ", ptch[j].f
        #~ print " "

    p_orig = p
    v_orig = v

    alpha = 0.05

    for k in range(itt):

        vn = laplacian_smooth(f,v,alpha,1)

        vel = v-vn
        vel[immobile_verts,:] = 0.0
        stored_vert = v[immobile_verts,:]
        
        #~ print "vel\n", vel

        # Store the max, min, and mean velocity for each iteration.
        #  vel_l = sqrt(sum(vel.^2,2));
        #  mx_vel(k) = max(vel_l);
        #  mn_vel(k) = min(vel_l);
        #  m_vel(k) = mean(vel_l);

        v,p,ptch = move_in_surface(p,ptch,vel)
        
        v[immobile_verts,:] = stored_vert
        
        if np.mod(k+1,5) == 0:
            print(" finished with iteration #%d" % (k+1))
    
    v_new = v.copy()
    f_new = f.copy()
    
    fNlist = []
    for Array in f_new:
        fNlist.append(list(Array))
    for i in range(len(fNlist)):
        for j in range(len(fNlist[i])):
            fNlist[i][j] = int(fNlist[i][j])
    
    return fNlist, v_new

def move_node_to_random_connected_face(f, v):

#   - Called by regularize Surface
#
# Initializes the Hermite patch description of the element faces and 
# vertices. The vertices are assigned to be defined relative to one of it 
# connected faces.
#
# Input
# f - (num_faces, 4) array of the faces of the mesh
# v - (num_vertices,3) array of the vertex location
# 
# Output
# p - cell array of structures storing the location of the vertices
#     relative to the Hermite patches.
# ptch - cell array of structures storing the Hermite surface description
#       of each face.
# bndCrvs - cell array with the indices to vertices on the boundary for
#       each boundary curve.
#
# Greg Sturgeon
# December 2010

    num_v = np.shape(v)[0]
    num_f = np.shape(f)[0]

    #~ % Find Connected Faces
    #~ % connFs is a cell array storing all the faces connected to a given vertex
    # connFs = cell(num_v,1);
    connFs = list(range(num_v))
    for i in range(num_v):
        connFs[i] = []

    for i in range(num_f):

        if i not in connFs[f[i,0]]:
           connFs[f[i,0]].append(i)

        if i not in connFs[f[i,1]]:
           connFs[f[i,1]].append(i)

        if i not in connFs[f[i,2]]:
           connFs[f[i,2]].append(i)

        if i not in connFs[f[i,3]]:
           connFs[f[i,3]].append(i)

    #print statement for debug
    #~ for i in range(num_v):
        #~ print connFs[i]
    
    #~ % Create Hermite Patch
    #~ % Bx,By,Bz are cell arrays of the geometric coefficients of 
    #~ % Hermite patches created for each face.

    #~ % The patch is paramaterized such that the face edge 1->2 
    #~ % defines the u direction and 1->4 defines the w direction.

    #~ % [Bx By Bz] = subdivideQuad16PtHermite(f,v);
    Bx, By, Bz = subdivide_quad16_pt_hermite81(f,v)

    epsilon = 0

    # p = cell(num_v,1);
    #~ p = list(range(num_v))
    #~ for i in range(num_v):
        #~ p[i] = []
        
    class PClass():
        def __init__(self, xyz, u, w, Bx, By, Bz, f, ptch_num):
            self.xyz = xyz
            self.u = u
            self.w = w
            self.Bx = Bx
            self.By = By
            self.Bz = Bz
            self.f = f
            self.ptch_num = ptch_num
            self.isOnBnd = 0
            self.isOnTwoBnd = 0
            self.bnd_edge_dir = [0,0]
            self.neighbor = []
            self.bnd_edge = []
    
    p = []
    #~ for i in range(len(p)):
        #~ p[i] = PClass
    #~ p = [PClass() for val in range(num_v)]
    
    #~ print "connFs = ", connFs
    
    for i in range(num_v):
        
        
        rn = random.randrange(len(connFs[i]))
        ptch_num = connFs[i][rn]
        
        #print "ptch num = ", ptch_num
       
        # rotate the patch such that v{i} is f(_,1)
       
        Bx_i, By_i, Bz_i, f_rot  = rotate_patch(f[ptch_num,:],
                                                i,
                                                Bx[ptch_num,:,:],
                                                By[ptch_num,:,:],
                                                Bz[ptch_num,:,:])
       
        #print "f_rot = ", f_rot
        
        xyz = hermite_point(Bx_i, By_i, Bz_i,epsilon,epsilon)
        
        #print "p xyz = ", xyz
        
        p.append(PClass(xyz, 
                        epsilon, 
                        epsilon, 
                        Bx_i, 
                        By_i, 
                        Bz_i, 
                        f_rot, 
                        ptch_num))
       
        #~ p[i].u = epsilon
        #~ p[i].w = epsilon
       
        #~ p[i].Bx = Bx_i
        #~ p[i].By = By_i
        #~ p[i].Bz = Bz_i
       
        #~ p[i].f = f_rot
       
        #~ p[i].ptch_num = ptch_num
    
    # Create Adjacency Matrix

    # eM1 = zeros((num_v,num_v))
    # eM2 = zeros((num_v,num_v))

    eM1 = Matrix(num_v,num_v)
    eM2 = Matrix(num_v,num_v)

    # Mark an edge (vertex pair) in edgeMatrix1 (eM1) the first time an edge is 
    # used and edgeMatrix2 (eM2)the second time.  The edgeMatrix contains the 
    # face # that uses the edge.

    for i in range(num_f):
        
        if not eM1[np.minimum(f[i,0],f[i,1]),np.maximum(f[i,0],f[i,1])]:
            eM1[np.minimum(f[i,0],f[i,1]),np.maximum(f[i,0],f[i,1])]=i+1
        else:
            eM2[np.minimum(f[i,0],f[i,1]),np.maximum(f[i,0],f[i,1])]=i+1

        if not eM1[np.minimum(f[i,1],f[i,2]),np.maximum(f[i,1],f[i,2])]:
            eM1[np.minimum(f[i,1],f[i,2]),np.maximum(f[i,1],f[i,2])]=i+1
        else:
            eM2[np.minimum(f[i,1],f[i,2]),np.maximum(f[i,1],f[i,2])]=i+1

        if not eM1[np.minimum(f[i,2],f[i,3]),np.maximum(f[i,2],f[i,3])]:
            eM1[np.minimum(f[i,2],f[i,3]),np.maximum(f[i,2],f[i,3])]=i+1
        else:
            eM2[np.minimum(f[i,2],f[i,3]),np.maximum(f[i,2],f[i,3])]=i+1

        if not eM1[np.minimum(f[i,3],f[i,0]),np.maximum(f[i,3],f[i,0])]:
            eM1[np.minimum(f[i,3],f[i,0]),np.maximum(f[i,3],f[i,0])]=i+1
        else:
            eM2[np.minimum(f[i,3],f[i,0]),np.maximum(f[i,3],f[i,0])]=i+1
    
    #~ print "eM1 = "
    #~ print eM1
    #~ print "eM2 = "
    #~ print eM2

    #~ for i in range(num_v):
        #~ p[i].isOnBnd = 0
        #~ p[i].isOnTwoBnd = 0
        #~ p[i].bnd_edge_dir = [0, 0]

    # ptch = cell(num_f,1);
    #~ ptch = list(range(num_v))
    #~ for i in range(num_v):
        #~ ptch[i] = []
    
    class PatchClass:
        def __init__(self,Bx,By,Bz,f):
            self.neighbor = []
            self.bnd_edge = []
            self.Bx = Bx
            self.By = By
            self.Bz = Bz
            self.f = f
    
    ptch = [PatchClass(Bx[i,:,:],
                       By[i,:,:],
                       Bz[i,:,:],
                       f[i,:]) for i in range(num_f)]
        
    #~ for i in range(len(ptch)):
        #~ print "i = ", i
        #~ print "ptch neigh = ", ptch[i].neighbor
        #~ print "ptch bnd edge = ", ptch[i].bnd_edge
        #~ print "f = ", ptch[i].f
    
    #~ for i in range(len(ptch)):
        #~ ptch[i] = PClass
            
    #~ for i in range(num_f):
        #~ ptch[i].neighbor = []
        #~ ptch[i].bnd_edge = []
        #~ ptch[i].Bx = Bx[i]
        #~ ptch[i].By = By[i]
        #~ ptch[i].Bz = Bz[i]
        #~ ptch[i].f = f[i,:]

    e_boundary = np.zeros((num_f,2), np.int32)
    num_bnd_e = 0

    for i in range(num_f): 
        for j in range(4):
            # Store the patches that are neighbors to each patch
            nei_1 = eM1[np.minimum(f[i,j],f[i,(j+1)%4]),
                        np.maximum(f[i,j],f[i,(j+1)%4])] 
            if nei_1 and nei_1 != (i+1):
                ptch[i].neighbor.append(nei_1-1) #-1 for python?
            
            nei_2 = eM2[np.minimum(f[i,j],f[i,(j+1)%4]),
                        np.maximum(f[i,j],f[i,(j+1)%4])] #-1 for python?
            if nei_2 and nei_2 != (i+1):
                ptch[i].neighbor.append(nei_2-1) #-1 for python?
            
            # Identify which patches have boundaries and which of their edges 
            # are the boundary.

            if not nei_2: # Edge only has one face
                ptch[i].bnd_edge.append([np.minimum(f[i,j],f[i,(j+1)%4]), 
                                         np.maximum(f[i,j],f[i,(j+1)%4])])
                try:
                    e_boundary[num_bnd_e,:] = [np.minimum(f[i,j],f[i,(j+1)%4]), 
                                               np.maximum(f[i,j],f[i,(j+1)%4])]
                except IndexError:
                    e_boundary = np.vstack((e_boundary, 
                                         np.array([np.minimum(f[i,j],f[i,(j+1)%4]), 
                                                np.maximum(f[i,j],f[i,(j+1)%4])])))
                num_bnd_e = num_bnd_e+1

    #~ e_boundary[num_bnd_e:,:]=[]
    #slice it the other way since there is no "empty" in python
    e_boundary = e_boundary[:num_bnd_e,:]
    
    #~ print "e_boundary = ", e_boundary

    bndCrvs = find_boundary_crvs(e_boundary)
    
    #~ print "bndCrvs = ", bndCrvs
    
    for c in range(len(bndCrvs)):
        #print "c = ", c
        curCrv = bndCrvs[c]
        #print "curCrv = ", curCrv
        for k in range(len(curCrv)):
            if not p[curCrv[k]].isOnBnd:
                p[curCrv[k]].isOnBnd = c+1 #HMMM CHANGE BACK TO "c" LATER?
                
                # this "if" fixed a bug that caused an error when 
                # "extraordinary nodes on the boundary" exist in the mesh --
                # the indexing scheme has problems because the current 
                # data structure construction doesn't recognize that a 
                # face with only one vertex on the boundary, is on the 
                # boundary.  It seems that this "if" statement looking 
                # for an empty list and just skipping the steps works 
                # fine, but I only tested it on a couple of cases.  MJG 12/3/11
                if ptch[p[curCrv[k]].ptch_num].bnd_edge:
                    be = ptch[p[curCrv[k]].ptch_num].bnd_edge[0]
                    
                    #print "be = ", be
                    
                    Bx, By, Bz, face, u, w = rotate_patch_bnd(p[curCrv[k]].f,
                                                              be,
                                                              p[curCrv[k]].Bx,
                                                              p[curCrv[k]].By,
                                                              p[curCrv[k]].Bz,
                                                              p[curCrv[k]].u,
                                                              p[curCrv[k]].w)
                    
                    p[curCrv[k]].Bx = Bx
                    p[curCrv[k]].By = By
                    p[curCrv[k]].Bz = Bz
                    p[curCrv[k]].f = face
                    p[curCrv[k]].u = u
                    p[curCrv[k]].w = w
                
            else:
                print("A node is on two different boundary curves")

    return p, ptch
    #~ return p, ptch, bndCrvs

def rotate_patch(face, v_num, Bx, By, Bz):

    # Rotate the face counterclockwise until the vertex v_num is the first 
    # vertex of the face.  Reparamaterize the hermite patch with each rotation.

    rotate = 1
    while face[0] != v_num and rotate <=4:
        face = np.array([face[1], face[2], face[3], face[0]])
        rotate = rotate +1
        Bx = rotate_parameterization(Bx)
        By = rotate_parameterization(By)
        Bz = rotate_parameterization(Bz)
        
    return Bx, By, Bz, face

def rotate_patch_bnd(face, be, Bx, By, Bz, u, w):

# Rotate the face counterclockwise until the vertex v_num is the first 
# vertex of the face.  Reparamaterize the hermite patch with each rotation.

    rotate = 1
    while not (((face[0] == be[0] or face[0] == be[1]) and 
        (face[1] == be[0] or face[1] == be[1])) and (rotate <=5)):
            face = [face[1], face[2], face[3], face[0]]
            if rotate == 1 or rotate == 3:
                w = 1-w                            #### perhaps '-1' in python
            elif rotate == 2 or rotate == 4:
                u = 1-u
            
            rotate = rotate +1
            Bx = rotate_parameterization(Bx)
            By = rotate_parameterization(By)
            Bz = rotate_parameterization(Bz)
    
    if rotate ==6:
        print("Error: Rotated 360 degrees")
        
    return Bx, By, Bz, face, u, w

def rotate_parameterization(B):
    
    Bul = B[0:2,0:2]
    Bur = B[0:2,2:4]
    Bll = B[2:4,0:2]
    Blr = B[2:4,2:4]
    
    Bn = np.zeros([4,4], np.float32)
    
    Bn[0:2,0:2] = np.rot90(Bul,3).copy()
    Bn[0:2,2:4] = -np.rot90(Bll,3).copy()
    Bn[2:4,0:2] = np.rot90(Bur,3).copy()
    Bn[2:4,2:4] = -np.rot90(Blr,3).copy()

    return Bn

def find_boundary_crvs(e_boundary):

    # bndEConn = zeros(maximum(maximum(e_boundary)),1)
    
    #WARNING : e_boundary ALWAYS HAS SECOND DIMENSION = 2, RIGHT?
    
    bndEConn = [ [] for val in range(int(np.max(e_boundary))+1) ]
    bndCrvs = [ [] for val in range(int(np.max(e_boundary))+1) ]

    for k in range(int(np.shape(e_boundary)[0])):
        # bndEConn[e_boundary[k,0]](end+1) = e_boundary(k,2);
        # bndEConn[e_boundary[k,1]](end+1) = e_boundary(k,1);
        bndEConn[e_boundary[k,0]].append(e_boundary[k,1])
        bndEConn[e_boundary[k,1]].append(e_boundary[k,0])
    
    #~ print "bndEConn = "
    #~ print bndEConn
    
    # Group the connected curves together to form closed curve chains
    c = 0
    moreCrvs = 1
    # bndCrvs = cell(0);
    bndCrvs = []
    
    #~ next_node = None
    
    while moreCrvs:
        moreCrvs = 0
        for k in range(len(bndEConn)):
            if bndEConn[k] != []:
                cur_node = k
                moreCrvs = 1
                break
        
        #~ print "cur_node = ", k
        
        if moreCrvs:
            bndCrvs.append([])
            bndCrvs[c].append(cur_node)
            #~ bndCrvs[c][0] = cur_node
            next_node = bndEConn[cur_node][0]
            #~ bndEConn[cur_node][0] = []
            del bndEConn[cur_node][0]
            
            #~ print "next_node = ", next_node

            while next_node:
                prev_node = cur_node;
                cur_node = next_node; 
                #~ next_node = setdiff1d(bndEConn[next_node],bndCrvs[c])
                next_node = np.setdiff1d(bndEConn[next_node],bndCrvs[c])
                next_node = next_node[0] if next_node else next_node
                bndCrvs[c].append(cur_node)
                
                #~ print "isempty prev node = ", prev_node
                #~ print "isempty cur node = ", cur_node
                #~ print "isempty next node = ", next_node
                #~ print " "

                if bndEConn[cur_node][0] == prev_node:
                    #~ bndEConn[cur_node][0] = []
                    del bndEConn[cur_node][0]
                elif bndEConn[cur_node][1] == prev_node:
                    #~ bndEConn[cur_node][1] = []
                    del bndEConn[cur_node][1]

            for j in range(len(bndCrvs[c])):
                bndEConn[bndCrvs[c][j]] = []
                #~ del bndEConn[bndCrvs[c][j]]

            c=c+1
            #print "c = ", c
    
    #~ print "bndEConn = ", bndEConn
    #~ print "bndCrvs = ", bndCrvs
    
    return bndCrvs

def subdivide_quad16_pt_hermite81(f,v):

    # function [Bx By Bz ] = subdivideQuad16PtHermite(f,v)
    # Perform Interpolatory subdivision or Catmull-Clark subdivision surface  
    # on a quad mesh and interpolate the corners and 12 other points to obtain 
    # the Geometric Coefficients of the Hermite Surface.
    #
    # f = num_f x 4 array of the faces
    # v = num_v x 3 array of the vertex coordinates
    #
    # Bx, By, Bz are cell arrays (one for each face) of the Geometric
    # Coefficients defining the Hermite interpolation.
    #
    # Note if approximating subdivisions are used this interpolates the corners 
    # of the subdivided surface not the original surface. 
    #
    # 
    # Revision 1.0 
    # Initial Release May 2010
    #
    # Revision 1.1
    # Clean mesh to remove shared vertices and duplicate faces
    # September 2010
    #
    # Greg Sturgeon

    num_f = len(f)

    # Perform Subdivision twice (3 times?).
    fN,vN = quad_interp_subdiv(f,v)
    fN = np.array(fN)
    fN,vN = quad_interp_subdiv(fN,vN)
    fN = np.array(fN)
    fN,vN = quad_interp_subdiv(fN,vN)
    fN = np.array(fN)

    # Bx = cell(num_f,1);
    Bx = list(range(num_f))
    for i in range(num_f):
        Bx[i] = np.array([[0.,0.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.]])
    Bx = np.array(Bx)

    # By = cell(num_f,1);
    By = list(range(num_f))
    for i in range(num_f):
        By[i] = np.array([[0.,0.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.]])
    By = np.array(By)

    # Bz = cell(num_f,1);
    Bz = list(range(num_f))
    for i in range(num_f):
        Bz[i] = np.array([[0.,0.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.]])
    Bz = np.array(Bz)

    Mfinv = np.matrix([[0., 0., 0., 1.],[ 1., 1., 1., 1.],[0., 0., 1., 0.],[3., 2., 1., 0.]])
    
    for k in range(num_f):
        n=(k*64)-1
        
        # % Store the locations of the 25 points from subdividing each face.
        p81x =  np.array([[vN[fN[1+n,0],0],         
                        vN[fN[2+n,0],0],        
                        vN[fN[5+n,0],0],        
                        vN[fN[6+n,0],0],         
                        vN[fN[1+16+n,0],0],    
                        vN[fN[2+16+n,0],0],   
                        vN[fN[5+16+n,0],0],    
                        vN[fN[6+16+n,0],0],   
                        vN[fN[6+16+n,1],0]  ],
                       [vN[fN[4+n,0],0],         
                        vN[fN[3+n,0],0],        
                        vN[fN[8+n,0],0],        
                        vN[fN[7+n,0],0],         
                        vN[fN[4+16+n,0],0],   
                        vN[fN[3+16+n,0],0],    
                        vN[fN[8+16+n,0],0],    
                        vN[fN[7+16+n,0],0],   
                        vN[fN[7+16+n,1],0]  ],
                       [vN[fN[13+n,0],0],       
                        vN[fN[14+n,0],0],       
                        vN[fN[9+n,0],0],        
                        vN[fN[10+n,0],0],       
                        vN[fN[13+16+n,0],0],  
                        vN[fN[14+16+n,0],0],  
                        vN[fN[9+16+n,0],0],    
                        vN[fN[10+16+n,0],0], 
                        vN[fN[10+16+n,1],0] ],
                       [vN[fN[16+n,0],0],       
                        vN[fN[15+n,0],0],       
                        vN[fN[12+n,0],0],      
                        vN[fN[11+n,0],0],       
                        vN[fN[16+16+n,0],0],  
                        vN[fN[15+16+n,0],0],  
                        vN[fN[12+16+n,0],0],   
                        vN[fN[11+16+n,0],0], 
                        vN[fN[11+16+n,1],0] ],
                       [vN[fN[1+48+n,0],0],   
                        vN[fN[2+48+n,0],0],   
                        vN[fN[5+48+n,0],0],   
                        vN[fN[6+48+n,0],0],   
                        vN[fN[1+32+n,0],0],   
                        vN[fN[2+32+n,0],0],    
                        vN[fN[5+32+n,0],0],    
                        vN[fN[6+32+n,0],0],   
                        vN[fN[6+32+n,1],0]  ],
                       [vN[fN[4+48+n,0],0],   
                        vN[fN[3+48+n,0],0],   
                        vN[fN[8+48+n,0],0],   
                        vN[fN[7+48+n,0],0],  
                        vN[fN[4+32+n,0],0],    
                        vN[fN[3+32+n,0],0],   
                        vN[fN[8+32+n,0],0],    
                        vN[fN[7+32+n,0],0],    
                        vN[fN[7+32+n,1],0]  ],
                       [vN[fN[13+48+n,0],0],  
                        vN[fN[14+48+n,0],0], 
                        vN[fN[9+48+n,0],0],   
                        vN[fN[10+48+n,0],0], 
                        vN[fN[13+32+n,0],0],  
                        vN[fN[14+32+n,0],0],  
                        vN[fN[9+32+n,0],0],   
                        vN[fN[10+32+n,0],0],  
                        vN[fN[10+32+n,1],0] ],
                       [vN[fN[16+48+n,0],0],  
                        vN[fN[15+48+n,0],0], 
                        vN[fN[12+48+n,0],0], 
                        vN[fN[11+48+n,0],0], 
                        vN[fN[16+32+n,0],0],  
                        vN[fN[15+32+n,0],0],  
                        vN[fN[12+32+n,0],0],  
                        vN[fN[11+32+n,0],0],  
                        vN[fN[11+32+n,1],0] ],
                       [vN[fN[16+48+n,3],0],  
                        vN[fN[15+48+n,3],0], 
                        vN[fN[12+48+n,3],0], 
                        vN[fN[11+48+n,3],0], 
                        vN[fN[16+32+n,3],0],  
                        vN[fN[15+32+n,3],0],  
                        vN[fN[12+32+n,3],0],  
                        vN[fN[11+32+n,3],0],  
                        vN[fN[11+32+n,2],0] ]])
        
        #~ print "p81x = ", p81x
        
        p81y = np.array([[vN[fN[1+n,0],1], 
                       vN[fN[2+n,0],1], 
                       vN[fN[5+n,0],1], 
                       vN[fN[6+n,0],1], 
                       vN[fN[1+16+n,0],1], 
                       vN[fN[2+16+n,0],1], 
                       vN[fN[5+16+n,0],1], 
                       vN[fN[6+16+n,0],1], 
                       vN[fN[6+16+n,1],1]],
                      [vN[fN[4+n,0],1], 
                       vN[fN[3+n,0],1], 
                       vN[fN[8+n,0],1], 
                       vN[fN[7+n,0],1], 
                       vN[fN[4+16+n,0],1], 
                       vN[fN[3+16+n,0],1], 
                       vN[fN[8+16+n,0],1], 
                       vN[fN[7+16+n,0],1], 
                       vN[fN[7+16+n,1],1]],
                      [vN[fN[13+n,0],1], 
                       vN[fN[14+n,0],1], 
                       vN[fN[9+n,0],1], 
                       vN[fN[10+n,0],1], 
                       vN[fN[13+16+n,0],1], 
                       vN[fN[14+16+n,0],1], 
                       vN[fN[9+16+n,0],1], 
                       vN[fN[10+16+n,0],1], 
                       vN[fN[10+16+n,1],1]],
                      [vN[fN[16+n,0],1], 
                       vN[fN[15+n,0],1], 
                       vN[fN[12+n,0],1], 
                       vN[fN[11+n,0],1], 
                       vN[fN[16+16+n,1],2], 
                       vN[fN[15+16+n,0],1], 
                       vN[fN[12+16+n,0],1], 
                       vN[fN[11+16+n,0],1], 
                       vN[fN[11+16+n,1],1]],
                      [vN[fN[1+48+n,0],1], 
                       vN[fN[2+48+n,0],1], 
                       vN[fN[5+48+n,0],1], 
                       vN[fN[6+48+n,0],1], 
                       vN[fN[1+32+n,0],1], 
                       vN[fN[2+32+n,0],1], 
                       vN[fN[5+32+n,0],1], 
                       vN[fN[6+32+n,0],1], 
                       vN[fN[6+32+n,1],1]],
                      [vN[fN[4+48+n,0],1], 
                       vN[fN[3+48+n,0],1], 
                       vN[fN[8+48+n,0],1], 
                       vN[fN[7+48+n,0],1], 
                       vN[fN[4+32+n,0],1], 
                       vN[fN[3+32+n,0],1], 
                       vN[fN[8+32+n,0],1], 
                       vN[fN[7+32+n,0],1], 
                       vN[fN[7+32+n,1],1]],
                      [vN[fN[13+48+n,0],1], 
                       vN[fN[14+48+n,0],1], 
                       vN[fN[9+48+n,0],1], 
                       vN[fN[10+48+n,0],1], 
                       vN[fN[13+32+n,0],1], 
                       vN[fN[14+32+n,0],1], 
                       vN[fN[9+32+n,0],1], 
                       vN[fN[10+32+n,0],1], 
                       vN[fN[10+32+n,1],1]],
                      [vN[fN[16+48+n,0],1], 
                       vN[fN[15+48+n,0],1], 
                       vN[fN[12+48+n,0],1], 
                       vN[fN[11+48+n,0],1], 
                       vN[fN[16+32+n,0],1], 
                       vN[fN[15+32+n,0],1], 
                       vN[fN[12+32+n,0],1], 
                       vN[fN[11+32+n,0],1], 
                       vN[fN[11+32+n,1],1]],
                      [vN[fN[16+48+n,3],1], 
                       vN[fN[15+48+n,3],1], 
                       vN[fN[12+48+n,3],1], 
                       vN[fN[11+48+n,3],1], 
                       vN[fN[16+32+n,3],1], 
                       vN[fN[15+32+n,3],1], 
                       vN[fN[12+32+n,3],1], 
                       vN[fN[11+32+n,3],1], 
                       vN[fN[11+32+n,2],1]]])
        
        #~ print "p81y = ", p81y
        
        p81z = np.array([[vN[fN[1+n,0],2], 
                       vN[fN[2+n,0],2], 
                       vN[fN[5+n,0],2], 
                       vN[fN[6+n,0],2], 
                       vN[fN[1+16+n,0],2], 
                       vN[fN[2+16+n,0],2], 
                       vN[fN[5+16+n,0],2], 
                       vN[fN[6+16+n,0],2], 
                       vN[fN[6+16+n,1],2]],
                      [vN[fN[4+n,0],2], 
                       vN[fN[3+n,0],2], 
                       vN[fN[8+n,0],2], 
                       vN[fN[7+n,0],2], 
                       vN[fN[4+16+n,0],2], 
                       vN[fN[3+16+n,0],2], 
                       vN[fN[8+16+n,0],2], 
                       vN[fN[7+16+n,0],2], 
                       vN[fN[7+16+n,1],2]],
                      [vN[fN[13+n,0],2], 
                       vN[fN[14+n,0],2], 
                       vN[fN[9+n,0],2], 
                       vN[fN[10+n,0],2], 
                       vN[fN[13+16+n,0],2], 
                       vN[fN[14+16+n,0],2], 
                       vN[fN[9+16+n,0],2], 
                       vN[fN[10+16+n,0],2], 
                       vN[fN[10+16+n,1],2]],
                      [vN[fN[16+n,0],2], 
                       vN[fN[15+n,0],2], 
                       vN[fN[12+n,0],2], 
                       vN[fN[11+n,0],2], 
                       vN[fN[16+16+n,0],2], 
                       vN[fN[15+16+n,0],2], 
                       vN[fN[12+16+n,0],2], 
                       vN[fN[11+16+n,0],2], 
                       vN[fN[11+16+n,1],2]],
                      [vN[fN[1+48+n,0],2], 
                       vN[fN[2+48+n,0],2], 
                       vN[fN[5+48+n,0],2], 
                       vN[fN[6+48+n,0],2], 
                       vN[fN[1+32+n,0],2], 
                       vN[fN[2+32+n,0],2], 
                       vN[fN[5+32+n,0],2], 
                       vN[fN[6+32+n,0],2], 
                       vN[fN[6+32+n,1],2]],
                      [vN[fN[4+48+n,0],2], 
                       vN[fN[3+48+n,0],2], 
                       vN[fN[8+48+n,0],2], 
                       vN[fN[7+48+n,0],2], 
                       vN[fN[4+32+n,0],2], 
                       vN[fN[3+32+n,0],2], 
                       vN[fN[8+32+n,0],2], 
                       vN[fN[7+32+n,0],2], 
                       vN[fN[7+32+n,1],2]],
                      [vN[fN[13+48+n,0],2], 
                       vN[fN[14+48+n,0],2], 
                       vN[fN[9+48+n,0],2], 
                       vN[fN[10+48+n,0],2], 
                       vN[fN[13+32+n,0],2], 
                       vN[fN[14+32+n,0],2], 
                       vN[fN[9+32+n,0],2], 
                       vN[fN[10+32+n,0],2], 
                       vN[fN[10+32+n,1],2]],
                      [vN[fN[16+48+n,0],2], 
                       vN[fN[15+48+n,0],2], 
                       vN[fN[12+48+n,0],2], 
                       vN[fN[11+48+n,0],2], 
                       vN[fN[16+32+n,0],2], 
                       vN[fN[15+32+n,0],2], 
                       vN[fN[12+32+n,0],2], 
                       vN[fN[11+32+n,0],2], 
                       vN[fN[11+32+n,1],2]],
                      [vN[fN[16+48+n,3],2], 
                       vN[fN[15+48+n,3],2], 
                       vN[fN[12+48+n,3],2], 
                       vN[fN[11+48+n,3],2], 
                       vN[fN[16+32+n,3],2], 
                       vN[fN[15+32+n,3],2], 
                       vN[fN[12+32+n,3],2], 
                       vN[fN[11+32+n,3],2], 
                       vN[fN[11+32+n,2],2]]])
        
        #~ print "p81z = ", p81z

        # Determine the cumulative length along the curves    
        lu = np.sqrt((p81x[:,0:8]-p81x[:,1:9])**2 + (p81y[:,0:8]-p81y[:,1:9])**2 +(p81z[:,0:8]-p81z[:,1:9])**2)
        lu = np.cumsum(lu,axis=1)
        lw = np.sqrt((p81x[0:8,:]-p81x[1:9,:])**2 + (p81y[0:8,:]-p81y[1:9,:])**2 + (p81z[0:8,:]-p81z[1:9,:])**2)
        lw = np.cumsum(lw,axis=0) 
        
        #~ numpy.set_printoptions(linewidth=125)
        #~ numpy.set_printoptions(precision=4)
        #~ print "lu = "
        #~ print array(lu)
        #~ print "lw = "
        #~ print array(lw)
        
        u = np.zeros((lu.shape[0],9))
        w = np.zeros((9,lw.shape[1]))
        
        #~ print "shape u = ", shape(u)
        #~ print "shape w = ", shape(w)
        
        u[:,1:9] = lu/(np.array([lu[:,-1], lu[:,-1], lu[:,-1], lu[:,-1], lu[:,-1], lu[:,-1], lu[:,-1], lu[:,-1]])).transpose()
        w[1:9,:] = lw/(np.array([lw[-1,:], lw[-1,:], lw[-1,:], lw[-1,:], lw[-1,:], lw[-1,:], lw[-1,:], lw[-1,:]]))
        
        #~ print "u = "
        #~ print u
        #~ print "w = "
        #~ print w

        # U = u([1:2,8:9],[1:2,8:9]);
        # W = w([1:2,8:9],[1:2,8:9]);
        U = np.hstack((np.vstack((u[0:2,0:2],u[7:9,0:2])),np.vstack((u[0:2,7:9],u[7:9,7:9]))))
        W = np.hstack((np.vstack((w[0:2,0:2],w[7:9,0:2])),np.vstack((w[0:2,7:9],w[7:9,7:9]))))
        
        #~ print "U #1 = "
        #~ print U
        #~ print "W #1 = "
        #~ print W

        U = np.reshape(U.transpose(),16)
        W = np.reshape(W.transpose(),16)
        
        #~ print "U #2 = "
        #~ print U
        #~ print "W #2 = "
        #~ print W
        
        U = np.array([U**3, U**2, U, np.ones(U.shape,dtype=float)]).transpose()
        W = np.array([W**3, W**2, W, np.ones(W.shape,dtype=float)]).transpose()
        
        #~ print "U #3 = "
        #~ print U
        #~ print "W #3 = "
        #~ print W

        E = np.zeros((16,16))
        for i in range(16):
            UtW = np.reshape(U[i,:], (4,1))[:] * np.reshape(W[i,:], (1,4))[:]
            E_tmp = np.reshape(UtW.transpose(),(1,16))
            E[i,:] = E_tmp

        #Einv = np.linalg.pinv(E)
        Einv = scipy.linalg.pinv(E)
        
        #NEVER DOUBLE-CHECKED THESE AREAS ::> OJO
        
        px = np.reshape(np.hstack((np.vstack((p81x[0:2, 0:2],p81x[7:9, 0:2])),
                     np.vstack((p81x[0:2, 7:9],p81x[7:9, 7:9])))),(4,4))
        py = np.reshape(np.hstack((np.vstack((p81y[0:2, 0:2],p81y[7:9, 0:2])),
                     np.vstack((p81y[0:2, 7:9],p81y[7:9, 7:9])))),(4,4))
        pz = np.reshape(np.hstack((np.vstack((p81z[0:2, 0:2],p81z[7:9, 0:2])),
                     np.vstack((p81z[0:2, 7:9],p81z[7:9, 7:9])))),(4,4))

        px = np.reshape(px.transpose(), (16,1))
        py = np.reshape(py.transpose(), (16,1))
        pz = np.reshape(pz.transpose(), (16,1))

        ax = np.array(np.matrix(Einv)*px[:])
        ay = np.array(np.matrix(Einv)*py[:])
        az = np.array(np.matrix(Einv)*pz[:])

        # Algebraic Coefficients
        Ax = np.reshape(ax,(4,4)).transpose()
        Ay = np.reshape(ay,(4,4)).transpose()
        Az = np.reshape(az,(4,4)).transpose()

        # Geometric Coefficients
        Bx[k,:,:] = np.array(Mfinv * np.matrix(Ax) * Mfinv.transpose())
        By[k,:,:] = np.array(Mfinv * np.matrix(Ay) * Mfinv.transpose())
        Bz[k,:,:] = np.array(Mfinv * np.matrix(Az) * Mfinv.transpose())
        
        # Adjust for rounding errors in the matrix manipulation.
        
        Bx[k,0,0] = p81x[0,0]
        Bx[k,1,0] = p81x[0,8]
        Bx[k,0,1] = p81x[8,0]
        Bx[k,1,1] = p81x[8,8]

        By[k,0,0] = p81y[0,0]
        By[k,1,0] = p81y[0,8]
        By[k,0,1] = p81y[8,0]
        By[k,1,1] = p81y[8,8]

        Bz[k,0,0] = p81z[0,0]
        Bz[k,1,0] = p81z[0,8]
        Bz[k,0,1] = p81z[8,0]
        Bz[k,1,1] = p81z[8,8]
        
    return Bx, By, Bz

def laplacian_smooth(f,v,alpha,itt):
    # LaplacianSmooth(f,v,alpha,itt)
    # 
    # Perform Laplacian Smoothing on the mesh, moving v towards a weighted 
    # average of its neighbors.  The weighting is porportional to the distance
    # between the point and the neighboring point.
    #
    # Input
    # f - (num_faces, 4) array of the faces of the mesh
    # v - (num_vertices,3) array of the vertex location
    # alpha - small positive constant
    # itt - number of iterations
    #
    # Output
    # p - (num_vertices,3) array of the updated position of v
    #
    # Greg Sturgeon
    # December 2010

    p=v.copy()

    num_v = len(v)
    num_f = len(f)
    # edgeConn = cell(num_v,1);
    edgeConn = list(range(num_v))
    for i in range(num_v):
        edgeConn[i] = []

    for i in range(num_f):

       # edgeConn is a cell array storing all the vertices connected to a 
       # given vertex.  edgeConn{k} has an array of all the vertices that are 
       # connected to vertex k.  (There may be duplicates stored).
        
        edgeConn[f[i,0]].append(f[i,3])
        edgeConn[f[i,0]].append(f[i,1])
       
        edgeConn[f[i,1]].append(f[i,0])
        edgeConn[f[i,1]].append(f[i,2])
       
        edgeConn[f[i,2]].append(f[i,1])
        edgeConn[f[i,2]].append(f[i,3])
       
        edgeConn[f[i,3]].append(f[i,2])
        edgeConn[f[i,3]].append(f[i,0])


    for n in range(itt):
        q=p.copy()
        for i in range(num_v):
            # % Index to vertices that are edge connected to vertex k
            adj = np.unique(edgeConn[i])
            
            if adj != []:
                # % weight is proportional to distance.
                onesM = np.array([np.ones(np.shape(adj))],np.float32)
                onesM = onesM.transpose()
                #print "onesT\n", onesM
                #print "v[i,:]\n", v[i,:]
                w = np.sqrt(np.sum(((onesM*v[i,:])-q[adj,:])**2,1))
                #print "w\n", w
    
                p[i,:] = v[i,:] + alpha*(1.0/np.sum(w)) * np.sum(np.array([w, w, w]).transpose()*(q[adj,:]-onesM*v[i,:]),axis=0)
                #print "p[i,:]\n", p[i,:]

    return np.array(p)

def move_in_surface(p, ptch, vel):

    # MoveInSurface(p,ptch,vel,varargin)
    # 
    # Moves points p (defined parametrically relative to the Hermite patch 
    # (ptch) which contains the point) on the surface in the direction vel.  
    # 
    #
    # Input
    # p - cell array of structures storing the location of the vertices
    #     relative to the Hermite patches.
    # ptch - cell array of structures storing the Hermite surface description
    #       of each face.
    # vel - array of velocity vectors describing the direction of motion for 
    #       each point.
    # Optional argument 'Move_Nodes_On_Edge' constrains the points on the
    # boundary to remain on the boundary edges.
    # 
    # Output
    # v - updated vertex locations
    # p - cell array of structures storing the updated location of the vertices
    #     relative to the Hermite patches.
    # ptch - cell array of structures storing the Hermite surface description
    #       of each face.
    #
    #
    # Greg Sturgeon
    # December 2010
        
    num_p = len(p)
    
    v = np.zeros([num_p,3],np.float32)
    tau_u = np.zeros([num_p,1])
    tau_w = np.zeros([num_p,1])
    target = np.zeros([num_p,3])
    
    while True:
        for i in range(num_p):    
            target[i,:] = p[i].xyz - vel[i,:]
            tau_u[i], tau_w[i] = pullback_2uv(p[i],target[i,:])
        
        # Rescale vel to traverse no more than half the patch in either
        # direction.
        if any(abs(tau_u)>0.5) or any(abs(tau_w)>0.5):
            vel = vel*0.5
        else:
            break

    
    for i in range(num_p):

        du = tau_u[i]
        dw = tau_w[i]
        
        
        # Check if u or w would be <0 or >1

        if (p[i].u + du) > 1:
            # % Crosses edge 2 3 - subtract 1 for python
            dir = np.array([1, 2])
            p[i] = crossed_u_edge(p[i],ptch,dir,target[i,:])   
        elif (p[i].u + du) < 0:
            # % Crosses edge 1 4 - subtract 1 for python
            dir = np.array([0, 3])
            p[i] = crossed_u_edge(p[i],ptch,dir,target[i,:])
        else:
            p[i].u = p[i].u + du
        
        # Move nodes on an edge back to the edge.
        
        if p[i].isOnBnd:  # removed: "and any(strcmp(varargin,'Move_Nodes_On_Edge'))"
            p[i].w = 0
        
        elif (p[i].w + dw) > 1:
            # Crosses edge 3 4 - subtract 1 for python
            dir = np.array([3, 2])
            p[i] = crossed_w_edge(p[i],ptch,dir,target[i,:]) 
        elif (p[i].w + dw) < 0:
            # Crosses edge 1 2 - subtract 1 for python
            dir = np.array([0, 1])
            p[i] = crossed_w_edge(p[i],ptch,dir,target[i,:])
        else:
            p[i].w = p[i].w + dw
        
        v[i,:] = hermite_point(p[i].Bx,p[i].By,p[i].Bz,p[i].u,p[i].w)
        p[i].xyz = v[i,:]

    return v, p, ptch

def pullback_2uv(p, target):

    vel = target - p.xyz
    epsilon =  0.001
    
    x_u = ((hermite_point(p.Bx,p.By,p.Bz,p.u+epsilon,p.w) - 
            hermite_point(p.Bx,p.By,p.Bz,p.u-epsilon,p.w))/(2*epsilon))
    x_w = ((hermite_point(p.Bx,p.By,p.Bz,p.u,p.w+epsilon) - 
            hermite_point(p.Bx,p.By,p.Bz,p.u,p.w-epsilon))/(2*epsilon))

    MI = np.array([[np.dot(x_u,x_u), np.dot(x_u,x_w)],[np.dot(x_u,x_w), np.dot(x_w,x_w)]])

    S = np.cross(x_u,x_w.transpose())
    
    #~ print "x_u = ", x_u
    #~ print "x_w = ", x_w
    #~ print "MI = ", MI
    #~ print "S = ", S

    x_uu = ((hermite_point(p.Bx,p.By,p.Bz,p.u+epsilon,p.w)
        - 2*hermite_point(p.Bx,p.By,p.Bz,p.u,p.w)) +
        hermite_point(p.Bx,p.By,p.Bz,p.u-epsilon,p.w)/(epsilon**2));

    x_ww = ((hermite_point(p.Bx,p.By,p.Bz,p.u,p.w+epsilon)
        - 2*hermite_point(p.Bx,p.By,p.Bz,p.u,p.w)) +
        hermite_point(p.Bx,p.By,p.Bz,p.u,p.w-epsilon)/(epsilon**2));

    x_uw = ((hermite_point(p.Bx,p.By,p.Bz,p.u+epsilon,p.w+epsilon)
        - hermite_point(p.Bx,p.By,p.Bz,p.u-epsilon,p.w+epsilon))/(2*np.sqrt(2)*epsilon))

    MII_inuw = np.array([[np.dot(-S,x_uu), np.dot(-S,x_uw)],[np.dot(-S,x_uw), np.dot(-S,x_ww)]])
    
    #~ print "MII_inuw = ", MII_inuw
    
    #weirdly, eigh() routine agrees with matlab but eig() does not
    
    V = scipy.linalg.eigh(MII_inuw, MI) # this is the orig
    #V = numpy.linalg.eigh(MII_inuw, MI)
    V_vecs = V[1] #pluck out eigenvectors from tuple
    
    #~ print "V = ", V

    # Principal Coordinates in uv
    f1_uv = V_vecs[:,0].copy().transpose()
    f2_uv = V_vecs[:,1].copy().transpose()
    
    #~ print "f1_uv = ", f1_uv
    #~ print "f2_uv = ", f2_uv
    
    # Principal Coordinates in Tangent Plane
    f1 = hermite_point(p.Bx,p.By,p.Bz,p.u+0.01*f1_uv[0],p.w+0.01*f1_uv[1]) - p.xyz
    f2 = hermite_point(p.Bx,p.By,p.Bz,p.u+0.01*f2_uv[0],p.w+0.01*f2_uv[1]) - p.xyz
    
    # changed from np to scipy
    f1 = f1/scipy.linalg.norm(f1)
    f2 = f2/scipy.linalg.norm(f2)
    
    #~ print "f1 = ", f1
    #~ print "f2 = ", f2
    
    t1 = np.dot(vel,f1)
    t2 = np.dot(vel,f2)
    
    #~ print "t1 = ", t1
    #~ print "t2 = ", t2
    
    # changed from np to scipy
    MI_inv = scipy.linalg.inv(MI)
    
    #~ print "MI_inv = ", MI_inv
    
    f1_uv.shape = (2,1)
    f2_uv.shape = (2,1)
    
    # changed from np to scipy
    a = t1 * np.matrix(MI_inv)*(scipy.linalg.pinv(f1_uv).transpose())
    b = t2 * np.matrix(MI_inv)*(scipy.linalg.pinv(f2_uv).transpose())
    
    f1_uv.shape = (2,)
    f2_uv.shape = (2,)
    
    #~ print "a = ", a
    #~ print "b = ", b
    
    tau_u = a[0] + b[0]
    tau_w = a[1] + b[1]
    
    du = tau_u
    dw = tau_w
    
    #~ print "du = ", du
    #~ print "dw = ", dw
    #~ print " "
    
    return du, dw

def crossed_u_edge(p,ptch,dir,target):

    crossedEdge = np.array([p.f[dir[0]], p.f[dir[1]]])
    
    newPtch, Bx, By, Bz, face = find_new_patch(ptch,p,crossedEdge,dir)

    if newPtch == None:
        if (dir[0] == 1 and dir[1] == 2):
            p.u = 1
        elif (dir[0] == 0 and dir[1] == 3):
            p.u = 0
    
    else:
        p.ptch_num = newPtch
        p.Bx = Bx
        p.By = By
        p.Bz = Bz
        p.f = face

        if (dir[0] == 1 and dir[1] == 2):
            p.u = 0
        elif (dir[0] == 0 and dir[1] == 3):     
            p.u = 1
            
        p.xyz = hermite_point(p.Bx,p.By,p.Bz,p.u,p.w)
        
        
        du, dw = pullback_2uv(p,target)

        if abs(du) > 1:
            print("du -too large might cause the point to traverse multiple patches")

        if (dir[0] == 1 and dir[1] == 2):
            if du > 1:
                dir = np.array([1,2])
                p.u = 1
                p = crossed_u_edge(p,ptch,dir,target)
            elif (p.u + du) < 0:
                p.u = 0
            else:
                p.u = du
                
        elif (dir[0] == 0 and dir[1] == 3):
            if du < -1:
                dir = np.array([0, 3])
                p.u = 0
                p = crossed_u_edge(p,ptch,dir,target)
            elif (p.u + du) > 1:
                p.u = 1
            else:
                p.u = 1 + du

        p.xyz = hermite_point(p.Bx,p.By,p.Bz,p.u,p.w)

    return p

def crossed_w_edge(p, ptch, dir, target):

    crossedEdge = np.array([p.f[dir[0]], p.f[dir[1]]])
    
    newPtch, Bx, By, Bz, face = find_new_patch(ptch,p,crossedEdge,dir)

    if newPtch == None:
        if (dir[0] == 3 and dir[1] == 2):
            p.w = 1
        elif (dir[0] == 0 and dir[1] == 1):
            p.w = 0   

    else:
        p.ptch_num = newPtch
        p.Bx = Bx
        p.By = By
        p.Bz = Bz
        p.f = face

        if (dir[0] == 3 and dir[1] == 2):
            p.w = 0
        elif (dir[0] == 0 and dir[1] == 1):
            p.w = 1
        
        p.xyz = hermite_point(p.Bx,p.By,p.Bz,p.u,p.w)
        
        du, dw = pullback_2uv(p,target)

        if abs(dw) > 1:
            print("dw -too large might cause the point to traverse multiple patches")
        
        if (dir[0] == 3 and dir[1] == 2):
            if dw > 1:
                dir = np.array([3, 2])
                p.w = 1
                p = crossed_w_edge(p,ptch,dir,target)
            elif (p.w + dw) < 0:
                p.w = 0
            else:
                p.w = dw
                
        elif (dir[0] == 0 and dir[1] == 1):
            if dw < -1:
                dir = np.array([0, 1])
                p.w = 0
                p = crossed_w_edge(p,ptch,dir,target)
            elif (p.w + dw) > 1:
                p.w = 1
            else:
                p.w = 1 + dw

        p.xyz = hermite_point(p.Bx,p.By,p.Bz,p.u,p.w)
        
    return p

def find_new_patch(ptch, p, crossedEdge, dir):

    newPtch = None
    neighbors = ptch[p.ptch_num].neighbor
    for k in range(len(neighbors)):
        Q = True
        for num in crossedEdge:
            if num not in ptch[neighbors[k]].f:
                Q = False
        if Q:
            newPtch = neighbors[k]
            
    if newPtch == None:
        Bx = p.Bx
        By = p.By
        Bz = p.Bz
        face = p.f
    else:
        # Rotate the new patch
        Bx, By, Bz, face = rotate_patch2(ptch[newPtch].f, 
                                         crossedEdge, 
                                         dir, 
                                         ptch[newPtch].Bx, 
                                         ptch[newPtch].By, 
                                         ptch[newPtch].Bz)

    return newPtch, Bx, By, Bz, face

def rotate_patch2(face, crossedEdge, dir, Bx, By, Bz):

    # Rotate so that the paramater along the crossedEdge matches directions 
    # across the patches.  ie when the crossedEdge is 2 3 then the new patch 
    # should be rotated so edge 1 4 matches

    if (dir[0] == 1 and dir[1] == 2):
        rotate = 1
        while (not(face[0] == crossedEdge[0] and face[3] == crossedEdge[1]) and  rotate <=5):
            face = np.array([face[1], face[2], face[3], face[0]])
            rotate = rotate +1
            Bx = rotate_parameterization(Bx)
            By = rotate_parameterization(By)
            Bz = rotate_parameterization(Bz)
            
        if rotate == 6:
            print("ERROR: Rotated 360 degrees")
            
    elif (dir[0] == 0 and dir[1] == 3):
        rotate = 1
        while (not(face[1] == crossedEdge[0] and face[2] == crossedEdge[1]) and rotate <=5):
            face = np.array([face[1], face[2], face[3], face[0]])
            rotate = rotate +1
            Bx = rotate_parameterization(Bx)
            By = rotate_parameterization(By)
            Bz = rotate_parameterization(Bz)

        if rotate == 6:
            print("ERROR: Rotated 360 degrees")
            
    elif (dir[0] == 3 and dir[1] == 2):
        rotate = 1
        while (not(face[0] == crossedEdge[0] and face[1] == crossedEdge[1]) and  rotate <=5):
            face = np.array([face[1], face[2], face[3], face[0]])
            rotate = rotate +1
            Bx = rotate_parameterization(Bx)
            By = rotate_parameterization(By)
            Bz = rotate_parameterization(Bz)
            
        if rotate == 6:
            print("ERROR: Rotated 360 degrees")
      
    elif (dir[0] == 0 and dir[1] == 1):
        rotate = 1
        while (not(face[3] == crossedEdge[0] and face[2] == crossedEdge[1]) and rotate <=5):
            face = np.array([face[1], face[2], face[3], face[0]])
            rotate = rotate +1
            Bx = rotate_parameterization(Bx)
            By = rotate_parameterization(By)
            Bz = rotate_parameterization(Bz)

        if rotate == 6:
            print("ERROR: Rotated 360 degrees")

    return Bx, By, Bz, face

def hermite_point(Bx, By, Bz, u, w):

    Mf = np.matrix([[2., -2., 1., 1.,], [-3., 3., -2., -1.], [0., 0., 1., 0.], [1., 0., 0., 0.]])

    u_vector = np.array([u**3, u**2, u, 1.])
    w_vector = np.array([[w**3], [w**2], [w], [1.]])
    
    px = float(u_vector * Mf * np.matrix(Bx) * Mf.T * w_vector)
    py = float(u_vector * Mf * np.matrix(By) * Mf.T * w_vector)
    pz = float(u_vector * Mf * np.matrix(Bz) * Mf.T * w_vector)

    p = np.array([px, py, pz])
    #~ print "xyz = ", p
    
    return p

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
