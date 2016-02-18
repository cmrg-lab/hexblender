"""
hex_dual(e,n)
e - hexahedral element
n - nodal positions

Return p_new improved nodal positions

Improve element quality based on Vartziotis & Wipper which performs a geometric 
transformation on the dual to the hexahedra.  
1) an octahedron is created for each hexahedra 2) the normals of the tri-faces 
of the octahedron are scaled towards having a consistent length 3) a new hex is 
created from the tips of the normals 4) the hexahedra is scaled to have a 
consistent volume with the original hexahedra 5) relaxation 6) weighted 
averaging of each hexahedral vertex to determine the new node positions 
7) invalid element handling where the nodal positions are not updated if they 
create poorly shaped elements

This implementation allows for motion of the exterior nodes and includes a trust
region to provide an upward bound on the motion.
 


January 17, 2012
Greg Sturgeon
"""

# Improve Q

import numpy


def hex_dual(e,n):        
    if e.min()==1:
        e=e-1
    
    num_n = numpy.size(n,0)
    num_e = numpy.size(e,0)
    
    sigma_min = .5
    sigma_max = .6
    rho = .75
    eta = .5
    
    tr = 0.125
    
    fh = numpy.array([[0,0,1,2,0,4],[1,4,5,6,3,7],[2,5,6,7,7,6],[3,1,2,3,4,5]])
    
    fo = numpy.array([[0,0,0,0,5,5,5,5],[1,2,3,4,4,1,2,3],[4,1,2,3,1,2,3,4]])
    
    o = numpy.zeros([num_e,6,3])
    c = numpy.zeros([num_e,8,3])
    nk = numpy.zeros([num_e,8,3])
    pp = numpy.zeros([num_e,8,3])
    k=0
    
    for k in range(6):
        o[:,k,:] = numpy.sum(.25*(n[e[:,fh[:,k]],:]),1)
    
    
    for k in range(8):
        c[:,k,:] = numpy.sum(1.0/3.0*(o[:,fo[:,k],:]),1)
        nk[:,k,:] = numpy.cross(o[:,fo[1,k],:] - o[:,fo[0,k],:],o[:,fo[2,k],:] - o[:,fo[0,k],:])
    
    qM = QMetric(e,n)
    q = numpy.sum(abs(qM),1)
    q[numpy.any(qM<=0)] = (abs(q)).min()
    
    mtH = vol_hex_tets(n[e,:])
    
    sigma = sigma_min + (sigma_max - sigma_min)*(1-q)
    sigma[sigma>sigma_max] = sigma_max
    sigma[sigma<sigma_min] = sigma_min
    
    for k in range(8):
        # scaling for the normal from tet centroid
        sntc = sigma/(numpy.sqrt(numpy.sqrt(numpy.sum(nk[:,k,:]**2,1))))
        sntc = numpy.vstack((sntc,sntc,sntc)).T
        
        pp[:,k,:] = c[:,k,:] + sntc * nk[:,k,:]
    
    mtHp = vol_hex_tets(pp)
    
    L = (mtH/mtHp)**(1.0/3.0)
    L = numpy.vstack((L,L,L)).T
    
    Horig = n[e,:]
    H = pp
    
    C = .125*numpy.sum(pp,1)
    Hs = numpy.zeros(H.shape)
    
    for k in range(8):
       Hs[:,k,:] = C + L*(H[:,k,:]-C) 
    
    
    # Relax
    Hr = (1-rho)*Horig + rho*(Hs)

    ppij = numpy.zeros((num_n,3))
    wij = numpy.zeros((num_n,1))
    for ie in range(num_e):
        for k in range(8):
            ppij[e[ie,k],:] = ppij[e[ie,k],:] + ((1.0-q[ie])**eta) * Hr[ie,k,:]
            wij[e[ie,k],:] = wij[e[ie,k],:] + ((1.0-q[ie])**eta)
    
    wij = numpy.concatenate((wij,wij,wij),1)
    p_new = ppij/wij
    
    m = p_new - n
    m_mag = numpy.sqrt(numpy.sum(m**2,1)).reshape(m.shape[0],1)
    m_norm = NormalizeV3(m)
    
    m_mag[m_mag>tr] = tr
    m = m_mag*m_norm
    
    p_new = n + m
    
    # Invalid element handling    
    p_new = InvalidElementHandling(e,p_new,n)
    
    return p_new
    

def normalize_v3(arr):
    ''' Normalize a numpy array of 3 component vectors shape=(n,3) '''
    lens = numpy.sqrt( arr[:,0]**2 + arr[:,1]**2 + arr[:,2]**2 )
    # Correct for division by zero    
    lens[numpy.nonzero(lens==0)]=1
    arr[:,0] /= lens
    arr[:,1] /= lens
    arr[:,2] /= lens                
    return arr


def invalid_element_handling(e,p_new,n):
    minQ_Original = find_invalid_elements(e,n)
    
    print("%s Invalid Elements Originally" % str(len(numpy.nonzero(minQ_Original<=0)[0])))
    # Invalid element handling
    wereInvalidElements = 1  
    while wereInvalidElements:
        minQ = find_invalid_elements(e,p_new)
    
        print("%s Invalid Elements" % str(len(numpy.nonzero(minQ<=0)[0])))
        
        
        invalid = numpy.logical_and(minQ<=0,minQ < minQ_Original)
        invalid = numpy.nonzero(invalid)[0]
        
        
        if numpy.size(invalid)>0:
            p_new[e[invalid,:],:] = n[e[invalid,:],:]
            print("Fixing %s Invalid Elements" % str(len(invalid)))
        else:
            wereInvalidElements = 0    
    return p_new


def find_invalid_elements(e,n):
    #Based on determinates of the tetrahedrons at the corners of the hexahedrons
    q = q_metric(e,n)
    
    invalid = q.min(axis=1)
    
    # Return the minimum q value for each element
    return invalid



def q_metric(e,n):
    if e.min()==1:
        e=e-1
    
    num_e = numpy.size(e,0)
    
    eh = numpy.array([[0,1,2,3,4,5,6,7],[3,0,1,2,7,4,5,6],[4,5,6,7,5,6,7,4],[1,2,3,0,0,1,2,3]])
    q = numpy.zeros((num_e,8))
    for k in range(8):
        d = numpy.hstack(((n[e[:,eh[1,k]],:] - n[e[:,eh[0,k]],:]),(n[e[:,eh[2,k]],:] - n[e[:,eh[0,k]],:]),(n[e[:,eh[3,k]],:] - n[e[:,eh[0,k]],:])))
        detD = d[:,0]*d[:,4]*d[:,8] + d[:,1]*d[:,5]*d[:,6] + d[:,2]*d[:,3]*d[:,7] - d[:,2]*d[:,4]*d[:,6] - d[:,1]*d[:,3]*d[:,8] - d[:,0]*d[:,5]*d[:,7]
        
        tr = numpy.sum(d**2,axis=1)
        
        i_neg = numpy.nonzero(detD<0)[0]
        i_z = numpy.nonzero(detD==0)[0]
        i_pos = numpy.nonzero(detD>0)[0]
        q[i_neg,k] = -.125*(3.0*(abs(detD[i_neg]) )**(2.0/3.0) / tr[i_neg])
        q[i_z,k] = 0
        q[i_pos,k] = .125*(3.0*(abs(detD[i_pos]) )**(2.0/3.0) / tr[i_pos])
    
    return q

def vol_hex_tets(h):
    
    num_e = numpy.size(h,0)
    
    eh = numpy.array([[0,1,2,3,4,5,6,7],[3,0,1,2,7,4,5,6],[4,5,6,7,5,6,7,4],[1,2,3,0,0,1,2,3]])
    mtH = numpy.zeros((num_e,8))
    for k in range(8):
        d = numpy.hstack(((h[:,eh[1,k],:] - h[:,eh[0,k],:]),(h[:,eh[2,k],:] - h[:,eh[0,k],:]),(h[:,eh[3,k],:] - h[:,eh[0,k],:])))
        detD = d[:,0]*d[:,4]*d[:,8] + d[:,1]*d[:,5]*d[:,6] + d[:,2]*d[:,3]*d[:,7] - d[:,2]*d[:,4]*d[:,6] - d[:,1]*d[:,3]*d[:,8] - d[:,0]*d[:,5]*d[:,7]
        
        mtH[:,k] = 1.0/(8.0*6.0)*abs(detD)
    
    mtH = numpy.sum(mtH,1)
    return mtH
