"""
Improve element quality based on Vartziotis & Wipper which performs a geometric 
transformation on the dual to the hexahedra.

The normal vectors from the centroid of the octahedral faces are used to create a
new hexahedron.   Where the octahedron is the dual to the original hex.  The new
hex is then scaled based on the ratio of the principal axes of the new and 
original hexahedron.

Relaxation is preformed on the anisotropically scaled hex to regularize the 
lengths of its edges.


July 19, 2012
Greg Sturgeon
"""

import numpy

def hex_dual_non_iso(e,n):        
    if e.min()==1:
        e=e-1
    
    num_n = numpy.size(n,0)
    num_e = numpy.size(e,0)
    
#    sigma_min = .5
#    sigma_max = .6
#    eta = .5
    rho = .75
    tr = .25
    
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
    
    
    for k in range(8):
        pp[:,k,:] = c[:,k,:] + nk[:,k,:]

    
    opp = numpy.zeros([num_e,6,3])
    for k in range(6):
        opp[:,k,:] = numpy.sum(.25*(pp[:,fh[:,k],:]),1)
    
    
    Horig = n[e,:]
    H = pp
    
    C = .125*numpy.sum(pp,1)
    Hs = numpy.zeros(H.shape)
    ie = 0
    s = numpy.zeros([num_e,3,3])
    for ie in range(num_e):
        lpa_o = numpy.array([numpy.sqrt(numpy.sum((o[ie,0,:] - o[ie,5,:])**2)),numpy.sqrt(numpy.sum((o[ie,1,:] - o[ie,3,:])**2)),numpy.sqrt(numpy.sum((o[ie,2,:] - o[ie,4,:])**2))])
        lpa_opp = numpy.array([numpy.sqrt(numpy.sum((opp[ie,0,:] - opp[ie,5,:])**2)),numpy.sqrt(numpy.sum((opp[ie,1,:] - opp[ie,3,:])**2)),numpy.sqrt(numpy.sum((opp[ie,2,:] - opp[ie,4,:])**2))])
        
        s = lpa_o / lpa_opp
        pa_opp = numpy.array([opp[ie,0,:] - opp[ie,5,:],opp[ie,1,:] - opp[ie,3,:],opp[ie,2,:] - opp[ie,4,:]])
        
        # Normalize pa_opp
        pa_opp = pa_opp/numpy.tile(lpa_opp.reshape(3,1),(1,3))
        # Need to Orthonormalize pa_opp
        lpa = lpa_opp.copy()
        
        order = numpy.zeros((3,1),'int')
        order[0] = numpy.argmax(lpa)
        lpa[order[0]] = 0
        order[1] = numpy.argmax(lpa)
        lpa[order[1]] = 0
        order[2] = numpy.argmax(lpa)
        
        pa_on = numpy.zeros((3,3))
        pa_on[0,:] = numpy.cross(pa_opp[order[1],:],pa_opp[order[2],:])
        pa_on[1,:] = numpy.cross(pa_on[0,:],pa_opp[order[2],:])
        pa_on[2,:] = numpy.cross(pa_on[0,:],pa_on[1,:])
        
        # Check the negative signs if the order is changed
        for k in range(3):
            pa_on[k,:] = numpy.sign(numpy.dot(pa_opp[order[k],:],pa_on[k,:]))*pa_on[k,:]
        
        pa_on = normalize_v3(pa_on)
        pa_on = pa_on.T
            
        
        s[numpy.logical_and(abs(s)<0.5,s<0)] = -0.5
        s[numpy.logical_and(abs(s)<0.5,s>0)] = 0.5
        s[numpy.nonzero(s<-2.0)] = -2.0
        s[numpy.nonzero(s>2.0)] = 2.0
        
        s = s[order].reshape(3)
        
        Hs[ie,:,:] = numpy.dot(pa_on, numpy.dot(numpy.linalg.inv(pa_on),(H[ie,:,:]- numpy.tile(C[ie,:],(8,1))).T) * numpy.tile(s,(8,1)).T).T + numpy.tile(C[ie,:],(8,1))
    
    Hr = (1-rho)*Horig + rho*(Hs)
    # Equal weighting
    ppij = numpy.zeros((num_n,3))
    wij = numpy.zeros((num_n,1))
    for ie in range(num_e):
        for k in range(8):
            ppij[e[ie,k],:] = ppij[e[ie,k],:] + Hr[ie,k,:]
            wij[e[ie,k],:] = wij[e[ie,k],:] + 1.0
    
    wij = numpy.concatenate((wij,wij,wij),1)
    p_new = ppij/wij
    
    # Trust Region
    m = p_new - n
    m_mag = numpy.sqrt(numpy.sum(m**2,1))
    m_norm = normalize_v3(m)
    
    m_mag[m_mag>tr] = tr
#    m_mag = tr * m_mag / m_mag.max()
    m = m_mag.reshape(m.shape[0],1)*m_norm
    
    
    p_new = n + m
    
#    # Invalid element handling
#    p_new = InvalidElementHandling(e,p_new,n)
    
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

def invalid_element_isolation(e,p_new,n):
    minQ_Original = FindInvalidElements(e,n)
    
    print("%s Invalid Elements Originally" % str(len(numpy.nonzero(minQ_Original<=0)[0])))
    # Invalid element handling
    wereInvalidElements = 1  
    while wereInvalidElements:
        minQ = FindInvalidElements(e,p_new)
    
        print("%s Invalid Elements" % str(len(numpy.nonzero(minQ<=0)[0])))
        
        
        invalid = numpy.logical_and(minQ<=0,minQ < minQ_Original)
        invalid = numpy.nonzero(invalid)[0]
        
        
        if numpy.size(invalid)>0:
            p_new[e[invalid,:],:] = n[e[invalid,:],:]
            print("Fixing %s Invalid Elements" % str(len(invalid)))
        else:
            wereInvalidElements = 0    
    return p_new


def invalid_element_handling(e,p_new,n):
    minQ_Original = FindInvalidElements(e,n)
    
    print("%s Invalid Elements Originally" % str(len(numpy.nonzero(minQ_Original<=0)[0])))
    # Invalid element handling
    wereInvalidElements = 1  
    while wereInvalidElements:
        minQ = FindInvalidElements(e,p_new)
    
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
    q = QMetric(e,n)
    
#    q = numpy.sum(q,1)
#    q = q.reshape((q.shape[0],1))
    
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
