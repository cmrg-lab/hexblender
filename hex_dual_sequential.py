"""
Improve element quality based on Vartziotis & Wipper which performs a geometric 
transformation on the dual to the hexahedra.  
Improves the shape of the element that contains the worst tetrahedron.
The normal vectors from the centroid of the octahedral faces are used to create a
new hexahedron.   Where the octahedron is the dual to the original hex.  The new
hex is then scaled based on the ratio of the principal axes of the new and 
original hexahedron.
Relaxation is preformed on the anisotropically scaled hex to regularize the 
lengths of its edges.

The algorithm stops upon reaching a maximum number of itterations or when none 
of the elements have poorly shaped tetahedrons.


July 13, 2012
Greg Sturgeon
"""

import numpy

from hexblender.regularize_elements import normalize_v3


def hex_dual_sequential(e,n,iters=10):
    if e.min()==1:
        e=e-1
    
    num_e = numpy.size(e,0)
    
    tr = 0.25 # Trust region
    sigma = 0.5
    
    fh = numpy.array([[0,0,1,2,0,4],[1,4,5,6,3,7],[2,5,6,7,7,6],[3,1,2,3,4,5]])
    
    fo = numpy.array([[0,0,0,0,5,5,5,5],[1,2,3,4,4,1,2,3],[4,1,2,3,1,2,3,4]])
    
    noPoorlyShaped = 0
    k_itt=0
    #while k_itt<=num_e:
    while k_itt<=iters:
        print("k_itt = %s" % k_itt)
        qM = QMetricNormalize(e,n)
        q = qM.min(axis=1)
        j = q.argmin()
        
        print("%s \t %s \t" % (str(j), str(q.min()), str(sum(q<=0))))
        if sum(q<=0) == 0:
            noPoorlyShaped = 1
            break        
        
        # Apply Geometric Transform
        c = numpy.zeros([8,3])
        o = numpy.zeros([6,3])
        nk = numpy.zeros([8,3])
        pp = numpy.zeros([8,3])
        
        for k in range(6):
            o[k,:] = numpy.sum(.25*(n[e[j,fh[:,k]],:]),0)
        
        for k in range(8):
            c[k,:] = numpy.sum(1.0/3.0*(o[fo[:,k],:]),0)
            nk[k,:] = numpy.cross(o[fo[1,k],:] - o[fo[0,k],:],o[fo[2,k],:] - o[fo[0,k],:])
            pp[k,:] = c[k,:] + sigma/(numpy.sqrt(numpy.sqrt(sum(nk[k,:]**2))))*nk[k,:]
        
        
        opp = numpy.zeros([6,3])
        for k in range(6):
            opp[k,:] = numpy.sum(.25*(pp[fh[:,k],:]),0)
        
        H = pp
        # Centroid of new Hex H
        C = .125*numpy.sum(pp,0)
        
        # Anisotropic scaling
        # lpa - length of principal axes of the octahedron
        lpa_o = numpy.array([numpy.sqrt(numpy.sum((o[0,:] - o[5,:])**2)),numpy.sqrt(numpy.sum((o[1,:] - o[3,:])**2)),numpy.sqrt(numpy.sum((o[2,:] - o[4,:])**2))])
        lpa_opp = numpy.array([numpy.sqrt(numpy.sum((opp[0,:] - opp[5,:])**2)),numpy.sqrt(numpy.sum((opp[1,:] - opp[3,:])**2)),numpy.sqrt(numpy.sum((opp[2,:] - opp[4,:])**2))])
        
        s = lpa_o / lpa_opp
        # Principal axes of the new octahedron opp
        pa_opp = numpy.array([opp[0,:] - opp[5,:],opp[1,:] - opp[3,:],opp[2,:] - opp[4,:]])
        
        # Orthonormalize pa_opp
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
        pa_on = pa_on.T #Principal Axes OrthoNormalized
        
        # Ensure that the scaling is within a reasonable range
        s[numpy.logical_and(abs(s)<0.5,s<0)] = -0.5
        s[numpy.logical_and(abs(s)<0.5,s>0)] = 0.5
        s[numpy.nonzero(s<-2.0)] = -2.0
        s[numpy.nonzero(s>2.0)] = 2.0
        
        # Order the scaling to agree with the order of the principal axes
        s = s[order].reshape(3)
        
        Hs = numpy.dot(pa_on, numpy.dot(numpy.linalg.inv(pa_on),(H[:,:]- numpy.tile(C,(8,1))).T) * numpy.tile(s,(8,1)).T).T + numpy.tile(C,(8,1))
        
        Hr = Hs
        
        m = Hr - n[e[j],:]
        m_mag = numpy.sqrt(numpy.sum(m**2,1))
        m_norm = normalize_v3(m)
        
        m = (tr*m_mag/m_mag.max()).reshape(8,1)*m_norm
        
        
        n[e[j],:] = n[e[j],:] + m
        
        k_itt=k_itt+1
    return list([n, noPoorlyShaped])
    
def q_metric_normalize(e,n):
    if e.min()==1:
        e=e-1
    
    num_e = numpy.size(e,0)
    
    eh = numpy.array([[0,1,2,3,4,5,6,7],[3,0,1,2,7,4,5,6],[4,5,6,7,5,6,7,4],[1,2,3,0,0,1,2,3]])
    q = numpy.zeros((num_e,8))
    for k in range(8):
        d0 = normalize_v3(n[e[:,eh[1,k]],:] - n[e[:,eh[0,k]],:])
        d1 = normalize_v3(n[e[:,eh[2,k]],:] - n[e[:,eh[0,k]],:])
        d2 = normalize_v3(n[e[:,eh[3,k]],:] - n[e[:,eh[0,k]],:])
        d = numpy.hstack((d0,d1,d2))
        
        detD = d[:,0]*d[:,4]*d[:,8] + d[:,1]*d[:,5]*d[:,6] + d[:,2]*d[:,3]*d[:,7] - d[:,2]*d[:,4]*d[:,6] - d[:,1]*d[:,3]*d[:,8] - d[:,0]*d[:,5]*d[:,7]
        
        tr = numpy.sum(d**2,axis=1)
        
        i_neg = numpy.nonzero(detD<0)[0]
        i_z = numpy.nonzero(detD==0)[0]
        i_pos = numpy.nonzero(detD>0)[0]
        q[i_neg,k] = -.125*(3.0*(abs(detD[i_neg]) )**(2.0/3.0) / tr[i_neg])
        q[i_z,k] = 0
        q[i_pos,k] = .125*(3.0*(abs(detD[i_pos]) )**(2.0/3.0) / tr[i_pos])
    print("in qmetric normalize!")
    return q
