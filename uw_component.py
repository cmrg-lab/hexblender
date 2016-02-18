# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 15:53:07 2011

@author: sturgeon
"""
import numpy

from hexblender.regularize_elements_pullback_functions import hermite_derivative
from hexblender.regularize_elements_pullback_functions import dot_v3
from hexblender.regularize_elements_pullback_functions import normalize_v3
from hexblender.regularize_elements_pullback_functions import magnitude_v3

def uw_component(mf,nExt,p_Bx, p_By, p_Bz, p_u, p_w):
    vec = mf-nExt
    
    [x_u,x_w] = hermite_derivative(p_Bx, p_By, p_Bz, p_u, p_w)
    
    mI = numpy.array([[(x_u*x_u).sum(1),(x_u*x_w).sum(1)],[(x_w*x_u).sum(1),(x_w*x_w).sum(1)]])
    
    a = dot_v3(vec, normalize_v3(x_u.copy()))
    b = dot_v3(vec, normalize_v3(x_w.copy()))
    a = a.reshape((a.shape[0]))
    b = b.reshape((b.shape[0]))
    
    
    tmp = numpy.array([[magnitude_v3(x_u), magnitude_v3(x_u)],
                       [magnitude_v3(x_w), magnitude_v3(x_w)]],'float')
    tmp.shape = (2,2,tmp.shape[2],)
    mItmp = mI/tmp
    
    mItmpinv = numpy.zeros((2,2,x_u.shape[0]),'float')
    # Inverse of MI
    c = 1/((mItmp[0,0,:]*mItmp[1,1,:]) - (mItmp[0,1,:]*mItmp[1,0,:]))
    mItmpinv[0,0,:] = c*mItmp[1,1,:]
    mItmpinv[1,1,:] = c*mItmp[0,0,:]
    mItmpinv[1,0,:] = -c*mItmp[1,0,:]
    mItmpinv[0,1,:] = -c*mItmp[0,1,:]
    
    vTu = mItmpinv[0,0,:]*a + mItmpinv[0,1,:]*b
    vTw = mItmpinv[1,0,:]*a + mItmpinv[1,1,:]*b
    
    vTu[numpy.nonzero(numpy.isnan(vTu))]=0
    vTw[numpy.nonzero(numpy.isnan(vTw))]=0
    
    return [vTu,vTw]
