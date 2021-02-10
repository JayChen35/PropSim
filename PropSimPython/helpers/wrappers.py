# Helper Methods for PropSimPython
# Project Caelus, Aphlex 1C Engine
# 02 February, 2021

import numpy as np
import math


def fast_interp_1(X : list, V : list, x_q: int or float) -> np.ndarray:
    n_x = len(X)
    dx = (X[-1]-X[1])/(n_x-1) #number
    
    i_x = 1 + (x_q - X[1])/dx #number
    
    i_x_lo, i_x_hi = generate_samples(i_x, n_x)
    
    V_lo = np.zeros(x_q)
    V_hi = np.zeros(x_q)
    for ii in range(1, x_q):
        V_lo[ii] = V[i_x_lo]
        V_hi[ii] = V[i_x_hi]
    
    V_q = np.multiply((i_x_hi - i_x), V_lo) + np.multiply((i_x - i_x_lo), V_hi)
    
    return V_q

def generate_samples(i_x, n_x):
    i_x_lo = math.floor(i_x)
    i_x_hi = math.ceil(i_x)
    
    #Check for integers
    if i_x_hi == i_x_lo:
        i_x_hi += 1

    #Low range check
    if i_x_lo < 1:
        i_x_hi = 2
        i_x_lo = 1

    #High range check
    if i_x_hi > n_x:
        i_x_lo = n_x-1
        i_x_hi = n_x
        
    return i_x_lo, i_x_hi

def fast_interp_2():
    raise NotImplementedError


def fast_interp_3():
    raise NotImplementedError