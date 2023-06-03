"Claudiu Acsinte 233207"

import numpy as np
from numpy import linalg as la
from matplotlib import pyplot as plt
from time import time



def f_1d (x):
    return (x-1)**3

def g_1d (x):
    return (x**2)/10

def f_2d (x,y):
    return x**2-y**2

def g_2d (x,y):
    return (x**2+y**2)/10


"""

    1-d version:
    
"""

# This Gauss-Seidel solves only for internal nodes, leaving the boundaries 
# unchanged, so we can use it for both uh and eh
def gs_step_1d(uh, fh):
    n = len(uh)-1
    h = 4/n
    uh_ = np.copy(uh)
    for i in range(1,n):
        uh[i] = (fh[i]*(h**2)+uh[i-1]+uh[i+1])/(2+h**2)
    sres = la.norm(uh-uh_,np.inf)
    return sres

# From fine to coarse, this is the version with full weighting, 
# It creates homogeneous boundaries, so it must be used for residuals
def restr(vh):
    n = len(vh) - 1
    v2h = np.zeros(n//2+1)
    for i in range(1,n//2):
        v2h[i] = (vh[2*i-1]+2*vh[2*i]+vh[2*i+1])/4
    return v2h

# From fine to coarse, this is the injection version, we will use it later on
# for the full multigrid method.
def restr2(vh):
    n = len(vh)-1
    v2h = np.zeros(n//2+1)
    for i in range(0,n//2+1):
        v2h[i] = vh[2*i]
    return v2h

## From coarse to fine
def prol(v2h):
    n = 2*(len(v2h)-1)
    vh = np.zeros(n+1)
    
    for i in range(0,n//2):
        vh[2*i] = v2h[i]
        vh[2*i+1] = (v2h[i]+v2h[i+1])/2
    
    # manteining the original boundaries condition
    vh[0] = v2h[0]
    vh[n] = v2h[n//2]
    return vh

def twoGrids_step_1d(uh,fh):
    n = len(uh)-1
    h = 4/n
    
    gs_step_1d(uh,fh) #Pre smoothing
    
    resh = np.zeros_like(uh) 
    for i in range(1,n):
        resh[i] = fh [i] - (-(uh[i-1]/(h**2))-(uh[i+1]/(h**2))+((2/(h**2))+1)*uh[i])    #print(f"residual before restriction: resh = {resh}")
    res2h = restr(resh) # restriction on residuals
    
    e2h = np.zeros_like(res2h)
    
    for i in range(0,5): # we calculate e_2h by performing 5 GS steps, as requested
        gs_step_1d(e2h,res2h)
        
    eh = prol(e2h) # we prolungate the error to the fine grid
    
    #uh = uh + Prol(e2h) # correction on fine grid
    for i in range(1,n):
        uh[i] += eh[i]
        
    sres = gs_step_1d(uh,fh) # post smooting

    return sres
    
     
def v_cycle_step_1d(uh,fh,alpha1,alpha2):
    
    n = len(uh) - 1
    h = 4/n
    
    if n==2:
        uh[1] = (uh[0]+(h**2)*fh[1]+uh[2])/(h**2+2)
        return 0
    
    else:
        for _ in range(0,alpha1):
            gs_step_1d(uh,fh)
            
        resh = np.zeros(n+1)
        for i in range(1,n):
            resh[i] = fh[i] - (-uh[i-1]/h**2-uh[i+1]/h**2+(2/h**2+1)*uh[i])
        f2h = restr(resh)
        u2h = np.zeros_like(f2h)
        
        v_cycle_step_1d(u2h,f2h,alpha1,alpha2)
        
        u2h = prol(u2h)
        for i in range(1,n):
            uh[i] += u2h[i]
            
        for _ in range(0,alpha2):
            sres = gs_step_1d(uh,fh)
    return sres


def full_mg_1d(uh, fh, alpha1, alpha2, nu):
    
    n = len(uh) - 1
    h = 4/n
    
    if n==2:
        
        uh[1] = (fh[0]+(h**2)*fh[1]+fh[2])/(h**2+2)
        return 0
    
    else:
        f2h = restr2(fh)
        u2h = np.zeros_like(f2h) 
        u2h[0] = uh[0]
        u2h[n//2] = uh[n]

        full_mg_1d(u2h, f2h, alpha1, alpha2, nu)

        uh_ = prol(u2h)
        for i in range(1,n):
            uh[i] = uh_[i]
            
        for _ in range(0,nu):
           psres =  v_cycle_step_1d(uh,fh,alpha1,alpha2)
    
    return psres 

"""
    
    2-d version
    
"""

def gs_step_2d(uh, fh):
    n = uh.shape[0]-1
    h = 4/n
    maxsofar = 0
    uh_ = np.copy(uh)
    for i in range(1,n):
        for j in range(1,n):
            uh[i,j] = (fh[i,j]*(h**2)+uh[i-1,j]+uh[i+1,j]+uh[i,j-1]+uh[i,j+1])/(4+h**2)
            maxsofar = max(maxsofar,abs(uh[i,j]-uh_[i,j]))
    return maxsofar

# From fine to coarse, this is the version with full weighting, 
# Output has homogeneous boundaries, so it must be used for residuals and eh
def restr_2d(vh):
    n = vh.shape[0] - 1
    v2h = np.zeros((n//2+1,n//2+1))
    for i in range(1,n//2):
        for j in range(1,n//2):
            v2h[i,j] = (vh[2*i-1,2*j-1]+vh[2*i-1,2*j+1]+vh[2*i+1,2*j-1]+vh[2*i+1,2*j+1]+2*(vh[2*i-1,2*j]+vh[2*i+1,2*j]+vh[2*i,2*j-1]+vh[2*i,2*j+1])+4*vh[2*i,2*j])/16
    return v2h

# From fine to coarse, this is the injection version, we will use it later on
# for the full multigrid method.
def restr2_2d(vh):
    n = vh.shape[0] - 1
    v2h = np.zeros((n//2+1,n//2+1))
    for i in range(0,n//2+1):
        for j in range(0,n//2+1):
            v2h[i,j] = vh[2*i,2*j]
    return v2h

## From coarse to fine
def prol_2d(v2h):
    n = 2*(v2h.shape[0]-1)
    vh = np.zeros((n+1,n+1))
    for i in range(0,n//2):
        for j in range(0,n//2):
            vh[2*i,2*j] = v2h[i,j]
            vh[2*i+1,2*j] = (v2h[i,j]+v2h[i+1,j])/2
            vh[2*i,2*j+1] = (v2h[i,j]+v2h[i,j+1])/2
            vh[2*i+1,2*j+1] = (v2h[i,j]+v2h[i+1,j]+v2h[i,j+1]+v2h[i+1,j+1])/4
        vh[0,i] = v2h[0,i]
        vh[n,i] = v2h[n//2,i]
        vh[i,0] = v2h[i,0]
        vh[i,n] = v2h[i,n//2]
    return vh

def twoGrids_step_2d(uh,fh):
    n = uh.shape[0]-1
    h = 4/n
    
    gs_step_2d(uh,fh) #Pre smoothing
    
    resh = np.zeros_like(uh) 
    for i in range(1,n):
        for j in range(1,n):
            resh[i,j] = fh[i,j] - (-(uh[i-1,j]/(h**2))-(uh[i+1,j]/(h**2))-(uh[i,j-1]/h**2)-(uh[i,j+1]/h**2)+((4/(h**2))+1)*uh[i,j])
    
    
    res2h = restr_2d(resh) # restriction on residuals
    
    e2h = np.zeros_like(res2h)
    
    for _ in range(0,5): # we calculate e_2h by performing 5 GS steps, as requested
        gs_step_2d(e2h,res2h)
        
    eh = prol_2d(e2h) # we prolungate the error to the fine grid
    
    #uh = uh + Prol(e2h) # correction on fine grid
    for i in range(1,n):
        for j in range(1,n):
             uh[i,j] += eh[i,j]
    
    sres = gs_step_2d(uh,fh) # post smooting
    

    return sres
    
     
def v_cycle_step_2d(uh,fh,alpha1,alpha2):
    
    n = uh.shape[0] - 1
    h = 4/n
    
    if n==2:
        uh[1,1] = (fh[0,1] + fh[1,0] + fh[2,1] + fh[1,2] + (h**2)*fh[1,1])/(h**2+4)
        return 0
    
    else:
        for _ in range(0,alpha1):
            gs_step_2d(uh,fh)
            
        resh = np.zeros((n+1,n+1))
        for i in range(1,n):
            for j in range(1,n):
                resh[i,j] = fh[i,j] - (-(uh[i-1,j]/(h**2))-(uh[i+1,j]/(h**2))-(uh[i,j-1]/h**2)-(uh[i,j+1]/h**2)+((4/(h**2))+1)*uh[i,j])
        f2h = restr_2d(resh)
        u2h = np.zeros_like(f2h)
        
        v_cycle_step_2d(u2h,f2h,alpha1,alpha2)
        
        u2h = prol_2d(u2h)
        for i in range(1,n):
            for j in range(1,n):
                uh[i,j] += u2h[i,j]
                
        for _ in range(0,alpha2):
            sres = gs_step_2d(uh,fh)
    return sres


def full_mg_2d(uh, fh, alpha1, alpha2, nu):
    
    n = uh.shape[0] - 1
    h = 4/n
    
    if n==2:
        uh[1,1] = (fh[0,1] + fh[1,0] + fh[2,1] + fh[1,2] + (h**2)*fh[1,1])/(h**2+4)
        return 0
    
    else:
        f2h = restr2_2d(fh)
        u2h = np.zeros_like(f2h) 
        for i in range(0,n//2+1):
             u2h[0,i] = uh[0,2*i]
             u2h[n//2,i] = uh[n,2*i]
             u2h[i,0] = uh[2*i,0]
             u2h[i,n//2] = uh[2*i,n]

        full_mg_2d(u2h, f2h, alpha1, alpha2, nu)

        uh_ = prol_2d(u2h)
        for i in range(1,n):
            for j in range(1,n):
                uh[i,j] = uh_[i,j]
            
        for _ in range(0,nu):
           psres =  v_cycle_step_2d(uh,fh,alpha1,alpha2)
    
    return psres 

