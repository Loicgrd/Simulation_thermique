"================================Fonction matrices système linéaire=================================="

import numpy as np
from numpy.linalg import *
import matplotlib.pyplot as plt
from scipy.linalg import lu_factor
from scipy.linalg import lu_solve
import math as math
import scipy.linalg as sp


"===============Monocouche============"

"Matrice A"
def matriceA_mono(lbda, rho, cp, epaisseur, dx ,dt, h_ext,h_int):
    
    alpha = lbda / (rho * cp)
    a = alpha * dt / (dx ** 2)
    Nx = int(epaisseur / dx)
    A = np.zeros((Nx,Nx))
    
    
    for i in range(1,Nx-1):
        A[i,i] = a + 1
        A[i,i+1] = - a/2
        A[i,i-1] = - a/2
        
        
        
    #Coefficients convectifs
    a_ext = lbda / (h_ext * dx)
    a_int = lbda / (h_int * dx)
     
    A[0,0] = 1 + a_ext
    A[0,1] = - a_ext
    A[Nx-1,Nx-1] = 1 + a_int 
    A[Nx-1,Nx-2] = - a_int
        
    return A


"Matrice B"
def matriceB_mono(lbda, rho, cp, epaisseur, dx, dt):
    alpha = lbda / (rho * cp)
    a = alpha * dt / (dx ** 2)
    Nx = int(epaisseur / dx)
    A = np.zeros((Nx,Nx))
    
    # A[Nx-1,Nx-1] = 1
    # A[0,0] = 1
    
    for i in range(1,Nx-1):
        A[i,i] = 1 - a
        A[i,i+1] = a/2
        A[i,i-1] = a/2
    


    
    return A




"===========Double couche==========="
"Matrice A"
def matriceA_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, h_ext, h_int ):
    alpha1 = lbda1 / (rho1 * cp1)
    a1 = alpha1 * dt / (dx1 ** 2)
    Nx1 = int(epaisseur1 / dx1)
    
    alpha2 = lbda2 / (rho2 * cp2)
    a2 = alpha2 * dt / (dx2 ** 2)
    Nx2 = int(epaisseur2 / dx2)
    
    Nx12 = Nx1 + Nx2
    
    A = np.zeros((Nx12,Nx12))

    a12 = lbda1 / (lbda2 * dx1 / dx2 + lbda1)
    a21 = lbda2 / (lbda1 * dx2 / dx1 + lbda2)
     
    
    #Boucle première couche
    for i in range(1,Nx1):
        A[i,i] = a1 + 1
        A[i,i+1] = - a1/2
        A[i,i-1] = - a1/2
    
    #Boucle deuxieme couche
    for i in range(Nx1,Nx12-1):
        A[i,i] = a2 + 1
        A[i,i+1] = - a2/2
        A[i,i-1] = - a2/2
    
    #Point en contact couche 1 et 2
    A[Nx1-1,Nx1-2] = - a12
    A[Nx1-1,Nx1] = - a21
    A[Nx1-1,Nx1-1] = 1    
    

    
    #Coefficients convectifs
    a_ext = lbda1 / (h_ext * dx1)
    a_int = lbda2 / (h_int * dx2)
     
    A[0,0] = 1 + a_ext
    A[0,1] = - a_ext
    A[Nx12-1,Nx12-1] = 1 + a_int 
    A[Nx12-1,Nx12-2] = - a_int
    

    return A




"Matrice B"
def matriceB_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt):
    alpha1 = lbda1 / (rho1 * cp1)
    a1 = alpha1 * dt / (dx1 ** 2)
    Nx1 = int(epaisseur1 / dx1)
    
    alpha2 = lbda2 / (rho2 * cp2)
    a2 = alpha2 * dt / (dx2 ** 2)
    Nx2 = int(epaisseur2 / dx2)
    
    Nx12 = Nx1 + Nx2
    
    A = np.zeros((Nx12,Nx12))

    
    #Boucle première couche
    for i in range(1,Nx1-1):
        A[i,i] = 1 - a1
        A[i,i+1] = a1/2
        A[i,i-1] = a1/2
    
    #Boucle deuxieme couche
    for i in range(Nx1,Nx12-1):
        A[i,i] = 1 - a2
        A[i,i+1] = a2/2
        A[i,i-1] = a2/2
    
    #Point en contact couche 1 et 2
    A[Nx1-1,Nx1-2] = 0
    A[Nx1-1,Nx1] = 0
    A[Nx1-1,Nx1-1] = 0
     


    return A







