"=============================Fonctions calcul température parois========================================"

import numpy as np
from numpy.linalg import *
import matplotlib.pyplot as plt
from scipy.linalg import lu_factor
from scipy.linalg import lu_solve
import math as math
import scipy.linalg as sp

from Rayonnement import coordonnees, flux_dir, flux_diffus
from Matrices import matriceA_mono, matriceB_mono, matriceA_double, matriceB_double
from Temperature_int import temp_int, renouvellement_air


"==========Monocouche=========="
"Parois entre air exterieur et intérieur"

def conduc_mono(lbda, rho, cp, epaisseur, dx, dt, Nt, T_ext, T_int, h_ext, h_int):    

    
    Nx = int(epaisseur / dx)
    T_n = np.zeros((Nx))
    
    #Initialisation vecteur température
    T_last = np.zeros((Nx))
    # T_last[0] = T_ext
    # T_last[Nx-1] = T_int
    
    #Initialisation vecteur conditions limites
    vec = np.zeros((Nx))


    #LISTE DES H_EXT H_INT // T_EXT T_INT    
    
    Au = np.zeros((Nx-1))
    Al = np.zeros((Nx-1))
    Ad = np.zeros((Nx))
    
    for i in range (Nt-1):
        #Incrémentation temporelle des h et températures
        # h_ext = h_ext[i]
        # h_int = h_int[i]
        # T_ext = T_ext[i]
        # T_int = T_int[i]
        
        
        matA = matriceA_mono(lbda, rho, cp, epaisseur, dx ,dt, h_ext, h_int)
        matB = matriceB_mono(lbda, rho, cp, epaisseur, dx, dt)
        
        
        
        Au = np.diag(matA,1)
        Al = np.diag(matA,-1)
        Ad = np.diag(matA,0)
        
        #Actualisation températures intérieurs et exterieur pour les conditions limites
        vec[0] = T_ext
        vec[Nx-1] = T_int
        
        #Calcul vecteur memmbre gauche (vecteur CL + valeurs temps n-1)
        vec_end = np.add(vec,np.dot(matB,T_last))

        #Résolution système linéaire
        # T_n = np.linalg.solve(matA,vec_end)
        
        
        du, dd, dl, T_n, info = sp.lapack.dgtsv(Al,Ad,Au, vec_end)
        
        T_last = T_n
    
    return T_n

def solve_mono(lbda, rho, cp, epaisseur, dx, dt, Nt, T_ext, T_int, h_ext, h_int, T_last, flux_ray):    

    
    Nx = int(epaisseur / dx)
    T_n = np.zeros((Nx))
    

    
    #Initialisation vecteur conditions limites
    vec = np.zeros((Nx))


    #LISTE DES H_EXT H_INT // T_EXT T_INT    
    
    Au = np.zeros((Nx-1))
    Al = np.zeros((Nx-1))
    Ad = np.zeros((Nx))
    

        #Incrémentation temporelle des h et températures
        # h_ext = h_ext[i]
        # h_int = h_int[i]
        # T_ext = T_ext[i]
        # T_int = T_int[i]
        
        
    matA = matriceA_mono(lbda, rho, cp, epaisseur, dx ,dt, h_ext, h_int)
    matB = matriceB_mono(lbda, rho, cp, epaisseur, dx, dt)
    
    
    
    Au = np.diag(matA,1)
    Al = np.diag(matA,-1)
    Ad = np.diag(matA,0)
    
    #Actualisation températures intérieurs et exterieur pour les conditions limites
    vec[0] = T_ext + flux_ray/h_ext ########
    vec[Nx-1] = T_int ########
    
    #Calcul vecteur memmbre gauche (vecteur CL + valeurs temps n-1)
    vec_end = np.add(vec,np.dot(matB,T_last)) ########
    # vec_end = np.dot(matB,T_last)
    # vec_end = vec
    #Résolution système linéaire
    # T_n = np.linalg.solve(matA,vec_end)
    
    
    du, dd, dl, T_n, info = sp.lapack.dgtsv(Al,Ad,Au, vec_end)
    
    T_last = T_n
    
    return T_n

"Parois en contact avec le sol"
def sol_tp_mono(lbda, rho, cp, epaisseur, dx, dt, Nt, T_int, h_ext, h_int, T_sol):    

    
    Nx = int(epaisseur / dx)
    T_n = np.zeros((Nx))
    
    #Initialisation vecteur température
    T_last = np.zeros((Nx))
    # T_last[0] = T_ext
    # T_last[Nx-1] = T_int
    
    #Initialisation vecteur conditions limites
    vec = np.zeros((Nx))


    #LISTE DES H_EXT H_INT // T_EXT T_INT    
    
    Au = np.zeros((Nx-1))
    Al = np.zeros((Nx-1))
    Ad = np.zeros((Nx))
    
    for i in range (Nt-1):
        #Incrémentation temporelle des h et températures
        # h_ext = h_ext[i]
        # h_int = h_int[i]
        # T_ext = T_ext[i]
        # T_int = T_int[i]
        
        
        matA = matriceA_mono(lbda, rho, cp, epaisseur, dx ,dt, h_ext, h_int)
        matB = matriceB_mono(lbda, rho, cp, epaisseur, dx, dt)
        
        #Condition initiales 
        matA[0,0] = 1 
        matA[0,1] = 0        
        vec[0] = T_sol
        
        Au = np.diag(matA,1)
        Al = np.diag(matA,-1)
        Ad = np.diag(matA,0)
        
        #Actualisation températures intérieurs pour les conditions limites
        vec[Nx-1] = T_int
        
        #Calcul vecteur memmbre gauche (vecteur CL + valeurs temps n-1)
        # vec_end = np.add(vec,np.dot(matB,T_last))
        vec_end = np.dot(matB,T_last)
        #Résolution système linéaire
        # T_n = np.linalg.solve(matA,vec_end)
        
        
        du, dd, dl, T_n, info = sp.lapack.dgtsv(Al,Ad,Au, vec_end)
        
        T_last = T_n
    
    return T_n


def solve_tp_mono(lbda, rho, cp, epaisseur, dx, dt, Nt, T_int, h_ext, h_int, T_sol, T_last):    

    
    Nx = int(epaisseur / dx)
    T_n = np.zeros((Nx))
    

    
    #Initialisation vecteur conditions limites
    vec = np.zeros((Nx))


    #LISTE DES H_EXT H_INT // T_EXT T_INT    
    
    Au = np.zeros((Nx-1))
    Al = np.zeros((Nx-1))
    Ad = np.zeros((Nx))
    

    #Incrémentation temporelle des h et températures
    # h_ext = h_ext[i]
    # h_int = h_int[i]
    # T_ext = T_ext[i]
    # T_int = T_int[i]
    
    
    matA = matriceA_mono(lbda, rho, cp, epaisseur, dx ,dt, h_ext, h_int)
    matB = matriceB_mono(lbda, rho, cp, epaisseur, dx, dt)
    
    #Condition initiales 
    matA[0,0] = 1 
    matA[0,1] = 0        
    vec[0] = T_sol
    
    Au = np.diag(matA,1)
    Al = np.diag(matA,-1)
    Ad = np.diag(matA,0)
    
    #Actualisation températures intérieurs pour les conditions limites
    vec[Nx-1] = T_int
    
    #Calcul vecteur memmbre gauche (vecteur CL + valeurs temps n-1)
    vec_end = np.add(vec,np.dot(matB,T_last))
    
    #Résolution système linéaire
    # T_n = np.linalg.solve(matA,vec_end)
    
    
    du, dd, dl, T_n, info = sp.lapack.dgtsv(Al,Ad,Au, vec_end)
    
    T_last = T_n
    
    return T_n


"===============Double couchess==============="
"Parois entre air exterieur et intérieur"
def conduc_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_ext, T_int, Nt, h_ext, h_int):    
    
    # matB = matriceB_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt,h_ext,h_int,T_ext,T_int )
    
    #Nombre de points
    Nx1 = int(epaisseur1 / dx1)
    Nx2 = int(epaisseur2 / dx2)    
    Nx12 = Nx1 + Nx2
    
    T_n = np.zeros((Nx12))
    
   
    
    #Initialisation vecteur température
    T_last = np.zeros((Nx12))
    # T_last[0] = T_ext
    # T_last[Nx12-1] = T_int
    
    #Initialisation vecteur conditions limites
    vec = np.zeros((Nx12))


    #LISTE DES H_EXT H_INT // T_EXT T_INT
    
    
    Au = np.zeros((Nx12-1))
    Al = np.zeros((Nx12-1))
    Ad = np.zeros((Nx12))

    for i in range (Nt-1):
        
       
        
        #Actualisation matrices A et B
        matB = matriceB_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt)
        matA = matriceA_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, h_ext, h_int )

        Au = np.diag(matA,1)
        Al = np.diag(matA,-1)
        Ad = np.diag(matA,0)

        

        #Actualisation températures intérieurs et exterieur pour les conditions limites
        vec[0] = T_ext
        vec[Nx12-1] = T_int
        
        #Calcul vecteur membre gauche (vecteur CL + valeurs temps n-1)
        vec_end = np.add(vec,np.dot(matB,T_last))
        
 
        
        #Résolution système linéaire
        #T_n = sp.solve(matA,vec_end)
        
        du, dd, dl, T_n, info = sp.lapack.dgtsv(Al,Ad,Au, vec_end)
        

        
        T_last = T_n

    
    return T_n


def solve_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_ext, T_int, h_ext, h_int, T_last, flux_ray):    
    
    # matB = matriceB_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt,h_ext,h_int,T_ext,T_int )
    
    #Nombre de points
    Nx1 = int(epaisseur1 / dx1)
    Nx2 = int(epaisseur2 / dx2)    
    Nx12 = Nx1 + Nx2
    
    T_n = np.zeros((Nx12))
    
   

    
    #Initialisation vecteur conditions limites
    vec = np.zeros((Nx12))


    #LISTE DES H_EXT H_INT // T_EXT T_INT
    
    
    Au = np.zeros((Nx12-1))
    Al = np.zeros((Nx12-1))
    Ad = np.zeros((Nx12))

   
    #Actualisation matrices A et B
    matB = matriceB_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt)
    matA = matriceA_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, h_ext, h_int )

    Au = np.diag(matA,1)
    Al = np.diag(matA,-1)
    Ad = np.diag(matA,0)

    

    #Actualisation températures intérieurs et exterieur pour les conditions limites
    vec[0] = T_ext + flux_ray/h_ext
    vec[Nx12-1] = T_int
    
    #Calcul vecteur membre gauche (vecteur CL + valeurs temps n-1)
    vec_end = np.add(vec,np.dot(matB,T_last))
    
 
    
    #Résolution système linéaire
    #T_n = sp.solve(matA,vec_end)
    
    du, dd, dl, T_n, info = sp.lapack.dgtsv(Al,Ad,Au, vec_end)
    

    
    T_last = T_n

    
    return T_n



"Parois en contact avec le sol"

def sol_tp_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_int, Nt, h_int, h_ext, T_sol): #Fonction Calcul sol sur terre plein
    Nx1 = int(epaisseur1 / dx1)
    Nx2 = int(epaisseur2 / dx2)    
    Nx12 = Nx1 + Nx2
    
    T_n = np.zeros((Nx12))
    
   
    
    #Initialisation vecteur température
    T_last = np.zeros((Nx12))
    
    #Initialisation vecteur conditions limites
    vec = np.zeros((Nx12))


    #LISTE DES H_EXT H_INT // T_EXT T_INT   
    
    Au = np.zeros((Nx12-1))
    Al = np.zeros((Nx12-1))
    Ad = np.zeros((Nx12)) 

    for i in range(Nt-1):
        
        matB = matriceB_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt)
        matA = matriceA_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, h_ext, h_int )
        
        
        
        #Condition initiales 
        matA[0,0] = 1 
        matA[0,1] = 0        
        vec[0] = T_sol
        
        Au = np.diag(matA,1)
        Al = np.diag(matA,-1)
        Ad = np.diag(matA,0)
        
        
        
        
        
        
        #Actualisation températures intérieurs et exterieur pour les conditions limites
        vec[Nx12-1] = T_int
        
        #Calcul vecteur membre gauche (vecteur CL + valeurs temps n-1)
        vec_end = np.add(vec,np.dot(matB,T_last))
        
 
        
        #Résolution système linéaire
        #T_n = sp.solve(matA,vec_end)
        
        du, dd, dl, T_n, info = sp.lapack.dgtsv(Al,Ad,Au, vec_end)
        

        
        T_last = T_n
        
        

    return T_n



def solve_tp_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_int, h_int, h_ext, T_sol, T_last): #Fonction Calcul sol sur terre plein
    Nx1 = int(epaisseur1 / dx1)
    Nx2 = int(epaisseur2 / dx2)    
    Nx12 = Nx1 + Nx2
    
    T_n = np.zeros((Nx12))
    
   
    
    
    #Initialisation vecteur conditions limites
    vec = np.zeros((Nx12))


    #LISTE DES H_EXT H_INT // T_EXT T_INT   
    
    Au = np.zeros((Nx12-1))
    Al = np.zeros((Nx12-1))
    Ad = np.zeros((Nx12)) 

        
    matB = matriceB_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt)
    matA = matriceA_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, h_ext, h_int )
    
    
    
    #Condition initiales 
    matA[0,0] = 1 
    matA[0,1] = 0        
    vec[0] = T_sol
    
    Au = np.diag(matA,1)
    Al = np.diag(matA,-1)
    Ad = np.diag(matA,0)
    
    
    
    
    
    
    #Actualisation températures intérieurs et exterieur pour les conditions limites
    vec[Nx12-1] = T_int
    
    #Calcul vecteur membre gauche (vecteur CL + valeurs temps n-1)
    vec_end = np.add(vec,np.dot(matB,T_last))
    
 
    
    #Résolution système linéaire
    #T_n = sp.solve(matA,vec_end)
    
    du, dd, dl, T_n, info = sp.lapack.dgtsv(Al,Ad,Au, vec_end)
    

    
    T_last = T_n
        
        

    return T_n


