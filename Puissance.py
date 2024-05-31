import numpy as np
from numpy.linalg import *
import matplotlib.pyplot as plt
from scipy.linalg import lu_factor
from scipy.linalg import lu_solve
import math as math
import scipy.linalg as sp

from Rayonnement import coordonnees, flux_dir, flux_diffus
from Parois import solve_mono, solve_double, solve_tp_mono, solve_tp_double
from Matrices import matriceA_mono, matriceB_mono, matriceA_double, matriceB_double
from Temperature_int import temp_int, renouvellement_air
from Graphes import x_axe_mono, x_axe_double
from Modeles import m4murs_1soltp_1plafondplat





def nrj_mur(h_int,h_ext,T_mur,S_mur,dx1,rho1,cp1,epaisseur1,dx2=0,rho2=0,cp2=0,epaisseur2=0):
    nrj = 0
    # for i in range(len(T_mur)-1):
    #     if i*dx1<epaisseur1:
    #         nrj = nrj + rho1 * cp1 * dx1 * S_mur * (T_mur[i] + 273)
    #     else:
    #         nrj = nrj + rho2 * cp2 * dx2 * S_mur * (T_mur[i] + 273)
    nrj = nrj/3600 #Conversion j en wh
    return nrj

def R_cond(lbda,epaisseur):
    if lbda != 0:
        R = epaisseur/lbda
    else:
        R=0
    return R


def R_conv(h):
    R = 1/h
    return R

def R_tot(h_int,h_ext,lbda1,epaisseur1,lbda2=0,epaisseur2=0):
    R = R_cond(lbda1,epaisseur1) + R_cond(lbda2,epaisseur2) + R_conv(h_int) + R_conv(h_ext)    
    return R

def perte_th(h_int,h_ext,T_int,T_ext,lbda1,epaisseur1,lbda2=0,epaisseur2=0): #Si positif flux entrant (Référientiel maison)
    dep = 1 / R_tot(h_int,h_ext,lbda1,epaisseur1,lbda2,epaisseur2) * (T_ext - T_int)
    return dep
    














