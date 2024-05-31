'=======================Fonctions rayonnement=================================='


import numpy as np
from numpy.linalg import *
import matplotlib.pyplot as plt
from scipy.linalg import lu_factor
from scipy.linalg import lu_solve
import math as math
import scipy.linalg as sp




def coordonnees(nbr_days,lat):
    n = nbr_days*24
    h = np.zeros(n) #Hauteur soleil heure par heure
    a = np.zeros(n) #Asimut soleil heure par heure
    
    for y in range (nbr_days):
        delta = 23.45 * np.sin((2*np.pi/365 * (y+284))) 

        for t_sv in range (24):
            i = y*24 + t_sv
            
            w = 15*(12-t_sv) #Rotation soleil heure par heure

            h[i] = np.arcsin(np.sin(np.radians(lat)) * np.sin(np.radians(delta)) + np.cos(np.radians(lat)) *  np.cos(np.radians(delta) ) * np.cos(np.radians(w))) #déclinaison solaire
            

            if np.cos(h[i]) > np.tan(np.radians(delta))/np.tan(np.radians(lat)):
            
                a[i] = np.arcsin(np.cos(np.radians(delta)) * np.sin(np.radians(w)) / np.cos(np.radians(h[i])))
            
            elif np.cos(h[i]) < np.tan(delta)/np.tan(lat):
                a[i] = np.pi - np.arcsin(np.cos(np.radians(delta)) * np.sin(np.radians(w)) / np.cos(np.radians(h[i])))
                
    return h,a






def flux_dir(Cs,lat,h,a,beta,gamma): #Calcul des flux direct
    n = len(h)

    flux = np.zeros(n)

    for i in range (n):

        if h[i] > 0:

            flux[i] = max(Cs * (np.sin(beta) * np.cos(h[i]) * np.cos ( a[i] + gamma) + np.cos(beta) * np.sin(h[i])) * 0.8,0)
    
            
    return flux
    

def flux_atm(Cs,lat,h,a,beta,gamma): #Calcul des flux direct
    n = len(h)

    flux = np.zeros(n)

    for i in range (n):

        if h[i] > 0:

            flux[i] = max(Cs * (np.sin(beta) * np.cos(h[i]) * np.cos ( a[i] + gamma) + np.cos(beta) * np.sin(h[i])),0)
    
            
    return flux
    


def flux_diffus(Cs,lat,h,a,beta,gamma,alb_sol,P_diff):
    n = len(h)
    flux_diffus = np.zeros(n)
    flux_direct_sol = flux_dir(Cs,lat,h,a,0,0) #Flux arrivant au sol
    flux_atmm = flux_atm(Cs,lat,h,a,0,0)
    
    F_ciel = (1+np.cos(beta))/2 #Facteur de forme parois ciel
    F_sol = (1-np.cos(beta))/2 #Facteur de forme parois sol

    
    
    for i in range (n):

        flux_diffus[i] = F_ciel * P_diff * flux_atmm[i]  + F_sol * alb_sol * flux_direct_sol[i]

            

    return flux_diffus

