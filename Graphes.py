"============Fonction tracés de grphes=================="
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

from Modeles import m4murs_1soltp_1plafondplat, d4murs_1soltp_1plafondplat
from NRJ import perte_th,nrj_mur

"Incrémentation axe epaisseur"
def x_axe_mono(epaisseur, dx):
    Nx = int(epaisseur / dx)
    x = np.zeros(Nx)
    for i in range(Nx)     :
        x[i] = i*dx
    return x


def x_axe_double(epaisseur1, epaisseur2, dx1, dx2):
    #Nombre de points
    Nx1 = int(epaisseur1 / dx1)
    Nx2 = int(epaisseur2 / dx2)    
    Nx12 = (Nx1 + Nx2)
    x = np.zeros(Nx12)
    for i in range(Nx12):
        if i<Nx1:
            x[i] = i*dx1
        else:
            x[i] = x[i-(i-Nx1+1)] + (i-Nx1+1)*dx2
    return x

#Fonction coordonnées fonction du jour de l'année
def coordonnees_j(j,lat):
    n = 24
    h = np.zeros(n) #Hauteur soleil heure par heure
    a = np.zeros(n) #Asimut soleil heure par heure
    
    for y in range (1):
        delta = 23.45 * np.sin((2*np.pi/365 * (j+284))) 

        for t_sv in range (24):
            
            w = 15*(12-t_sv) #Rotation soleil heure par heure

            h[t_sv] = np.arcsin(np.sin(np.radians(lat)) * np.sin(np.radians(delta)) + np.cos(np.radians(lat)) *  np.cos(np.radians(delta) ) * np.cos(np.radians(w))) #déclinaison solaire
            

            if np.cos(h[t_sv]) > np.tan(np.radians(delta))/np.tan(np.radians(lat)):
            
                a[t_sv] = np.arcsin(np.cos(np.radians(delta)) * np.sin(np.radians(w)) / np.cos(np.radians(h[t_sv])))
            
            elif np.cos(h[t_sv]) < np.tan(delta)/np.tan(lat):
                a[t_sv] = np.pi - np.arcsin(np.cos(np.radians(delta)) * np.sin(np.radians(w)) / np.cos(np.radians(h[t_sv])))
                
    return h,a


"============Tracé flux solaire=================="
lat = 43
nbr_years = 1
nbr_days = nbr_years * 365
Cs = 1350

one_day = 1
nbr_days = one_day

alb_sol = 0.2 #Albédo sol (herbe)
P_diff = 0.2 #Diffusion atmosphere par temps clair (0.7 temps nuage)

h,a  = coordonnees_j(180,lat)

print(h)
print(a)

beta_mur =  np.pi /2
beta_horizontal = 0

gamma_horizontal = 0
gamma_sud = 0
gamma_ouest = np.pi/2
gamma_nord = np.pi
gamma_est = np.pi * 3/2

flux_dir_plafond = flux_dir(Cs,lat,h,a,beta_horizontal,gamma_horizontal)
flux_dif_plafond = flux_diffus(Cs, lat, h, a, beta_horizontal, gamma_horizontal, alb_sol, P_diff)
flux_ray_plafond = np.add(flux_dir_plafond,flux_dif_plafond)

flux_dir_sud = flux_dir(Cs,lat,h,a,beta_mur,gamma_sud)
flux_dif_sud = flux_diffus(Cs, lat, h, a, beta_mur, gamma_sud, alb_sol, P_diff)
flux_ray_sud = np.add(flux_dir_sud,flux_dif_sud)

flux_dir_ouest = flux_dir(Cs,lat,h,a,beta_mur,gamma_ouest)
flux_dif_ouest = flux_diffus(Cs, lat, h, a, beta_mur, gamma_ouest, alb_sol, P_diff)
flux_ray_ouest = np.add(flux_dir_ouest,flux_dif_ouest)

flux_dir_nord = flux_dir(Cs,lat,h,a,beta_mur,gamma_nord)
flux_dif_nord = flux_diffus(Cs, lat, h, a, beta_mur, gamma_nord, alb_sol, P_diff)
flux_ray_nord = np.add(flux_dir_nord,flux_dif_nord)

flux_dir_est = flux_dir(Cs,lat,h,a,beta_mur,gamma_est)
flux_dif_est = flux_diffus(Cs, lat, h, a, beta_mur, gamma_est, alb_sol, P_diff)
flux_ray_est = np.add(flux_dir_est,flux_dif_est)

x_t1 = np.linspace(0,one_day*23,one_day*24)

plt.figure(1)
plt.title("Flux solaire direct")
plt.plot(x_t1,flux_dir_plafond,label="plafond")
plt.plot(x_t1,flux_dir_sud,label="sud")
plt.plot(x_t1,flux_dir_est,label="est")
plt.plot(x_t1,flux_dir_nord,label="nord")
plt.plot(x_t1,flux_dir_ouest,label="ouest")
plt.legend()
plt.show()

plt.figure(2)
plt.title("Flux solaire diffus")
plt.plot(x_t1,flux_dif_plafond,label="plafond")
plt.plot(x_t1,flux_dif_sud,label="sud")
plt.plot(x_t1,flux_dif_est,label="est")
plt.plot(x_t1,flux_dif_nord,label="nord")
plt.plot(x_t1,flux_dif_ouest,label="ouest")
plt.legend()
plt.show()



"=============================================================="

print((1+np.cos(0))/2)
print((1+np.cos(np.pi/2))/2)

"Tracé analytique"
# print(y_anal)
# plt.figure(2)
# plt.title("Profile de température mur mono couche")
# plt.plot(y_anal, x_axe_mono(epaisseur1, dx1),label="Température (°C)")
# plt.legend()


