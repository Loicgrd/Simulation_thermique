"===============Fonction principale==============="

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
from Graphes import x_axe_mono
from Modeles import m4murs_1soltp_1plafondplat, d4murs_1soltp_1plafondplat
from NRJ import perte_th,nrj_mur,  perte_ventil, perte_menuis

#Valeur des températures exterieure sur 1 an
# T_ext = [1,1.8,1.9,0.2,-1.2,-1.8,-1.7,-2,-0.4,1.2,2,2.6,4.9,6,6.2,7.2,7.3,7,5.1,2.1,-1,-1.2,-0.3,0]
# T_ext =[21,20.8,21.2,21,20,19.2,19.7,18.9,20,20.8,23,25.2,27.4,29.3,30.5,32,32.8,33.4,33.4,33.2,33,31.2,
#            27.3,25]


"Variables temporelles"
nbr_jours = 10 #Nombre jours de simulation 
Nt = int(nbr_jours*24*4) #Nombre de point temporels du maillage
Nt_h = int(nbr_jours*24)
dt = 60*15 #Pas de temps
temps = nbr_jours*24*60*60

T_ext=[]
for j in range (nbr_jours):
    # T_ext.extend([21,20.8,21.2,21,20,19.2,19.7,18.9,20,20.8,23,25.2,27.4,29.3,30.5,32,32.8,33.4,33.4,33.2,33,31.2,27.3,25])
    T_ext.extend([1,1.8,1.9,0.2,-1.2,-1.8,-1.7,-2,-0.4,1.2,2,2.6,4.9,6,6.2,7.2,7.3,7,5.1,2.1,-1,-1.2,-0.3,0])




"Dimension batiment"
hauteur = 2.5 #Hauteur sous plafond (m)
shab = 100 #Surface habitable
V_int = hauteur * shab #Volume interieur
S_est = 25 #Surface mur est
S_nord = 25 #Surface mur nord
S_ouest = 25 #Surface mur ouest
S_sud = 25 #Surface mur sud
S_plafond = 100 #Surface du plafond
S_sol = 100 #Surface du plancher bas

#Toit pente
S_toit_nord = 0 #Surface toiture nord (m²)
angle_nord = 0 #Angle de la toiture nord (°)

S_toit_est = 0 #Surface toiture est (m²)
angle_est = 0 #Angle de la toiture est (°)

S_toit_sud = 0 #Surface toiture sud (m²)
angle_sud = 0 #Angle de la toiture sud (°)

S_toit_ouest = 0 #Surface toiture ouest (m²)
angle_ouest = 0 #Angle de la toiture ouest (°)

"Caractéristique flux solaire"
lat = 43
Cs = 1350
alb_sol = 0.2 #Albédo sol (herbe)
P_diff = 0.2 #Diffusion atmosphere par temps clair (0.7 temps nuage)

"Angle murs"
beta_mur =  np.pi /2
beta_horizontal = 0

gamma_horizontal = 0
gamma_sud = 0
gamma_ouest = np.pi/2
gamma_nord = np.pi
gamma_est = np.pi * 3/2


"Caractéristique murs"
#Béton plein
lbda1 = 1.8 #(W.m-1.K-1)
rho1 = 2300 #(kg.m-3)
cp1 = 1000  #(J.kg-1.K-1)
epaisseur1 = 0.5 #(m)
dx1 = 0.01 #Pas d'espace

#Laine de verre
lbda2 = 0.03
rho2 = 20
cp2 = 1030
epaisseur2 = 0.2
dx2 = 0.01 #Pas d'espace


"Caractéristiques fenetres"
S_fenetre = 0.750*0.600#Surface fenetre
S_porte = 2*0.90 #Surface porte
nbr_porte = 1 #Nombre de porte
nbr_fenetre = int(shab/12) #Nombre de fenetre
print(nbr_fenetre)

"Caractéristique thermique fenetre"
Uw_fenetre = 1.3 #Fenetre PVC double vitrage isolation renforcé
Uw_fenetre = 2 #Fenetre PVC double vitrage
Uw_fenetre = 4 #Fenetre bois simple vitrage

Uw_porte = 1.3 #Porte PVC


#Températures
T_sol = 13.5 #(°C) Température du sol à 1 mètre de profondeur
# T_ext = 30
T_int = 20 #Température intérieur


#Coefficients convectifs
h_ext = 15 #W/(m².K)
h_int = 5 #W/(m².K)


lat = 43 #Latitude du lieux
Cs = 1350 #Constante solaire
alb_sol = 0.2 #Albédo sol (herbe)
P_diff = 0.2 #Diffusion atmosphere par temps clair (0.7 temps nuage)

#Volume de ventilation horaire (m^3/h)
V_ventil_horaire = V_int
Q_varep = 1.3 #Coefficicient du systeme de ventilation

#Température de consigne (°C)
T_need = 21

# M_sud = solve_double(lbda1, rho1, cp1, epaisseur1, lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_ext, T_int, Nt, h_ext, h_int)
# M_est = conduc_double(lbda1, rho1, cp1, epaisseur1, lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_ext, T_int, Nt, h_ext, h_int)
# M_nord = conduc_double(lbda1, rho1, cp1, epaisseur1, lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_ext, T_int, Nt, h_ext, h_int)
# M_ouest = conduc_double(lbda1, rho1, cp1, epaisseur1, lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_ext, T_int, Nt, h_ext, h_int)

# Plafond = conduc_mono(lbda1, rho1, cp1, epaisseur1, dx1, dt, Nt, T_ext, T_int, h_ext, h_int)

# Sol = sol_tp_mono(lbda1, rho1, cp1, epaisseur1, dx1, dt, Nt, T_int, h_ext, h_int, T_sol)


delta_sud, perte_nrj, P, T_int_without, T_sud, T_est, T_nord, T_ouest, T_plafond, T_pb = m4murs_1soltp_1plafondplat(lbda1, rho1, cp1, epaisseur1, dx1, dt, Nt, Nt_h, T_ext, T_int, T_need, T_sol, h_ext, h_int, V_int, V_ventil_horaire, S_sud,S_est,S_nord,S_ouest, S_plafond, S_sol,shab,Q_varep, S_fenetre, S_porte, nbr_fenetre, nbr_porte, Uw_porte,Uw_fenetre, lat, Cs, alb_sol, P_diff, nbr_jours)

# print(T_int_without)

# print(delta_sud)
# print(P)
x = x_axe_mono(epaisseur1, dx1)

# plt.figure()
# plt.plot(x,T_sud)

# plt.figure()
# plt.plot(x,T_pb)


# plt.figure()
# plt.plot(np.linspace(0,nbr_jours,nbr_jours*24),delta_sud)




delta_sud2, perte_nrj2, P2, T_int_without2, T_sud2, T_est2, T_nord2, T_ouest2, T_plafond2, T_pb2 = d4murs_1soltp_1plafondplat(lbda1, rho1, cp1, epaisseur1, dx1, lbda2, rho2, cp2, epaisseur2, dx2, dt, Nt, Nt_h, T_ext, T_int, T_need, T_sol, h_ext, h_int, V_int, V_ventil_horaire, S_sud,S_est,S_nord,S_ouest, S_plafond, S_sol, shab,Q_varep, S_fenetre, S_porte, nbr_fenetre, nbr_porte,Uw_porte, Uw_fenetre, lat, Cs, alb_sol, P_diff, nbr_jours)

plt.figure()
plt.title("Perte énergétique au cours d'une journée")
plt.ylabel("Energie intérieure (kW)")
plt.xlabel("Temps (h)")
plt.plot(np.linspace(0,24,nbr_jours*24),perte_nrj,label="Sans isolation")
plt.plot(np.linspace(0,24,nbr_jours*24),perte_nrj2,label="Avec isolation")
plt.legend()
# plt.figure()
# plt.title("Perte énergétique au cours d'une journée")
# plt.plot(np.linspace(0,nbr_jours,nbr_jours*24),perte_nrj)




u1 = epaisseur1 + epaisseur2
u2 = int(epaisseur1/dx1 + epaisseur2/dx2)
x_db = np.linspace(0,u1,u2)


tps_simulation = 10

# "Tracé boucle température mur mono couche"
# plt.figure(10)
# plt.title("Evolution température mur béton")
# plt.xlabel("epaisseur (m)")
# plt.ylabel("Température (°C)")
# for jour_simul in range (tps_simulation):
#     Nt = int(jour_simul*24*4) #Nombre de point temporels du maillage
#     Nt_h = int(jour_simul*24)
#     dt = 60*15 #Pas de temps
#     temps = jour_simul*24*60*60
#     delta_sud, perte_nrj, P, T_int_without, T_sud, T_pb = m4murs_1soltp_1plafondplat(lbda1, rho1, cp1, epaisseur1, dx1, dt, Nt, Nt_h, T_ext, T_int, T_need, T_sol, h_ext, h_int, V_int, V_ventil_horaire, S_sud,S_est,S_nord,S_ouest, S_plafond, S_sol,shab,Q_varep, S_fenetre, S_porte, nbr_fenetre, nbr_porte, Uw_porte,Uw_fenetre, lat, Cs, alb_sol, P_diff, jour_simul)
#     plt.plot(x,T_sud,label=jour_simul)
# plt.legend()




"Tracé boucle température mur double couche"




"Tracé double couche"
plt.figure(4)
plt.title("Profil mur double couche avec isolant")
plt.ylabel("Temperature (Degres C)")
plt.xlabel("Epaisseur (m)")
plt.plot(x_db,T_sud2,label="Sud")
plt.plot(x_db,T_est2, label = "Est")
plt.plot(x_db,T_nord2, label = "Nord")
plt.plot(x_db,T_ouest2, label = "Ouest")
plt.plot(x_db,T_plafond2, label = "Plafond")
plt.plot(x_db,T_pb2, label = "Plancher bas")
plt.legend()




"Tracé mono couche"
u1 = epaisseur1 + epaisseur2
u2 = int(epaisseur1/dx1 + epaisseur2/dx2)
x_db = np.linspace(0,u1,u2)
plt.figure(5)
plt.title("Profil mur monocouche sans isolant")
plt.ylabel("Temperature (Degres C)")
plt.xlabel("Epaisseur (m)")
plt.plot(x,T_sud,label="Sud")
plt.plot(x,T_est, label = "Est")
plt.plot(x,T_nord, label = "Nord")
plt.plot(x,T_ouest, label = "Ouest")
plt.plot(x,T_plafond, label = "Plafond")
plt.plot(x,T_pb, label = "Plancher bas")
plt.legend()


