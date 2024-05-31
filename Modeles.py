"================Fonctions pour plusieurs modÃ¨les de maison============="
import numpy as np
from numpy.linalg import *
import matplotlib.pyplot as plt
from scipy.linalg import lu_factor
from scipy.linalg import lu_solve
import math as math
import scipy.linalg as sp

from NRJ import perte_th, nrj_mur,  perte_ventil, perte_menuis
from Rayonnement import coordonnees, flux_dir, flux_diffus
from Parois import solve_mono, solve_double, solve_tp_mono, solve_tp_double
from Matrices import matriceA_mono, matriceB_mono, matriceA_double, matriceB_double
from Temperature_int import temp_int, renouvellement_air, puissance_chauffe



# 4 murs non isolÃ© 1 plafond plat non isolÃ© 1 plancher non isolÃ© sur terre plein isolÃ©

def m4murs_1soltp_1plafondplat(lbda, rho, cp, epaisseur, dx, dt, Nt, Nt_h, T_ext, T_int, T_need, T_sol, h_ext, h_int, V_int, V_ventil_horaire, S_sud,S_est,S_nord,S_ouest, S_plafond, S_sol, shab,Q_varep,S_fenetre,S_porte,nbr_fenetre,nbr_porte, Uw_porte, Uw_fenetre, lat, Cs, alb_sol, P_diff, nbr_days):



    P_use = np.zeros(Nt_h)    
    T_int_without = np.zeros(Nt_h) #Température sans chauffage   
    delta_nrj_sud = np.zeros(Nt_h)
    perte_nrj_mur = np.zeros(Nt_h)
    
    Nx = int(epaisseur / dx)
    T_sud = np.zeros((Nx))
    T_est = np.zeros((Nx))
    T_nord = np.zeros((Nx))
    T_ouest = np.zeros((Nx))
    T_plafond = np.zeros((Nx))
    T_pb = np.zeros((Nx))
    
    
    
    
    "Calcul flux solaire parois"



    beta_mur =  np.pi /2
    beta_horizontal = 0

    gamma_horizontal = 0
    gamma_sud = 0
    gamma_ouest = np.pi/2
    gamma_nord = np.pi
    gamma_est = np.pi * 3/2
    
    h,a = coordonnees(nbr_days,lat)
    
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
    
    
    
    
    "=========================="
    
    
    # T_sud[0] = T_ext[0]
    # T_sud[Nx-1] = T_int
    # T_est[0] = T_ext[0]
    # T_est[Nx-1] = T_int
    # T_nord[0] = T_ext[0]
    # T_nord[Nx-1] = T_int
    # T_ouest[0] = T_ext[0]
    # T_ouest[Nx-1] = T_int
    # T_plafond[0] = T_ext[0]
    # T_plafond[Nx-1] = T_int
    # T_pb[0] = T_ext[0]
    # T_pb[Nx-1] = T_int
    
    # print("lol")
    # print(T_sud[0])
    
    
    S_murs = S_sud + S_est + S_ouest + S_nord
    nrj_sud_last = 0
    
    for i in range (Nt_h):

         
  
        
        
        for j in range (4): #Boucle 4 fois 15 min
        

        
        
        
        
            T_sud = solve_mono(lbda, rho, cp, epaisseur, dx, dt, Nt, T_ext[i], T_int, h_ext, h_int,T_sud,flux_ray_sud[i]) 
            T_est = solve_mono(lbda, rho, cp, epaisseur, dx, dt, Nt, T_ext[i], T_int, h_ext, h_int, T_est,flux_ray_est[i]) 
            T_nord = solve_mono(lbda, rho, cp, epaisseur, dx, dt, Nt, T_ext[i], T_int, h_ext, h_int, T_nord,flux_ray_nord[i])  
            T_ouest = solve_mono(lbda, rho, cp, epaisseur, dx, dt, Nt, T_ext[i], T_int, h_ext, h_int, T_ouest,flux_ray_ouest[i]) 
            T_plafond = solve_mono(lbda, rho, cp, epaisseur, dx, dt, Nt, T_ext[i], T_int, h_ext, h_int, T_plafond, flux_ray_plafond[i]) 
            T_pb = solve_tp_mono(lbda, rho, cp, epaisseur, dx, dt, Nt, T_int, h_ext, h_int, T_sol, T_pb) 
            
            T_int_without_p = temp_int(h_int,V_int,S_sud,S_est,S_nord,S_ouest, S_plafond, S_sol, T_int, T_sud, T_est, T_nord, T_ouest, T_plafond, T_pb, dt)

            # T_sud = solve_mono(lbda, rho, cp, epaisseur, dx, dt, Nt, T_ext, T_int, h_ext, h_int)
            # T_est = solve_mono(lbda, rho, cp, epaisseur, dx, dt, Nt, T_ext, T_int, h_ext, h_int)
            # T_nord = solve_mono(lbda, rho, cp, epaisseur, dx, dt, Nt, T_ext, T_int, h_ext, h_int)
            # T_ouest = solve_mono(lbda, rho, cp, epaisseur, dx, dt, Nt, T_ext, T_int, h_ext, h_int)
            # T_plafond = solve_mono(lbda, rho, cp, epaisseur, dx, dt, Nt, T_ext, T_int, h_ext, h_int) 
            # T_sol = solve_tp_mono(lbda, rho, cp, epaisseur, dx, dt, Nt, T_ext, T_int, h_ext, h_int)
        
        
            # T_int = renouvellement_air(T_ext[i],T_int,V_int,V_ventil_horaire) #Température après ventilation
        
            # print(T_nord[Nx-1])
        
        
        T_int_without[i] = T_int_without_p
    
    
        # P_use[i] = puissance_chauffe(T_ext[i],T_int,T_need,V_int) #Puissance requise pour augmenter la température à T consigne

        # P_use[i] = h_int * (S_sud * T_sud[Nx-1] + S_est*T_est[Nx-1] + S_nord*T_nord[Nx-1] + S_ouest*T_ouest[Nx-1] + S_plafond*T_plafond[Nx-1] + S_sol*T_sol[Nx-1]) / 1000
    
        # T_int = T_need #Après chauffage à température de consigne



        #NRJ MURS
        nrj_sud = nrj_mur(h_int,h_ext,T_sud,S_sud,dx,rho,cp,epaisseur)
        
        if i!= 0:
            
            delta_nrj_sud[i] = nrj_sud - nrj_sud_last #Variation énergie du mur sur 1h
        
        nrj_sud_last = nrj_sud
        #===========================
        #Calcul surface de fenetre
        
        #Calcul surface porte
        
        
        #Energie à compenser (Murs, plafond, plancher, ventilation, porte, fenetres)
        perte_nrj_mur[i] = perte_th(S_murs,h_int,h_ext,T_int,T_ext[i],lbda,epaisseur) + perte_th(S_plafond,h_int,h_ext,T_int,T_ext[i],lbda,epaisseur) + perte_th(S_murs,h_int,0,T_int,T_ext[i],lbda,epaisseur) + perte_ventil(T_int,T_ext[i],shab,Q_varep) + perte_menuis(Uw_porte,nbr_porte,S_porte, T_ext[i], T_int) + perte_menuis(Uw_fenetre,nbr_fenetre,S_fenetre, T_ext[i], T_int)
        
        
    return delta_nrj_sud,perte_nrj_mur,P_use, T_int_without, T_sud, T_est, T_nord, T_ouest, T_plafond, T_pb


def d4murs_1soltp_1plafondplat(lbda1, rho1, cp1, epaisseur1, dx1, lbda2, rho2, cp2, epaisseur2, dx2, dt, Nt, Nt_h, T_ext, T_int, T_need, T_sol, h_ext, h_int, V_int, V_ventil_horaire, S_sud,S_est,S_nord,S_ouest, S_plafond, S_sol, shab,Q_varep, S_fenetre,S_porte,nbr_fenetre,nbr_porte, Uw_porte, Uw_fenetre, lat, Cs, alb_sol, P_diff, nbr_days):

    P_use = np.zeros(Nt_h)    
    T_int_without = np.zeros(Nt_h*4) #Température sans chauffage   
    delta_nrj_sud = np.zeros(Nt_h)
    perte_nrj_mur = np.zeros(Nt_h)
    
    Nx = int(epaisseur1 / dx1) + int(epaisseur2 / dx2)
    T_sud = np.zeros((Nx))
    T_est = np.zeros((Nx))
    T_nord = np.zeros((Nx))
    T_ouest = np.zeros((Nx))
    T_plafond = np.zeros((Nx))
    T_pb = np.zeros((Nx))
    
    
    "Calcul flux solaire parois"



    beta_mur =  np.pi /2
    beta_horizontal = 0

    gamma_horizontal = 0
    gamma_sud = 0
    gamma_ouest = np.pi/2
    gamma_nord = np.pi
    gamma_est = np.pi * 3/2
    
    h,a = coordonnees(nbr_days,lat)
    
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
    
    
    
    
    
    
    S_murs = S_sud + S_est + S_ouest + S_ouest
    nrj_sud_last = 0
    
    for i in range (Nt_h):
        
        
        
        for j in range (4): #Boucle 4 fois 15 min
            T_sud = solve_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_ext[i], T_int, h_ext, h_int, T_sud, flux_ray_sud[i])
            T_est = solve_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_ext[i], T_int, h_ext, h_int, T_est, flux_ray_est[i])
            T_nord = solve_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_ext[i], T_int, h_ext, h_int, T_nord, flux_ray_nord[i])
            T_ouest = solve_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_ext[i], T_int, h_ext, h_int, T_ouest, flux_ray_ouest[i])
            T_plafond = solve_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_ext[i], T_int, h_ext, h_int, T_plafond, flux_ray_plafond[i])
            T_pb = solve_tp_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_int, h_int, h_ext, T_sol, T_pb)
            
            T_int_without_p = temp_int(h_int,V_int,S_sud,S_est,S_nord,S_ouest, S_plafond, S_sol, T_int, T_sud, T_est, T_nord, T_ouest, T_plafond, T_pb, dt)

            # T_sud = solve_mono(lbda, rho, cp, epaisseur, dx, dt, Nt, T_ext, T_int, h_ext, h_int)
            # T_est = solve_mono(lbda, rho, cp, epaisseur, dx, dt, Nt, T_ext, T_int, h_ext, h_int)
            # T_nord = solve_mono(lbda, rho, cp, epaisseur, dx, dt, Nt, T_ext, T_int, h_ext, h_int)
            # T_ouest = solve_mono(lbda, rho, cp, epaisseur, dx, dt, Nt, T_ext, T_int, h_ext, h_int)
            # T_plafond = solve_mono(lbda, rho, cp, epaisseur, dx, dt, Nt, T_ext, T_int, h_ext, h_int) 
            # T_sol = solve_tp_mono(lbda, rho, cp, epaisseur, dx, dt, Nt, T_ext, T_int, h_ext, h_int)
        
        
            # T_int = renouvellement_air(T_ext[i],T_int,V_int,V_ventil_horaire) #Température après ventilation
        
            # print(T_nord[Nx-1])
        
        
        T_int_without[i] = T_int_without_p
    
    
        # P_use[i] = puissance_chauffe(T_ext[i],T_int,T_need,V_int) #Puissance requise pour augmenter la température à T consigne

        # P_use[i] = h_int * (S_sud * T_sud[Nx-1] + S_est*T_est[Nx-1] + S_nord*T_nord[Nx-1] + S_ouest*T_ouest[Nx-1] + S_plafond*T_plafond[Nx-1] + S_sol*T_sol[Nx-1]) / 1000
    
        # T_int = T_need #Après chauffage à température de consigne
        
        
        
        #NRJ MURS
        nrj_sud = nrj_mur(h_int,h_ext,T_sud,S_sud,dx1,rho1,cp1,epaisseur1,dx2,rho2,cp2,epaisseur2)
        
        
        if i != 0:
        
            delta_nrj_sud[i] = nrj_sud - nrj_sud_last #Variation énergie du mur sur 1h
        
        nrj_sud_last = nrj_sud
        #===========================
        
        #Energie à compenser (Murs, plafond, plancher, ventilation)
        perte_nrj_mur[i] = perte_th(S_murs,h_int,h_ext,T_int,T_ext[i],lbda1,epaisseur1,lbda2,epaisseur2) + perte_th(S_plafond,h_int,h_ext,T_int,T_ext[i],lbda1,epaisseur1,lbda2,epaisseur2) + perte_th(S_sol,h_int,0,T_int,T_sol,lbda1,epaisseur1,lbda2,epaisseur2) + perte_ventil(T_int,T_ext[i],shab,Q_varep) + perte_menuis(Uw_porte,nbr_porte,S_porte, T_ext[i], T_int) + perte_menuis(Uw_fenetre,nbr_fenetre,S_fenetre, T_ext[i], T_int)
        
        
        
    return delta_nrj_sud,perte_nrj_mur,P_use, T_int_without, T_sud, T_est, T_nord, T_ouest, T_plafond, T_pb
