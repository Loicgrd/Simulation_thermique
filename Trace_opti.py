# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 18:02:35 2022

@author: loicg
"""
import numpy as np
from numpy.linalg import *
import matplotlib.pyplot as plt
from scipy.linalg import lu_factor
from scipy.linalg import lu_solve
import math as math
import sys
import scipy.linalg as sp
import time


nbr_jours = 10 #Nombre jours de simulation 
Nt = int(nbr_jours*24) #Nombre de point temporels du maillage
dt = 60*60 #Pas de temps
temps = nbr_jours*24*60*60


#Températures intérieur et exterieur
T_ext = 30
T_int = 20


#Coefficients convectifs
h_ext = 15 #W/(m².K)
h_int = 5 #W/(m².K)


class couche:
    name = str
    lbda = float
    rho = float
    cp = float
    epaisseur = float

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

# def pas_temps(a, h, lbda, dx):
#     #c=0.35
#     #temp =int( c*(lbda/(h*math.sqrt(a)))**2)
#     temp = int(dx**2/alpha)
#     return temp

# def pas_espace(h, lbda):
#     c=0.35
#     temp = float(c * lbda / h)
#     return temp

# dx1 = pas_espace(h_ext, lbda1)


# alpha = lbda1 / (rho1 * cp1)
# dt = pas_temps(alpha, h_ext, lbda1, dx1)


"===================Fonctions monocouches========================="

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

def matriceB_mono(lbda, rho, cp, epaisseur, dx, dt):
    alpha = lbda / (rho * cp)
    a = alpha * dt / (dx ** 2)
    Nx = int(epaisseur / dx)
    A = np.zeros((Nx,Nx))
    
    
    
    
    for i in range(1,Nx-1):
        A[i,i] = 1 - a
        A[i,i+1] = a/2
        A[i,i-1] = a/2
    


    
    return A


    
    
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
    

    
    for i in range (Nt-1):
        #Incrémentation temporelle des h et températures
        # h_ext = h_ext[i]
        # h_int = h_int[i]
        # T_ext = T_ext[i]
        # T_int = T_int[i]
        
        
        matA = matriceA_mono(lbda, rho, cp, epaisseur, dx ,dt, h_ext, h_int)
        matB = matriceB_mono(lbda, rho, cp, epaisseur, dx, dt)
        
        #Actualisation températures intérieurs et exterieur pour les conditions limites
        vec[0] = T_ext
        vec[Nx-1] = T_int
        
        #Calcul vecteur memmbre gauche (vecteur CL + valeurs temps n-1)
        vec_end = np.add(vec,np.dot(matB,T_last))
        
        #Résolution système linéaire
        T_n = np.linalg.solve(matA,vec_end)
        
        
        T_last = T_n
    
    return T_n

def x_axe_mono(epaisseur, dx):
    Nx = int(epaisseur / dx)
    x = np.zeros(Nx)
    for i in range(Nx)     :
        x[i] = i*dx
    return x




"========================Fonctions double couches======================="



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
    for i in range(1,Nx1-1):
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




    
def conduc_double_sp(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_ext, T_int, Nt, h_ext, h_int):    
    
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


    for i in range (Nt-1):
        
       
        
        #Actualisation matrices A et B
        matB = matriceB_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt)
        matA = matriceA_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, h_ext, h_int )

        #Actualisation températures intérieurs et exterieur pour les conditions limites
        vec[0] = T_ext
        vec[Nx12-1] = T_int
        
        #Calcul vecteur membre gauche (vecteur CL + valeurs temps n-1)
        vec_end = np.add(vec,np.dot(matB,T_last))
        
        #Résolution système linéaire
        T_n = sp.solve(matA,vec_end)
        
        
        T_last = T_n

    
    return T_n




def conduc_double_np(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_ext, T_int, Nt, h_ext, h_int):    
    
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


    for i in range (Nt-1):
        
       
        
        #Actualisation matrices A et B
        matB = matriceB_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt)
        matA = matriceA_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, h_ext, h_int )

        #Actualisation températures intérieurs et exterieur pour les conditions limites
        vec[0] = T_ext
        vec[Nx12-1] = T_int
        
        #Calcul vecteur membre gauche (vecteur CL + valeurs temps n-1)
        vec_end = np.add(vec,np.dot(matB,T_last))
        
        #Résolution système linéaire
        T_n = np.linalg.solve(matA,vec_end)
        
        
        T_last = T_n

    
    return T_n



def conduc_double_opti(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_ext, T_int, Nt, h_ext, h_int):    
    
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


"==============Solution analytique==================="
def anal_sol(lbda,rho,cp,epaisseur,dx,temps):
    
    # Dt = np.zeros((Nt))
    Nx = int(epaisseur / dx)
    # for i in range(Nt):
    #     Dt[i] = i*dt
    temp = np.zeros((Nx))
    temp[0] = 20
    temp[Nx-1] = 30
    alpha = lbda/(rho*cp)
    
    Dx = x_axe_mono(epaisseur, dx)
    
    T0 = 0
    teta0 = (T0 - temp[0])/(temp[0] - temp[Nx-1]) 
    for i in range (1,len(Dx)-1):
        
        summ=0

        
        for j in range (1,100):
            
            # summ = 1/n*(np.exp(-(n*np.pi/(2*epaisseur)*alpha*temps)))*np.sin(n*np.pi*i*dx/(2*epaisseur))
            
            # n=n+2
            Cn = 2 * ((-1)**j*(1-teta0) + teta0)/( j * np.pi)
            
            summ = summ + Cn * np.sin(j * np.pi * i*dx / epaisseur) * np.exp(-(j*np.pi/epaisseur)**2 * alpha * temps)
            
        temp[i] = temp[0] + i*dx*(temp[Nx-1] - temp[0])/epaisseur + (temp[Nx-1] - temp[0])*summ
    
    
    
    return temp

# y_anal = anal_sol(lbda1,rho1,cp1,epaisseur1,dx1,temps)

"=======================Tracé graphiques======================"

# x_plot = x_axe_double(epaisseur1, epaisseur2, dx1, dx2)


# YY = anal_sol(lbda1,rho1,cp1,epaisseur1,dx1,temps)



"Tracé boucle double couche"
# plt.figure(3)
# plt.title("Evolution toutes les 8h du profile de température mur deux couches")
# for i in range (0,168,8) :
    
#     Nt = int(i) #Nombre de point temporels du maillage
    
#     y_plot = conduc_double(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_ext, T_int, Nt, h_ext, h_int)
#     plt.plot(x_plot,y_plot)
# plt.ylabel("Température (°C)")
# plt.xlabel("Epaisseur (m)")
# plt.legend()
    

"Tracé boucle mono"
# plt.figure(4)
# plt.title("Evolution toutes les 8h du profile de température mur monocouche")
# for i in range (0,168,8) :
    
#     Nt = int(i) #Nombre de point temporels du maillage
    
#     y_plot = conduc_mono(lbda1, rho1, cp1, epaisseur1, dx1, dt, Nt, T_ext, T_int, h_ext, h_int)
#     plt.plot(x_axe_mono(epaisseur1, dx1),y_plot)
# plt.ylabel("Température (°C)")
# plt.xlabel("Epaisseur (m)")
# plt.legend()

"Tracé mono"
# tc=conduc_mono(lbda1, rho1, cp1, epaisseur1, dx1, dt, Nt, T_ext, T_int, h_ext, h_int)
# plt.figure(2)
# plt.title("Profile de température mur mono couche")
# plt.plot( x_axe_mono(epaisseur1, dx1),tc,label="Température (°C)")
# plt.ylabel("Température (°C)")
# plt.xlabel("Epaisseur (m)")
# plt.legend()

"Tracé analytique"
# print(y_anal)
# plt.figure(2)
# plt.title("Profile de température mur mono couche")
# plt.plot(y_anal, x_axe_mono(epaisseur1, dx1),label="Température (°C)")
# plt.legend()

"Tracé compraison numérique analytique"
# plt.figure(21)
# plt.title()
# opti_plot = conduc_double_opti(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_ext, T_int, Nt, h_ext, h_int)
# y_plot = conduc_double_np(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_ext, T_int, Nt, h_ext, h_int)
"=========================Comparaison incrémentation temporelle================================"
# j_max = 365

# TTopt = np.zeros((j_max))

# TTnoptnp = np.zeros((j_max))

# TTnoptsp = np.zeros((j_max))

# for nbr_jours in range (j_max):
    
    
#     Nt = int(nbr_jours*24) #Nombre de point temporels du maillage
#     dt = 60*60 #Pas de temps
#     # temps = nbr_jours*24*60*60


#     tps1=time.clock()
#     opti_plot = conduc_double_opti(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_ext, T_int, Nt, h_ext, h_int)
#     tps2=time.clock()
#     TTopt[nbr_jours] = tps2-tps1
#     print("Temps opti : ", tps2-tps1)
    
    
    
#     tp1 = time.clock()
#     y_plot = conduc_double_sp(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_ext, T_int, Nt, h_ext, h_int)
#     tp2= time.clock()
#     TTnoptsp[nbr_jours] = tp2-tp1
#     print("Temps non opti : ",tp2-tp1)
    
    
#     t1 = time.clock()
#     y_plot = conduc_double_np(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_ext, T_int, Nt, h_ext, h_int)
#     t2= time.clock()
#     TTnoptnp[nbr_jours] = t2-t1
#     print("Temps non opti : ",t2-t1)

# XX = np.linspace(0,j_max-1,num = j_max)

# plt.figure(6)
# plt.title("Temps d'éxecution")
# plt.plot(XX,TTopt,label="Temps dgtsv")
# plt.plot(XX,TTnoptsp,label="Temps solveur numpy")
# plt.plot(XX,TTnoptnp,label="Temps solveur scipy")
# plt.ylabel("Temps d'exécution (s)")
# plt.xlabel("Nombre de jours de simulation")
# plt.legend()
"============================================================================================="



"=================Comparaison incrémentation spatial======================"

nbr_jours = 1 #Nombre jours de simulation 
Nt = int(nbr_jours*24) #Nombre de point temporels du maillage
dt = 60*60 #Pas de temps
temps = nbr_jours*24*60*60



DDopt = np.zeros((100))

DDnoptnp = np.zeros((100))

DDnoptsp = np.zeros((100))

XXD = np.zeros((100))

for nbr_DD in range (100):
    
    XXD[nbr_DD] = (epaisseur1 + epaisseur2) / dx1


    tps1=time.clock()
    opti_plot = conduc_double_opti(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_ext, T_int, Nt, h_ext, h_int)
    tps2=time.clock()
    DDopt[nbr_DD] = tps2-tps1
    # print("Temps opti : ", tps2-tps1)
    
    
    
    tp1 = time.clock()
    y_plot = conduc_double_np(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_ext, T_int, Nt, h_ext, h_int)
    tp2= time.clock()
    DDnoptnp[nbr_DD] = tp2-tp1
    # print("Temps non opti : ",tp2-tp1)
    
    
    t1 = time.clock()
    y_plot = conduc_double_sp(lbda1, rho1, cp1, epaisseur1,lbda2, rho2, cp2, epaisseur2, dx1, dx2, dt, T_ext, T_int, Nt, h_ext, h_int)
    t2= time.clock()
    DDnoptsp[nbr_DD] = t2-t1
    
    dx1 = dx1 - 0.0001
    dx2 = dx2 - 0.0001


plt.figure(8)
plt.title("Temps d'éxecution")
plt.plot(XXD,DDopt,label="Temps dgtsv")
plt.plot(XXD,DDnoptnp,label="Temps numpy")
plt.plot(XXD,DDnoptsp,label="Temps scipy")
plt.ylabel("Temps d'exécution (s)")
plt.xlabel("N matrice")
plt.legend()
"=================================================================================="

# plt.figure(7)
# plt.title("Profile de température mur deux couches OPTI")
# # plt.ylim(T_ext -10,T_ext+10)
# plt.ylabel("Température (°C)")
# plt.xlabel("Epaisseur (m)")
# plt.plot(x_plot,opti_plot, label="opti")
# plt.plot(x_plot,y_plot, label="non opti")
# plt.legend()




