"=======================Fonctions calculant les températures intérieurs======================="

"==========Température intérieure=========="
def temp_int(h_int,V_int,S_sud,S_est,S_nord,S_ouest, S_plafond, S_sol, T_int, T_sud, T_est, T_nord, T_ouest, T_plafond, T_sol, dt):
    rho_air = 1.204 #Masse volumique air (20 °C)
    cp_air = 1256 #Capacité thermique air (m−3.K)
    
    g_sud = (h_int * S_sud * dt) / (V_int * rho_air * cp_air)
    g_est = (h_int * S_est * dt) / (V_int * rho_air * cp_air)
    g_nord = (h_int * S_nord * dt) / (V_int * rho_air * cp_air)
    g_ouest = (h_int * S_ouest * dt) / (V_int * rho_air * cp_air)
    g_plafond = (h_int * S_plafond * dt) / (V_int * rho_air * cp_air)
    g_sol = (h_int * S_sol * dt) / (V_int * rho_air * cp_air)
    
    n_sud = len(T_sud)-1
    n_est = len(T_est)-1
    n_nord = len(T_nord)-1
    n_ouest = len(T_ouest)-1
    n_plafond = len(T_plafond)-1
    n_sol = len(T_sol)-1
    
    temp = T_int + g_sud*(T_sud[n_sud] - T_int) + g_est*(T_est[n_est] - T_int) + g_nord*(T_nord[n_nord] - T_int) + g_ouest*(T_ouest[n_ouest] - T_int) + g_plafond*(T_plafond[n_plafond] - T_int) + g_sol*(T_sol[n_sol] - T_int)
    
    temp = 1 / (1 + g_sud + g_est + g_nord + g_ouest + g_plafond + g_sol) * (T_int + g_sud * T_sud[n_sud] + g_est*T_est[n_est] + g_nord*T_nord[n_nord] + g_ouest*T_ouest[n_ouest] + g_plafond*T_plafond[n_plafond] + g_sol*T_sol[n_sol])
    
    return temp


def puissance_chauffe(T_ext,T_int,T_need,V_int):
    roh_air = 1.204 #Masse volumique air (20 °C)
    Cp_air = 1256 #Capacité thermique air (J.m−3.K)
    cp_air_wh = 0.3489 #Capacité thermique air en Wh.m-3.K
    if T_int < T_need:
        # P_use = roh_air * V_int * Cp_air * (T_need - T_int)  /( 60*60*1000)
        P_use = V_int * cp_air_wh * roh_air * (T_need - T_int)  / 1000
    elif T_int >= T_need:
        P_use = 0
    return P_use


def renouvellement_air(T_ext,T_int,V_int,V_ventil_horaire):
    temp = (V_ventil_horaire/4 * T_ext + V_int * T_int) / (V_ventil_horaire/4 + V_int)
    
    
    return temp

