import numpy as np
from sympy import symbols, Matrix
# General

w = 800
h = 300

Initial_values_dict = {
    ## RFP Pathwy
    "D_u_f_R" : 21 ,
    "D_b_f_R": 0 ,
    "D_u_e_R" : 0 ,
    "D_b_e_R": 0 ,
    "N_R": 0 ,
    "m_u_f_R" : 0 ,
    "m_b_f_R": 0 ,
    "m_u_e_R" : 0 ,
    "m_b_e_R": 0 ,
    "R_R": 0 ,
    "X_uf_R": 0 ,
    "X_bf_R": 0 ,
    "X_ue_R": 0 ,
    "X_be_R": 0 ,
    "P_R" : 0 ,
    "R_X_R": 0,
    
    ## GFP Pathwy
    "D_u_f_G" : 21 ,
    "D_b_f_G": 0 ,
    "D_u_e_G" : 0 ,
    "D_b_e_G": 0 ,
    "N_G": 0 ,
    "m_u_f_G" : 0 ,
    "m_b_f_G": 0 ,
    "m_u_e_G" : 0 ,
    "m_b_e_G": 0 ,
    "R_G": 0 ,
    "X_uf_G": 0 ,
    "X_bf_G": 0 ,
    "X_ue_G": 0 ,
    "X_be_G": 0 ,
    "P_G" : 0 ,
    "R_X_G": 0,
    
    ## Background Proteins
    "D_u_N" : 1000 ,
    "D_b_N": 0 ,
    "N_N": 0 ,
    "m_u_N" : 0 ,
    "m_b_N" : 0 ,
    "R_N"   : 0 ,
    "X_u_N" : 0 ,
    "X_b_N" : 0 ,
    "P_N"   : 0 ,
    
    ## Resources
    "N" : 150 , # 100000
    "R" : 150,  # 100000
    "X" : 180 ,# 5000
    
    ## 3' miRN
    "D_u_f_3" : 20,
    "D_b_f_3" : 0,
    "D_u_e_3" : 0,
    "D_b_e_3" : 0,
    "N_3"     : 0,
    "m_u_f_3" : 0,
    
    "mi_uf3_R" : 0,
    "mi_bf3_R" : 0,
    "mi_ue3_R" : 0,
    "mi_be3_R" : 0,
    "X_uf3_R" : 0,
    "X_bf3_R" : 0,
    "X_ue3_R" : 0,
    "X_be3_R" : 0,
    "mi_uf3_G" : 0,
    "mi_bf3_G" : 0,
    "mi_ue3_G" : 0,
    "mi_be3_G" : 0,
    "X_uf3_G" : 0,
    "X_bf3_G" : 0,
    "X_ue3_G" : 0,
    "X_be3_G" : 0,
    
    "R_3_R" : 0,
    "R_X3_R" : 0,
    
    ## 5' miRN
    "D_u_f_5"  : 20,
    "D_b_f_5" : 0,
    "D_u_e_5" : 0,
    "D_b_e_5" : 0,
    "N_5" : 0,
    "m_u_f_5"  : 0,
    
    "mi_uf5_R" : 0,
    "mi_ue5_R" : 0,
    "X_uf5_R" : 0,
    "X_ue5_R" : 0,
    "mi_uf5_G" : 0,
    "mi_ue5_G" : 0,
    "X_uf5_G" : 0,
    "X_ue5_G" : 0,
    
    "R_5_R" : 0,
    "R_X5_R" : 0,
}


Rates_values_dict= {
#     "kd_p" : 1 ,
#     "kd_m" : 1 ,
#     "kd"   : 0.12/60 ,
#     "kd_c" : 0.14/60 ,
    
#     "km_p" : (2.63*10)/(60) ,
#     "km_m" : (2.63*1)/(60) ,
#     "km"   : (2.63*10)/(60) ,
#     "km_c" : 2.63/(60*1) ,
    
#     "kbind"  : 0.00001, 
#     "kX_m"   : 0.0001,
#     "kX_3"   : 0.001,
#     "kX_5"   : 0.001,
#     "kdeg"   : 1/60,
#     "lam_P"  : 0.00044/60,
#     "lam_mi" : 0.001/60,
    
    
    "kd_p" : 1 ,
    "kd_m" : 1 ,
    "kd"   : 10/60 ,
    "kd_c" : 10/60 ,
    
#     "kd_p" : 0 ,
#     "kd_m" : 0 ,
#     "kd"   : 0.2/60 ,
#     "kd_c" : 0.2/60 ,
    
    ##
    "km_p" : (2.63*0.1)/60,
    # "km_p" : 1, 2/60,
    "km_m" : (2.63*10)/60 ,
    # "km_m" : 1, 2/60,
    "km"   : (2.63*1)/60 ,
    # "km" : 1, 
    "km_c" : 2.63/(10*60),
    # "km_c" : 0.1, 2.63/60,
    
    "kbind"  : 0.00001, 
    "kX_m"   : 0.006,
    # "kX_m" : 0.00002, 0.00023/60,
    "kX_3"   : 0.008,
    "kX_5"   : 0.008,
    "kdeg"   : 0.1/60,
    # "kdeg" : 0.05, 0.015/60,
    "lam_P"  : 0.00044/60,
    # "lam_P"  : 0.0001, 0.00044/60,
    "lam_mi" : 0.023/60,
    # "lam_mi" : 0.001, 0.023/60,
}



delta1 = 0.1
delta2 = 0.1


Rates_bounds_dict = {
    "kd_p" : [Rates_values_dict["kd_p"] -(delta2*Rates_values_dict["kd_p"]), Rates_values_dict["kd_p"] +(delta2*Rates_values_dict["kd_p"])] ,
    "kd_m" : [Rates_values_dict["kd_m"] -(delta2*Rates_values_dict["kd_m"]), Rates_values_dict["kd_m"] +(delta2*Rates_values_dict["kd_m"])] ,
    "kd"   : [Rates_values_dict["kd"]   -(delta2*Rates_values_dict["kd"]),   Rates_values_dict["kd"]   +(delta2*Rates_values_dict["kd"])] ,
    "kd_c" : [Rates_values_dict["kd_c"] -(delta2*Rates_values_dict["kd_c"]), Rates_values_dict["kd_c"] +(delta2*Rates_values_dict["kd_c"])] ,
    
    "km_p" : [Rates_values_dict["km_p"] -(delta2*Rates_values_dict["km_p"]), Rates_values_dict["km_p"] +(delta2*Rates_values_dict["km_p"])] ,
    "km_m" : [Rates_values_dict["km_m"] -(delta2*Rates_values_dict["km_m"]), Rates_values_dict["km_m"] +(delta2*Rates_values_dict["km_m"])] ,
    "km"   : [Rates_values_dict["km"]   -(delta2*Rates_values_dict["km"]),   Rates_values_dict["km"]   +(delta2*Rates_values_dict["km"])] ,
    "km_c" : [Rates_values_dict["km_c"] -(delta2*Rates_values_dict["km_c"]), Rates_values_dict["km_c"] +(delta2*Rates_values_dict["km_c"])] ,
    
    "kbind"  : [Rates_values_dict["kbind"] -(delta1*Rates_values_dict["kbind"]), Rates_values_dict["kbind"] +(delta1*Rates_values_dict["kbind"])] ,
    "kX_m"   : [Rates_values_dict["kX_m"]  -(delta1*Rates_values_dict["kX_m"]),  Rates_values_dict["kX_m"]  +(delta1*Rates_values_dict["kX_m"])] ,
    "kX_3"   : [Rates_values_dict["kX_3"]  -(delta1*Rates_values_dict["kX_3"]),  Rates_values_dict["kX_3"]  +(delta1*Rates_values_dict["kX_3"])] ,
    "kX_5"   : [Rates_values_dict["kX_5"]  -(delta1*Rates_values_dict["kX_5"]),  Rates_values_dict["kX_5"]  +(delta1*Rates_values_dict["kX_5"])] ,
    "kdeg"   : [Rates_values_dict["kdeg"]  -(delta1*Rates_values_dict["kdeg"]),  Rates_values_dict["kdeg"]  +(delta1*Rates_values_dict["kdeg"])] ,
    "lam_P"  : [Rates_values_dict["lam_P"] -(delta1*Rates_values_dict["lam_P"]), Rates_values_dict["lam_P"] +(delta1*Rates_values_dict["lam_P"])] ,
    "lam_mi" : [Rates_values_dict["lam_mi"] -(delta1*Rates_values_dict["lam_mi"]), Rates_values_dict["lam_mi"] +(delta1*Rates_values_dict["lam_mi"])] ,
}


Rates_to_optimise = {
    "kX_3"   : [Rates_values_dict["kX_3"]  -(delta1*Rates_values_dict["kX_3"]),  Rates_values_dict["kX_3"]  +(delta1*Rates_values_dict["kX_3"])] ,
    "kX_5"   : [Rates_values_dict["kX_5"]  -(delta1*Rates_values_dict["kX_5"]),  Rates_values_dict["kX_5"]  +(delta1*Rates_values_dict["kX_5"])] ,
}
scale = 20

experimental_data = {
    "Control" : {
        "P_G" : list(np.array([0, 1, 1, 1, 1, 1, 1])*scale),
        "P_R" : list(np.array([0, 1, 1, 1, 1, 1, 1])*scale),
    },
    "1TS3" : {
        "P_G" : list(np.array([0, 2.058437718, 2.058437718, 2.058437718, 2.058437718, 2.058437718, 2.058437718])*scale),
        "P_R" : list(np.array([0, 0.320775135, 0.320775135, 0.320775135, 0.320775135, 0.320775135, 0.320775135])*scale),
    },
    "3TS3" : {
        "P_G" : list(np.array([0, 3.470857127, 3.470857127, 3.470857127, 3.470857127, 3.470857127, 3.470857127])*scale),
        "P_R" : list(np.array([0, 0.199718713, 0.199718713, 0.199718713, 0.199718713, 0.199718713, 0.199718713])*scale),
    },
    "1TS5" : {
        "P_G" : list(np.array([0, 1.776078657, 1.776078657, 1.776078657, 1.776078657, 1.776078657, 1.776078657])*scale),
        "P_R" : list(np.array([0, 0.159768743, 0.159768743, 0.159768743, 0.159768743, 0.159768743, 0.159768743])*scale),
    },
    "3TS5" : {
        "P_G" : list(np.array([0, 5.658702805, 5.658702805, 5.658702805, 5.658702805, 5.658702805, 5.658702805])*scale),
        "P_R" : list(np.array([0, 0.060305119, 0.060305119, 0.060305119, 0.060305119, 0.060305119, 0.060305119])*scale),
    },
}


time_points = np.array([0, 4, 20, 21, 22, 23, 24]) # Time in hrs
time_points = time_points*60*60
# kd_p_0 = 1 #
# kd_m_0 = 1 #
# kd_0   = 1 #
# kd_c_0 = 0.01

# km_p_0 = 1 #
# km_m_0 = 1 #
# km_0    = 1 #
# km_c_0 = 0.01

# kX_m_0 = 0.0001 #

# lam_P_0 = 0.0025

# lam_mi_0 = 0.0025
# kbind_0  = 1 #
# kX_3_0   = 0.1
# kX_5_0   = 0.1

# kdeg_0 = 0.1 #
# beta_0 = 0.1

# # Initial concentrations with '#' are estimated
# D_r_0  = 1
# Cd_r_0 = 0
# Nd_r_0 = 0
# m_r_0  = 0
# Cm_r_0 = 0
# Nm_r_0 = 0
# Ym_r_0 = 0
# YC_r_0 = 0
# Bm_r_0 = 0
# P_r_0  = 0

# D_g_0  = 1
# Cd_g_0 = 0
# Nd_g_0 = 0
# m_g_0  = 0
# Cm_g_0 = 0
# Nm_g_0 = 0
# Ym_g_0 = 0
# YC_g_0 = 0
# Bm_g_0 = 0
# P_g_0  = 0

# D_n_0  = 100
# Cd_n_0 = 0
# Nd_n_0 = 0
# m_n_0  = 0
# Cm_n_0 = 0
# Nm_n_0 = 0
# Ym_n_0 = 0
# YC_n_0 = 0
# Bm_n_0 = 0
# P_n_0  = 0

# Pol_0  = 100000
# R_0    = 100000
# X_0    = 50000 

# D_3_0  = 1
# Cd_3_0 = 0
# Nd_3_0 = 0
# mi3_0  = 0
# L3C_r_0 = 0
# Y3C_r_0 = 0
# L3m_r_0 = 0
# Y3m_r_0 = 0

# D_5_0  = 1
# Cd_5_0 = 0
# Nd_5_0 = 0
# mi5_0  = 0
# L5_0   = 0
# Y5_0   = 0

# #   Rates with '#' are estimated 
# kd_p_0 = 1 #
# kd_m_0 = 1 #
# kd_0   = 1 #
# kd_c_0 = 0.01

# km_p_0 = 1 #
# km_m_0 = 1 #
# km_0    = 1 #
# km_c_0 = 0.01

# kX_m_0 = 0.0001 #

# lam_P_0 = 0.0025

# lam_mi_0 = 0.0025
# kbind_0  = 1 #
# kX_3_0   = 0.1
# kX_5_0   = 0.1

# kdeg_0 = 0.1 #
# beta_0 = 0.1

# initial_values_noregulation = [D_r_0, Cd_r_0, Nd_r_0, m_r_0, Cm_r_0, Nm_r_0, Ym_r_0, YC_r_0, Bm_r_0, P_r_0, D_g_0, Cd_g_0, Nd_g_0, m_g_0, Cm_g_0, Nm_g_0, Ym_g_0, YC_g_0, Bm_g_0, P_g_0, D_n_0, Cd_n_0, Nd_n_0, m_n_0, Cm_n_0, Nm_n_0, Ym_n_0, YC_n_0, Bm_n_0, P_n_0, Pol_0, R_0, X_0]
# args_noregulation =  km_p_0, km_m_0, km_0, km_c_0, kd_p_0, kd_m_0, kd_0, kd_c_0, kX_m_0, kdeg_0, lam_P_0


# initial_values_noregulation_SIMP = [D_r_0, Cd_r_0, Nd_r_0, m_r_0, Cm_r_0, Nm_r_0, Ym_r_0, YC_r_0, P_r_0, D_g_0, Cd_g_0, Nd_g_0, m_g_0, Cm_g_0, Nm_g_0, Ym_g_0, YC_g_0, P_g_0, D_n_0, Cd_n_0, Nd_n_0, m_n_0, Cm_n_0, Nm_n_0, Ym_n_0, YC_n_0, P_n_0, Pol_0, R_0, X_0]
# args_noregulation_SIMP =  km_p_0, km_m_0, km_0, km_c_0, kd_p_0, kd_m_0, kd_0, kd_c_0, kX_m_0, kdeg_0, lam_P_0


# initial_values_threeregulation = [D_r_0, Cd_r_0, Nd_r_0, m_r_0, Cm_r_0, Nm_r_0, Ym_r_0, YC_r_0, Bm_r_0, P_r_0, D_g_0, Cd_g_0, Nd_g_0, m_g_0, Cm_g_0, Nm_g_0, Ym_g_0, YC_g_0, Bm_g_0, P_g_0, D_n_0, Cd_n_0, Nd_n_0, m_n_0, Cm_n_0, Nm_n_0, Ym_n_0, YC_n_0, Bm_n_0, P_n_0, Pol_0, R_0, X_0, D_3_0, Cd_3_0, Nd_3_0, mi3_0, L3m_r_0, L3C_r_0, Y3m_r_0, Y3C_r_0, B3_r_0, K3_r_0]
# args_threeregulation =  km_p_0, km_m_0, km_0, km_c_0, kd_p_0, kd_m_0, kd_0, kd_c_0, kbind_0, kX_3_0, kX_m_0, kdeg_0, lam_P_0, lam_mi_0

# initial_values_threeregulation_SIMP = [D_r_0, Cd_r_0, Nd_r_0, m_r_0, Cm_r_0, Nm_r_0, Ym_r_0, YC_r_0, P_r_0, D_g_0, Cd_g_0, Nd_g_0, m_g_0, Cm_g_0, Nm_g_0, Ym_g_0, YC_g_0, P_g_0, D_n_0, Cd_n_0, Nd_n_0, m_n_0, Cm_n_0, Nm_n_0, Ym_n_0, YC_n_0, P_n_0, Pol_0, R_0, X_0, D_3_0, Cd_3_0, Nd_3_0, mi3_0, L3m_r_0, L3C_r_0, Y3m_r_0, Y3C_r_0]
# args_threeregulation_SIMP =  km_p_0, km_m_0, km_0, km_c_0, kd_p_0, kd_m_0, kd_0, kd_c_0, kbind_0, kX_3_0, kX_m_0, kdeg_0, lam_P_0, lam_mi_0


# S_D_r_0  = 0
# S_Cd_r_0 = 0
# S_Nd_r_0 = 0
# S_m_r_0  = 0
# S_Cm_r_0 = 0
# S_Nm_r_0 = 0
# S_Ym_r_0 = 0
# S_YC_r_0 = 0
# S_Bm_r_0 = 0
# S_P_r_0  = 0

# S_D_g_0  = 0
# S_Cd_g_0 = 0
# S_Nd_g_0 = 0
# S_m_g_0  = 0
# S_Cm_g_0 = 0
# S_Nm_g_0 = 0
# S_Ym_g_0 = 0
# S_YC_g_0 = 0
# S_Bm_g_0 = 0
# S_P_g_0  = 0

# S_D_n_0  = 0
# S_Cd_n_0 = 0
# S_Nd_n_0 = 0
# S_m_n_0  = 0
# S_Cm_n_0 = 0
# S_Nm_n_0 = 0
# S_Ym_n_0 = 0
# S_YC_n_0 = 0
# S_Bm_n_0 = 0
# S_P_n_0  = 0

# S_Pol_0 = 0
# S_R_0   = 0
# S_X_0   = 0

# S_D_3_0   = 0
# S_Cd_3_0  = 0
# S_Nd_3_0  = 0
# S_mi3_0   = 0
# S_L3C_r_0 = 0
# S_Y3C_r_0 = 0
# S_L3m_r_0 = 0
# S_Y3m_r_0 = 0

# S_B3_r_0 = 0
# S_K3_r_0 = 0

# S_D_5_0  = 0
# S_Cd_5_0 = 0
# S_Nd_5_0 = 0
# S_mi5_0  = 0
# S_L5_0   = 0
# S_Y5_0   = 0

# initial_values_sensitivity_noregulation = [D_r_0, Cd_r_0, Nd_r_0, m_r_0, Cm_r_0, Nm_r_0, Ym_r_0, YC_r_0, Bm_r_0, P_r_0, D_g_0, Cd_g_0, Nd_g_0, m_g_0, Cm_g_0, Nm_g_0, Ym_g_0, YC_g_0, Bm_g_0, P_g_0, D_n_0, Cd_n_0, Nd_n_0, m_n_0, Cm_n_0, Nm_n_0, Ym_n_0, YC_n_0, Bm_n_0, P_n_0, Pol_0, R_0, X_0, S_D_r_0, S_Cd_r_0, S_Nd_r_0, S_m_r_0, S_Cm_r_0, S_Nm_r_0, S_Ym_r_0, S_YC_r_0, S_Bm_r_0, S_P_r_0, S_D_g_0, S_Cd_g_0, S_Nd_g_0, S_m_g_0, S_Cm_g_0, S_Nm_g_0, S_Ym_g_0, S_YC_g_0, S_Bm_g_0, S_P_g_0, S_D_n_0, S_Cd_n_0, S_Nd_n_0, S_m_n_0, S_Cm_n_0, S_Nm_n_0, S_Ym_n_0, S_YC_n_0, S_Bm_n_0, S_P_n_0, S_Pol_0, S_R_0, S_X_0]

# initial_values_sensitivity_noregulation_SIMP = [D_r_0, Cd_r_0, Nd_r_0, m_r_0, Cm_r_0, Nm_r_0, Ym_r_0, YC_r_0, P_r_0, D_g_0, Cd_g_0, Nd_g_0, m_g_0, Cm_g_0, Nm_g_0, Ym_g_0, YC_g_0, P_g_0, D_n_0, Cd_n_0, Nd_n_0, m_n_0, Cm_n_0, Nm_n_0, Ym_n_0, YC_n_0, P_n_0, Pol_0, R_0, X_0, S_D_r_0, S_Cd_r_0, S_Nd_r_0, S_m_r_0, S_Cm_r_0, S_Nm_r_0, S_Ym_r_0, S_YC_r_0, S_P_r_0, S_D_g_0, S_Cd_g_0, S_Nd_g_0, S_m_g_0, S_Cm_g_0, S_Nm_g_0, S_Ym_g_0, S_YC_g_0, S_P_g_0, S_D_n_0, S_Cd_n_0, S_Nd_n_0, S_m_n_0, S_Cm_n_0, S_Nm_n_0, S_Ym_n_0, S_YC_n_0, S_P_n_0, S_Pol_0, S_R_0, S_X_0]

# initial_values_sensitivity_threeregulation = [D_r_0, Cd_r_0, Nd_r_0, m_r_0, Cm_r_0, Nm_r_0, Ym_r_0, YC_r_0, Bm_r_0, P_r_0, D_g_0, Cd_g_0, Nd_g_0, m_g_0, Cm_g_0, Nm_g_0, Ym_g_0, YC_g_0, Bm_g_0, P_g_0, D_n_0, Cd_n_0, Nd_n_0, m_n_0, Cm_n_0, Nm_n_0, Ym_n_0, YC_n_0, Bm_n_0, P_n_0, Pol_0, R_0, X_0, D_3_0, Cd_3_0, Nd_3_0, mi3_0, L3m_r_0, L3C_r_0, Y3m_r_0, Y3C_r_0, B3_r_0, K3_r_0, S_D_r_0, S_Cd_r_0, S_Nd_r_0, S_m_r_0, S_Cm_r_0, S_Nm_r_0, S_Ym_r_0, S_YC_r_0, S_Bm_r_0, S_P_r_0, S_D_g_0, S_Cd_g_0, S_Nd_g_0, S_m_g_0, S_Cm_g_0, S_Nm_g_0, S_Ym_g_0, S_YC_g_0, S_Bm_g_0, S_P_g_0, S_D_n_0, S_Cd_n_0, S_Nd_n_0, S_m_n_0, S_Cm_n_0, S_Nm_n_0, S_Ym_n_0, S_YC_n_0, S_Bm_n_0, S_P_n_0, S_Pol_0, S_R_0, S_X_0, S_D_3_0, S_Cd_3_0, S_Nd_3_0, S_mi3_0, S_L3m_r_0, S_L3C_r_0, S_Y3m_r_0, S_Y3C_r_0, S_B3_r_0, S_K3_r_0]

# initial_values_sensitivity_threeregulation_SIMP = [D_r_0, Cd_r_0, Nd_r_0, m_r_0, Cm_r_0, Nm_r_0, Ym_r_0, YC_r_0, P_r_0, D_g_0, Cd_g_0, Nd_g_0, m_g_0, Cm_g_0, Nm_g_0, Ym_g_0, YC_g_0, P_g_0, D_n_0, Cd_n_0, Nd_n_0, m_n_0, Cm_n_0, Nm_n_0, Ym_n_0, YC_n_0, P_n_0, Pol_0, R_0, X_0, D_3_0, Cd_3_0, Nd_3_0, mi3_0, L3m_r_0, L3C_r_0, Y3m_r_0, Y3C_r_0, S_D_r_0, S_Cd_r_0, S_Nd_r_0, S_m_r_0, S_Cm_r_0, S_Nm_r_0, S_Ym_r_0, S_YC_r_0, S_P_r_0, S_D_g_0, S_Cd_g_0, S_Nd_g_0, S_m_g_0, S_Cm_g_0, S_Nm_g_0, S_Ym_g_0, S_YC_g_0, S_P_g_0, S_D_n_0, S_Cd_n_0, S_Nd_n_0, S_m_n_0, S_Cm_n_0, S_Nm_n_0, S_Ym_n_0, S_YC_n_0, S_P_n_0, S_Pol_0, S_R_0, S_X_0, S_D_3_0, S_Cd_3_0, S_Nd_3_0, S_mi3_0, S_L3m_r_0, S_L3C_r_0, S_Y3m_r_0, S_Y3C_r_0]


# sym_D_r, sym_Cd_r, sym_Nd_r, sym_m_r, sym_Ym_r, sym_YC_r, sym_B_r, sym_Cm_r, sym_Nm_r, sym_P_r  = symbols('D_r, Cd_r, Nd_r, m_r, Ym_r, YC_r, B_r, Cm_r, Nm_r, P_r')
# sym_D_g, sym_Cd_g, sym_Nd_g, sym_m_g, sym_Ym_g, sym_YC_g, sym_B_g, sym_Cm_g, sym_Nm_g, sym_P_g  = symbols('D_g, Cd_g, Nd_g, m_g, Ym_g, YC_g, B_g, Cm_g, Nm_g, P_g')
# sym_D_n, sym_Cd_n, sym_Nd_n, sym_m_n, sym_Ym_n, sym_YC_n, sym_B_n, sym_Cm_n, sym_Nm_n, sym_P_n  = symbols('D_n, Cd_n, Nd_n, m_n, Ym_n, YC_n, B_n, Cm_n, Nm_n, P_n')
# sym_D_3, sym_Cd_3, sym_Nd_3, sym_mi3, sym_K_3, sym_L3C_r, sym_L3m_r, sym_Y3C_r, sym_Y3m_r, sym_B3_r, sym_K3_r = symbols('D_3, Cd_3, Nd_3, mi3, K3_r, L3C_r, L3m_r, Y3C_r, Y3m_r, B3_r, K3_r')
# sym_D_5, sym_Cd_5, sym_Nd_5, sym_mi5, sym_K_5, sym_L_5, sym_Y_5 = symbols('D_5, Cd_5, Nd_5, mi5, K_5, L_5, Y_5')
# sym_Pol, sym_R, sym_X = symbols('Pol, R, X')

# S_sym_D_r, S_sym_Cd_r, S_sym_Nd_r, S_sym_m_r, S_sym_Ym_r, S_sym_YC_r, S_sym_B_r, S_sym_Cm_r, S_sym_Nm_r, S_sym_P_r  = symbols('S_D_r, S_Cd_r, S_Nd_r, S_m_r, S_Ym_r, S_YC_r, S_B_r, S_Cm_r, S_Nm_r, S_P_r')
# S_sym_D_g, S_sym_Cd_g, S_sym_Nd_g, S_sym_m_g, S_sym_Ym_g, S_sym_YC_g, S_sym_B_g, S_sym_Cm_g, S_sym_Nm_g, S_sym_P_g  = symbols('S_D_g, S_Cd_g, S_Nd_g, S_m_g, S_Ym_g, S_YC_g, S_B_g, S_Cm_g, S_Nm_g, S_P_g')
# S_sym_D_n, S_sym_Cd_n, S_sym_Nd_n, S_sym_m_n, S_sym_Ym_n, S_sym_YC_n, S_sym_B_n, S_sym_Cm_n, S_sym_Nm_n, S_sym_P_n  = symbols('S_D_n, S_Cd_n, S_Nd_n, S_m_n, S_Ym_n, S_YC_n, S_B_n, S_Cm_n, S_Nm_n, S_P_n')
# S_sym_D_3, S_sym_Cd_3, S_sym_Nd_3, S_sym_mi3, S_sym_K_3, S_sym_L3C_r, S_sym_L3m_r, S_sym_Y3C_r, S_sym_Y3m_r, S_sym_B3_r, S_sym_K3_r = symbols('S_D_3, S_Cd_3, S_Nd_3, S_mi3, S_K_3, S_L3C_r, S_L3m_r, S_Y3C_r, S_Y3m_r, S_B3_r, S_K3_r')
# S_sym_D_5, S_sym_Cd_5, S_sym_Nd_5, S_sym_mi5, S_sym_K_5, S_sym_L_5, S_sym_Y_5 = symbols('S_D_5, S_Cd_5, S_Nd_5, S_mi5, S_K_5, S_L_5, S_Y_5')
# S_sym_Pol, S_sym_R, S_sym_X = symbols('S_Pol, S_R, S_X')

# sym_km_p, sym_km_m, sym_km, sym_km_c, sym_kd_p, sym_kd_m, sym_kd, sym_kd_c, sym_kX_m, sym_lam_P, sym_lam_mi, sym_kbind, sym_kdeg, sym_kX_3, sym_kX_5 = symbols("km_p, km_m, km, km_c, kd_p, kd_m, kd, kd_c, kX_m, lam_P, lam_mi, kbind, kdeg, kX_3, kX_5")

# # No Regulation
# sym_k_noreg = sym_km_p, sym_km_m, sym_km, sym_km_c, sym_kd_p, sym_kd_m, sym_kd, sym_kd_c, sym_kX_m, sym_kdeg, sym_lam_P
# sym_y_noreg = sym_D_r, sym_Cd_r, sym_Nd_r, sym_m_r, sym_Cm_r, sym_Nm_r, sym_Ym_r, sym_YC_r, sym_B_r, sym_P_r, sym_D_g, sym_Cd_g, sym_Nd_g, sym_m_g, sym_Cm_g, sym_Nm_g, sym_Ym_g, sym_YC_g, sym_B_g, sym_P_g, sym_D_n, sym_Cd_n, sym_Nd_n, sym_m_n, sym_Cm_n, sym_Nm_n, sym_Ym_n, sym_YC_n, sym_B_n, sym_P_n, sym_Pol, sym_R, sym_X
# S_sym_y_noreg = (S_sym_D_r, S_sym_Cd_r, S_sym_Nd_r, S_sym_m_r, S_sym_Cm_r, S_sym_Nm_r, S_sym_Ym_r, S_sym_YC_r, S_sym_B_r, S_sym_P_r, S_sym_D_g, S_sym_Cd_g, S_sym_Nd_g, S_sym_m_g, S_sym_Cm_g, S_sym_Nm_g, S_sym_Ym_g, S_sym_YC_g, S_sym_B_g, S_sym_P_g, S_sym_D_n, S_sym_Cd_n, S_sym_Nd_n, S_sym_m_n, S_sym_Ym_n, S_sym_YC_n, S_sym_B_n, S_sym_Cm_n, S_sym_Nm_n, S_sym_P_n, S_sym_Pol, S_sym_R, S_sym_X)

# y_noreg = Matrix(list(sym_y_noreg))
# S_noreg = Matrix(list(S_sym_y_noreg))
# P_noreg = list(sym_k_noreg)
# rates_vals_noreg = km_p_0, km_m_0, km_0, km_c_0, kd_p_0, kd_m_0, kd_0, kd_c_0, kX_m_0, kdeg_0, lam_P_0

# # No Regulation SIMPLIFIED
# sym_k_noreg_SIMP = sym_km_p, sym_km_m, sym_km, sym_km_c, sym_kd_p, sym_kd_m, sym_kd, sym_kd_c, sym_kX_m, sym_kdeg, sym_lam_P
# sym_y_noreg_SIMP = sym_D_r, sym_Cd_r, sym_Nd_r, sym_m_r, sym_Cm_r, sym_Nm_r, sym_Ym_r, sym_YC_r, sym_P_r, sym_D_g, sym_Cd_g, sym_Nd_g, sym_m_g, sym_Cm_g, sym_Nm_g, sym_Ym_g, sym_YC_g, sym_P_g, sym_D_n, sym_Cd_n, sym_Nd_n, sym_m_n, sym_Cm_n, sym_Nm_n, sym_Ym_n, sym_YC_n, sym_P_n, sym_Pol, sym_R, sym_X
# S_sym_y_noreg_SIMP = (S_sym_D_r, S_sym_Cd_r, S_sym_Nd_r, S_sym_m_r, S_sym_Cm_r, S_sym_Nm_r, S_sym_Ym_r, S_sym_YC_r, S_sym_P_r, S_sym_D_g, S_sym_Cd_g, S_sym_Nd_g, S_sym_m_g, S_sym_Cm_g, S_sym_Nm_g, S_sym_Ym_g, S_sym_YC_g, S_sym_P_g, S_sym_D_n, S_sym_Cd_n, S_sym_Nd_n, S_sym_m_n, S_sym_Ym_n, S_sym_YC_n, S_sym_Cm_n, S_sym_Nm_n, S_sym_P_n, S_sym_Pol, S_sym_R, S_sym_X)

# y_noreg_SIMP = Matrix(list(sym_y_noreg_SIMP))
# S_noreg_SIMP = Matrix(list(S_sym_y_noreg_SIMP))
# P_noreg_SIMP = list(sym_k_noreg_SIMP)
# rates_vals_noreg_SIMP = km_p_0, km_m_0, km_0, km_c_0, kd_p_0, kd_m_0, kd_0, kd_c_0, kX_m_0, kdeg_0, lam_P_0

# # 3' Regulation
# sym_k_threereg = sym_km_p, sym_km_m, sym_km, sym_km_c, sym_kd_p, sym_kd_m, sym_kd, sym_kd_c, sym_kbind, sym_kX_3, sym_kX_m, sym_kdeg, sym_lam_P, sym_lam_mi
# sym_y_threereg = sym_D_r, sym_Cd_r, sym_Nd_r, sym_m_r, sym_Cm_r, sym_Nm_r, sym_Ym_r, sym_YC_r, sym_B_r, sym_P_r, sym_D_g, sym_Cd_g, sym_Nd_g, sym_m_g, sym_Cm_g, sym_Nm_g, sym_Ym_g, sym_YC_g, sym_B_g, sym_P_g, sym_D_n, sym_Cd_n, sym_Nd_n, sym_m_n, sym_Cm_n, sym_Nm_n, sym_Ym_n, sym_YC_n, sym_B_n, sym_P_n, sym_Pol, sym_R, sym_X, sym_D_3, sym_Cd_3, sym_Nd_3, sym_mi3, sym_L3m_r, sym_L3C_r, sym_Y3m_r, sym_Y3C_r, sym_B3_r, sym_K3_r
# S_sym_y_threereg = (S_sym_D_r, S_sym_Cd_r, S_sym_Nd_r, S_sym_m_r, S_sym_Cm_r, S_sym_Nm_r, S_sym_Ym_r, S_sym_YC_r, S_sym_B_r, S_sym_P_r, S_sym_D_g, S_sym_Cd_g, S_sym_Nd_g, S_sym_m_g, S_sym_Cm_g, S_sym_Nm_g, S_sym_Ym_g, S_sym_YC_g, S_sym_B_g, S_sym_P_g, S_sym_D_n, S_sym_Cd_n, S_sym_Nd_n, S_sym_m_n, S_sym_Ym_n, S_sym_YC_n, S_sym_B_n, S_sym_Cm_n, S_sym_Nm_n, S_sym_P_n, S_sym_Pol, S_sym_R, S_sym_X, S_sym_D_3, S_sym_Cd_3, S_sym_Nd_3, S_sym_mi3, S_sym_L3m_r, S_sym_L3C_r, S_sym_Y3m_r, S_sym_Y3C_r, S_sym_B3_r, S_sym_K3_r)

# y_threereg = Matrix(list(sym_y_threereg))
# S_threereg = Matrix(list(S_sym_y_threereg))
# P_threereg = list(sym_k_threereg)
# rates_vals_threereg = km_p_0, km_m_0, km_0, km_c_0, kd_p_0, kd_m_0, kd_0, kd_c_0, kbind_0, kX_3_0, kX_m_0, kdeg_0, lam_P_0, lam_mi_0

# # 3' Regulation SIMP
# sym_k_threereg_SIMP = sym_km_p, sym_km_m, sym_km, sym_km_c, sym_kd_p, sym_kd_m, sym_kd, sym_kd_c, sym_kbind, sym_kX_3, sym_kX_m, sym_kdeg, sym_lam_P, sym_lam_mi
# sym_y_threereg_SIMP = sym_D_r, sym_Cd_r, sym_Nd_r, sym_m_r, sym_Cm_r, sym_Nm_r, sym_Ym_r, sym_YC_r, sym_P_r, sym_D_g, sym_Cd_g, sym_Nd_g, sym_m_g, sym_Cm_g, sym_Nm_g, sym_Ym_g, sym_YC_g, sym_P_g, sym_D_n, sym_Cd_n, sym_Nd_n, sym_m_n, sym_Cm_n, sym_Nm_n, sym_Ym_n, sym_YC_n, sym_P_n, sym_Pol, sym_R, sym_X, sym_D_3, sym_Cd_3, sym_Nd_3, sym_mi3, sym_L3m_r, sym_L3C_r, sym_Y3m_r, sym_Y3C_r
# S_sym_y_threereg_SIMP = (S_sym_D_r, S_sym_Cd_r, S_sym_Nd_r, S_sym_m_r, S_sym_Cm_r, S_sym_Nm_r, S_sym_Ym_r, S_sym_YC_r, S_sym_P_r, S_sym_D_g, S_sym_Cd_g, S_sym_Nd_g, S_sym_m_g, S_sym_Cm_g, S_sym_Nm_g, S_sym_Ym_g, S_sym_YC_g, S_sym_P_g, S_sym_D_n, S_sym_Cd_n, S_sym_Nd_n, S_sym_m_n, S_sym_Ym_n, S_sym_YC_n, S_sym_Cm_n, S_sym_Nm_n, S_sym_P_n, S_sym_Pol, S_sym_R, S_sym_X, S_sym_D_3, S_sym_Cd_3, S_sym_Nd_3, S_sym_mi3, S_sym_L3m_r, S_sym_L3C_r, S_sym_Y3m_r, S_sym_Y3C_r)

# y_threereg_SIMP = Matrix(list(sym_y_threereg_SIMP))
# S_threereg_SIMP = Matrix(list(S_sym_y_threereg_SIMP))
# P_threereg_SIMP = list(sym_k_threereg_SIMP)
# rates_vals_threereg_SIMP = km_p_0, km_m_0, km_0, km_c_0, kd_p_0, kd_m_0, kd_0, kd_c_0, kbind_0, kX_3_0, kX_m_0, kdeg_0, lam_P_0, lam_mi_0