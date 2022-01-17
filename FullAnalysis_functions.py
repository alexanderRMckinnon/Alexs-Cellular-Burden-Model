import numpy as np

def transcription_or_translation(target_u_f, target_b_f, target_u_e, target_b_e, machinery, elongation, product, kp, km, k, kc):
    kc = kc*(target_u_e+target_b_e)
    dtarget_u_f_dt = -target_u_f*machinery*kp +target_b_f*km  
    dtarget_b_f_dt =  target_u_f*machinery*kp -target_b_f*km                                         -target_b_f*k
    dtarget_u_e_dt =                                         -target_u_e*machinery*kp +target_b_e*km +target_b_f*k +target_b_e*k
    dtarget_b_e_dt =                                          target_u_e*machinery*kp -target_b_e*km               -target_b_e*k
    dmachinery_dt  = -target_u_f*machinery*kp +target_b_f*km -target_u_e*machinery*kp +target_b_e*km                             +elongation*kc
    delongation_dt =                                                                                  target_b_f*k +target_b_e*k -elongation*kc
    dproduct_dt =                                                                                                                 elongation*kc
    return np.array([dtarget_u_f_dt, dtarget_b_f_dt, dtarget_u_e_dt, dtarget_b_e_dt, dmachinery_dt, delongation_dt, dproduct_dt])

def scale_calc(m, m_inv, R):
    return R/(m+m_inv +(1*10**(-50)))

def mRNA_3miRNAdegradation_fullpathway(m_u_f, m_b_f, m_u_e, m_b_e, mi_uf3, mi_bf3, mi_ue3, mi_be3, X_uf3, X_bf3, X_ue3, X_be3, m_3, X, R, R_j, R_3j, R_X3j, P, kbind, kX3, kdeg, kp, km, k, kc, TS):
    
    ## m_u_f + m_3 ->[k] mi_uf3
    dm_u_f_dt, dm_3_dt_1 , dR_j_dt_1, dmi_uf3_dt_1 , dR_3j_dt_1  = mRNA_binding(m_u_f, m_3, 0, kbind*TS, mi_uf3)
    ## m_b_f + m_3 ->[k] mi_bf3
    dm_b_f_dt, dm_3_dt_2 , dR_j_dt_2, dmi_bf3_dt_1 , dR_3j_dt_2  = mRNA_binding(m_b_f, m_3, 0, kbind*TS, mi_bf3)
    ## m_u_e + m_3 + {scale*R_j} ->[k] mi_ue3 + {scale*R_3j}
    dm_u_e_dt, dm_3_dt_3 , dR_j_dt_3, dmi_ue3_dt_1 , dR_3j_dt_3  = mRNA_binding(m_u_e, m_3, scale_calc(m_u_e, m_b_e, R_j), kbind*TS, mi_ue3)
    ## m_b_e + m_3 + {scale*R_j} ->[k] mi_be3 + {scale*R_3j}
    dm_b_e_dt, dm_3_dt_4 , dR_j_dt_4, dmi_be3_dt_1 , dR_3j_dt_4  = mRNA_binding(m_b_e, m_3, scale_calc(m_b_e, m_u_e, R_j), kbind*TS, mi_be3)
    
    
    ## mi_uf3 + X ->[k] X_uf3
    dmi_uf3_dt_2, dX_dt_1 , dR_3j_dt_5, dX_uf3_dt_1 , dR_X3j_dt_1  = mRNA_binding(mi_uf3, X, 0, kX3, 0)
    ## mi_bf3 + X ->[k] X_bf3
    dmi_bf3_dt_2, dX_dt_2 , dR_3j_dt_6, dX_bf3_dt_1 , dR_X3j_dt_2  = mRNA_binding(mi_bf3, X, 0, kX3, 0)
    ## mi_ue3 + X + {scale*R_3j} ->[k] X_ue3 + {scale*R_X3j}
    dmi_ue3_dt_2, dX_dt_3 , dR_3j_dt_7, dX_ue3_dt_1 , dR_X3j_dt_3  = mRNA_binding(mi_ue3, X, scale_calc(mi_ue3, mi_be3, R_3j), kX3, 0)
    ## mi_be3 + X + {scale*R_3j} ->[k] X_be3 + {scale*R_X3j}
    dmi_be3_dt_2, dX_dt_4 , dR_3j_dt_8, dX_be3_dt_1 , dR_X3j_dt_4  = mRNA_binding(mi_be3, X, scale_calc(mi_be3, mi_ue3, R_3j), kX3, 0)
    
    ## mi_uf3 + R <=> mi_bf3 -> mi_ue3 + R_j3
    ## mi_ue3 + R <=> mi_be3 -> mi_ue3 + R_j3
    ## R_j3 -> R + P
    dmi_uf3_dt_3, dmi_bf3_dt_3, dmi_ue3_dt_3, dmi_be3_dt_3, dR_dt_1, dR_3j_dt_9, dP_dt = transcription_or_translation(mi_uf3, mi_bf3, mi_ue3, mi_be3, R, R_3j, P, kp, km, k, kc)
    
    ## X_uf3 ->[k] X 
    dX_uf3_dt_2, dR_X3j_dt_5, dX_dt_5, dR_dt_2 =  mRNA_degradation(X_uf3, 0, kdeg, 0)
    ## X_bf3 ->[k] X + (1)*Ribosome
    dX_bf3_dt_2, dR_X3j_dt_6, dX_dt_6, dR_dt_3 =  mRNA_degradation(X_bf3, 0, kdeg, 1)
    ## X_ue3 + {scale*R_X3j}->[k] X + scale*Ribosomes
    dX_ue3_dt_2, dR_X3j_dt_7, dX_dt_7, dR_dt_4 =  mRNA_degradation(X_ue3, scale_calc(X_ue3, X_be3, R_X3j), kdeg, 0)
    ## X_be3 + {scale*R_X3j}->[k] X + (scale+1)*Ribosomes
    dX_be3_dt_2, dR_X3j_dt_8, dX_dt_8, dR_dt_5 =  mRNA_degradation(X_be3, scale_calc(X_be3, X_ue3, R_X3j), kdeg, 1)
    
    return np.array([dm_u_f_dt, dm_b_f_dt, dm_u_e_dt, dm_b_e_dt, dmi_uf3_dt_1+dmi_uf3_dt_2+dmi_uf3_dt_3, dmi_bf3_dt_1+dmi_bf3_dt_2+dmi_bf3_dt_3,  dmi_ue3_dt_1+dmi_ue3_dt_2+dmi_ue3_dt_3, dmi_be3_dt_1+dmi_be3_dt_2+dmi_be3_dt_3, dX_uf3_dt_1+dX_uf3_dt_2, dX_bf3_dt_1+dX_bf3_dt_2, dX_ue3_dt_1+dX_ue3_dt_2, dX_be3_dt_1+dX_be3_dt_2, dm_3_dt_1+dm_3_dt_2+dm_3_dt_3+dm_3_dt_4, dX_dt_1+dX_dt_2+dX_dt_3+dX_dt_4+dX_dt_5+dX_dt_6+dX_dt_7+dX_dt_8, dR_dt_1+dR_dt_2+dR_dt_3+dR_dt_4+dR_dt_5, dR_j_dt_1+dR_j_dt_2+dR_j_dt_3+dR_j_dt_4, dR_3j_dt_1+dR_3j_dt_2+dR_3j_dt_3+dR_3j_dt_4+dR_3j_dt_5+dR_3j_dt_6+dR_3j_dt_7+dR_3j_dt_8+dR_3j_dt_9, dR_X3j_dt_1+dR_X3j_dt_2+dR_X3j_dt_3+dR_X3j_dt_4+dR_X3j_dt_5+dR_X3j_dt_6+dR_X3j_dt_7+dR_X3j_dt_8, dP_dt])


def mRNA_5miRNAdegradation_fullpathway(m_u_f, m_u_e, m_b_e, mi_uf5, mi_ue5, X_uf5, X_ue5, m_5, X, R, R_j, R_5j, R_X5j, P, kbind, kX5, kdeg, kc, TS):
    ## m_u_f + m_5 ->[k] mi_uf5 
    dm_u_f_dt, dm_5_dt_1, dR_j_dt_1, dmi_uf5_dt_1 , dR_5j_dt_1  = mRNA_binding(m_u_f, m_5, 0, kbind*TS, mi_uf5)
    ## m_u_e + m_5 + {scale*R_j} ->[k] mi_ue5 + {scale*R_5j}
    dm_u_e_dt, dm_5_dt_2, dR_j_dt_2, dmi_ue5_dt_1 , dR_5j_dt_2  = mRNA_binding(m_u_e, m_5, scale_calc(m_u_e, m_b_e, R_j), kbind*TS, mi_ue5)
    
    ## mi_uf5 + X ->[k] X_uf5
    dmi_uf5_dt_2, dX_dt_1 , dR_5j_dt_3, dX_uf5_dt_1 , dR_X5j_dt_1  = mRNA_binding(mi_uf5, X, 0, kX5, 0)   
    ## mi_ue5 + X + {scale*R_5j} ->[k] X_ue5 + {scale*R_X5j}
    dmi_ue5_dt_2, dX_dt_2 , dR_5j_dt_4, dX_ue5_dt_1 , dR_X5j_dt_2  = mRNA_binding(mi_ue5, X, scale_calc(mi_ue5, 0, R_5j), kX5, 0)
    
    ## R_5j -> R + P
    _, _, _, _, dR_dt_1, dR_5j_dt_5, dP_dt = transcription_or_translation(0, mi_uf5, 0, mi_ue5, 0, R_5j, 0, 0, 0, 0, kc)
    
    ## X_uf5 ->[k] X
    dX_uf5_dt_2, dR_X5j_dt_3, dX_dt_3, dR_dt_2 =  mRNA_degradation(X_uf5, 0, kdeg, 0)
    ## X_ue5 + {scale*R_X5j}->[k] X + scale*Ribosomes
    dX_ue5_dt_2, dR_X5j_dt_4, dX_dt_4, dR_dt_3 =  mRNA_degradation(X_ue5, scale_calc(X_ue5, 0, R_X5j), kdeg, 0)
    
    return np.array([dm_u_f_dt, dm_u_e_dt, 0, dmi_uf5_dt_1+dmi_uf5_dt_2, dmi_ue5_dt_1+dmi_ue5_dt_2, dX_uf5_dt_1+dX_uf5_dt_2, dX_ue5_dt_1+dX_ue5_dt_2, dm_5_dt_1+dm_5_dt_2, dX_dt_1+dX_dt_2+dX_dt_3+dX_dt_4, dR_dt_1+dR_dt_2+dR_dt_3, dR_j_dt_1+dR_j_dt_2, dR_5j_dt_1+dR_5j_dt_2+dR_5j_dt_3+dR_5j_dt_4+dR_5j_dt_5, dR_X5j_dt_1+dR_X5j_dt_2+dR_X5j_dt_3+dR_X5j_dt_4, dP_dt])


######### Full mRNA production and normaldegradation pathway
def production_full_pathway(D_u_f, D_b_f, D_u_e, D_b_e, N_j, m_u_f, m_b_f, m_u_e, m_b_e, R_j, P, X_ufj, X_bfj, X_uej, X_bej, R_Xj, N, R, X, kd_p, kd_m, kd, kd_c, km_p, km_m, km, km_c, kX_m, kdeg, lam_P):    
    # Transcription
#     dD_u_f_dt, dD_b_f_dt, dD_u_e_dt, dD_b_e_dt, dN_dt, dN_j_dt, dm_u_f_dt_1 = transcription_or_translation(D_u_f, D_b_f, D_u_e, D_b_e, N, N_j, m_u_f, kd_p, kd_m, kd, kd_c)
    dD_u_f_dt, dD_b_f_dt, dD_u_e_dt, dD_b_e_dt, dN_dt, dN_j_dt, dm_u_f_dt_1 = 0, 0, 0, 0, 0, 0, kd_c
    # Translation
    dm_u_f_dt_2, dm_b_f_dt_1, dm_u_e_dt_1, dm_b_e_dt_1, dR_dt_1, dR_j_dt_1, dP_dt_1 = transcription_or_translation(m_u_f, m_b_f, m_u_e, m_b_e, R, R_j, P, km_p, km_m, km, km_c)
    # Unbound Free mRNA normaldegradation
    dm_u_f_dt_3 , dX_dt_1, dX_ufj_dt, dR_j_dt_2,dR_Xj_dt_1, dR_dt_2 = mRNA_selfdegradation_fullpathway_SIMP(m_u_f, X, X_ufj, 0, 0, kX_m, kdeg, 0)
    # Bound Free mRNA normaldegradation
    dm_b_f_dt_2 , dX_dt_2, dX_bfj_dt, dR_j_dt_3,dR_Xj_dt_2, dR_dt_3 = mRNA_selfdegradation_fullpathway_SIMP(m_b_f, X, X_bfj, 0, 0, kX_m, kdeg, 1)
    # Unbound Elongating mRNA normaldegradation
    dm_u_e_dt_2 , dX_dt_3, dX_uej_dt, dR_j_dt_4,dR_Xj_dt_3, dR_dt_4 = mRNA_selfdegradation_fullpathway_SIMP(m_u_e, X, X_uej, scale_calc(m_u_e, m_b_e, R_j), scale_calc(X_uej, X_bej, R_Xj), kX_m, kdeg, 0)
    # Bound Elongating mRNA normaldegradation
    dm_b_e_dt_2 , dX_dt_4, dX_bej_dt, dR_j_dt_5,dR_Xj_dt_4, dR_dt_5 = mRNA_selfdegradation_fullpathway_SIMP(m_b_e, X, X_bej, scale_calc(m_b_e, m_u_e, R_j), scale_calc(X_bej, X_uej, R_Xj), kX_m, kdeg, 1)
    
    # Protein Dilution
    dP_dt_2= -lam_P*P
    return np.array([dD_u_f_dt, dD_b_f_dt, dD_u_e_dt, dD_b_e_dt, dN_j_dt, dm_u_f_dt_1+dm_u_f_dt_2+dm_u_f_dt_3, dm_b_f_dt_1+dm_b_f_dt_2, dm_u_e_dt_1+dm_u_e_dt_2, dm_b_e_dt_1+dm_b_e_dt_2, dR_j_dt_1+dR_j_dt_2+dR_j_dt_3+dR_j_dt_4+dR_j_dt_5, dP_dt_1+dP_dt_2, dX_ufj_dt, dX_bfj_dt, dX_uej_dt, dX_bej_dt, dR_Xj_dt_1+dR_Xj_dt_2+dR_Xj_dt_3+dR_Xj_dt_4, dN_dt, dR_dt_1+dR_dt_2+dR_dt_3+dR_dt_4+dR_dt_5, dX_dt_1+dX_dt_2+dX_dt_3+dX_dt_4 ])

def production_mRNA_pathway_1TS(D_u_f, D_b_f, D_u_e, D_b_e, N_j, m_u_f, N, R, X, kd_p, kd_m, kd, kd_c, lam_mi):
    # Transcription
    
#     dD_u_f_dt, dD_b_f_dt, dD_u_e_dt, dD_b_e_dt, dN_dt, dN_j_dt, dm_u_f_dt_1 = transcription_or_translation(D_u_f, D_b_f, D_u_e, D_b_e, N, N_j, m_u_f, kd_p, kd_m, kd, kd_c)
    dD_u_f_dt, dD_b_f_dt, dD_u_e_dt, dD_b_e_dt, dN_dt, dN_j_dt, dm_u_f_dt_1 = 0, 0, 0, 0, 0, 0, kd
    # Dilution
    dm_u_f_dt_2 = -m_u_f*lam_mi
    return np.array([dD_u_f_dt, dD_b_f_dt, dD_u_e_dt, dD_b_e_dt, dN_j_dt, (dm_u_f_dt_1)+dm_u_f_dt_2, dN_dt])

def production_mRNA_pathway_3TS(D_u_f, D_b_f, D_u_e, D_b_e, N_j, m_u_f, N, R, X, kd_p, kd_m, kd, kd_c, lam_mi):
    # Transcription
    
#     dD_u_f_dt, dD_b_f_dt, dD_u_e_dt, dD_b_e_dt, dN_dt, dN_j_dt, dm_u_f_dt_1 = transcription_or_translation(D_u_f, D_b_f, D_u_e, D_b_e, N, N_j, m_u_f, kd_p, kd_m, kd, kd_c)
    dD_u_f_dt, dD_b_f_dt, dD_u_e_dt, dD_b_e_dt, dN_dt, dN_j_dt, dm_u_f_dt_1 = 0, 0, 0, 0, 0, 0, kd
    # Dilution
    dm_u_f_dt_2 = -m_u_f*lam_mi
    return np.array([dD_u_f_dt, dD_b_f_dt, dD_u_e_dt, dD_b_e_dt, dN_j_dt, (dm_u_f_dt_1)+dm_u_f_dt_2, dN_dt])

def mRNA_selfdegradation_fullpathway_SIMP(target, X, Y, scale1, scale2, k1, k2, m_or_C):
    # target + X + {scale1*elongating} ->[k1] Y + scale1*blocked
    dtarget_dt, dX_dt_1, delongating_dt, dY_dt_1, dblocked_dt_1 = mRNA_binding(target, X, scale1, k1, 0)
    # Y + {scale2*blocked} ->[k2] X + scale2*freeRib
    dY_dt_2, dblocked_dt_2, dX_dt_2, dR_dt = mRNA_degradation(Y, scale2, k2, m_or_C)
    return np.array([dtarget_dt, dX_dt_1+dX_dt_2, dY_dt_1+dY_dt_2, delongating_dt, dblocked_dt_1+dblocked_dt_2, dR_dt])

######### General mRNA binding function
## target + binder + {scale*elongating_Ribosomes} ->[k] product + {scale*blocked_Ribosomes}
## target + binder + {scale*elongating_Ribosomes} <=>[k] product + {scale*blocked_Ribosomes}
def mRNA_binding(target, binder, scale, k, product):
    K = 0
    dtarget_dt =                         -target*binder*k  + product*K
    dbinder_dt =                         -target*binder*k  + product*K
    delongating_Ribosomes_dt =     -scale*target*binder*k  + scale*product*K
    dproduct_dt =                                                target*binder*k  - product*K
    dblocked_Ribosomes_dt =                                scale*target*binder*k  - scale*product*K
    return np.array([dtarget_dt, dbinder_dt, delongating_Ribosomes_dt, dproduct_dt, dblocked_Ribosomes_dt])
## target + binder + {scale*elongating_Ribosomes} ->[k] product
def mRNA_binding_SIMP(target, binder, scale, k):
    dtarget_dt =                         -target*binder*k
    dbinder_dt =                         -target*binder*k
    delongating_Ribosomes_dt =     -scale*target*binder*k
    dproduct_dt =                                                target*binder*k
    return np.array([dtarget_dt, dbinder_dt, delongating_Ribosomes_dt, dproduct_dt])
## target + binder ->[k] product
def mRNA_binding_SIMP2(target, binder, scale, k):
    dtarget_dt =                         -target*binder*k
    dbinder_dt =                         -target*binder*k
    dproduct_dt =                                                target*binder*k
    return np.array([dtarget_dt, dbinder_dt, dproduct_dt])


######### General mRNA degradation
## Y + {scale*blocked_Ribosomes}->[k] X + scale*Ribosomes
def mRNA_degradation(Y, scale, k, m_or_C):
    dY_dt =                                -Y*k
    dblocked_Ribosomes_dt =          -scale*Y*k
    dX_dt =                                 Y*k
    dRibosomes_dt =          (scale+m_or_C)*Y*k
    return np.array([dY_dt, dblocked_Ribosomes_dt, dX_dt, dRibosomes_dt])

# target -> X + scale*Ribosomes
def mRNA_degradation_SIMP(Y, scale, k, m_or_C):
    dY_dt =                                -Y*k
    dX_dt =                                 Y*k
    dRibosomes_dt =          (scale+m_or_C)*Y*k
    return np.array([dY_dt, dX_dt, dRibosomes_dt])