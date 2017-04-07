# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 19:16:51 2017

@author: sportnoy
"""
import numpy as np

def fetch_model_pars(param_file_dict, blood_type='Fetal', field_strength='3T'):
    """After reading 'par_values.txt' file into dictionary 
      ('param_file_dict'), make dictionary of two-compartment model 
       parameters appropriate for blood type (Adult or Fetal) and field
       strength (1.5T or 3T)"""
       
    key_str = [k for k in param_file_dict.keys() if blood_type and field_strength in k][0]
    string_pars=param_file_dict[key_str]   
    par_dict={}
    for item in string_pars:
        key=item.split('=')[0]
        value=float(item.split('=')[1])
        par_dict[key]=value   
    return par_dict
    
def calc_coeffs(par_dict, tau_180=0.016):
    """Calculate list of intermediate coeffients, based on model parameters 
       (in par_dict) and refocusing interval, tau_180 (in seconds)"""   
    
    tau=par_dict['tau']
    R1ery0=par_dict['R1ery0']
    R1plas=par_dict['R1plas']
    R2_de_oxy=par_dict['R2_deoxy_minus_R2_oxy']    
    R2_dia_oxy=par_dict['R2_dia_plus_R2_oxy']
    R2plas=par_dict['R2plas']
    r1dHb=par_dict['r1_prime_dHb']
    w_dia_oxy=par_dict['w_dia_plus_w_oxy']*42.576*2*np.pi  #convert from ppm to rad/sec
    w_de_oxy=par_dict['wdeo_minus_w_oxy']*42.576*2*np.pi
    
    mu=tau*(1-(2*tau/tau_180)*np.tanh(tau_180/(2*tau)))
    
    K0 = R2plas
    K1 = R2_dia_oxy + R2_de_oxy + mu*(w_dia_oxy + w_de_oxy)**2
    K2 = -mu*(w_dia_oxy + w_de_oxy)**2
    K3 = -R2_de_oxy - 2*mu*w_de_oxy*(w_dia_oxy + w_de_oxy)
    K4 = 2*mu*w_de_oxy*(w_dia_oxy + w_de_oxy)
    K5 = mu*w_de_oxy**2
    M0 = R1plas
    M1 = R1ery0 - R1plas + r1dHb
    M2 = -r1dHb  
    
    coeff_list=[K0,K1,K2,K3,K4,K5,M0,M1,M2]
    
    return coeff_list
    
def sO2_from_Hct_T2(coeff_list, Hct=0.45, T2=0.15):
    """Calculate oxygen saturation from T2 relaxation time (in seconds),
       given a specific Hct value"""
    K0,K1,K2,K3,K4,K5,M0,M1,M2=[i for i in coeff_list]
    R2=1.0/T2
    root_term=(K4*Hct**2+K3*Hct)**2-4*((K5*Hct-K5*Hct**2)*(K0+K1*Hct+K2*Hct**2-R2))
    sO2=((-K4*Hct**2-K3*Hct)+np.sqrt(root_term))/(2*(K5*Hct-K5*Hct**2))
    sO2_alt=((-K4*Hct**2-K3*Hct)-np.sqrt(root_term))/(2*(K5*Hct-K5*Hct**2))
    sols=[sO2, sO2_alt]
    return sols    

def Hct_from_sO2_T1(coeff_list, sO2=0.85, T1=1.3):
     """Calculate hematocrit from T1 relaxation time (in seconds),
       given a specific sO2 value"""
     K0,K1,K2,K3,K4,K5,M0,M1,M2=[i for i in coeff_list]  
     R1=1.0/T1
     Hct=(R1-M0)/(M1+M2*sO2)
     return Hct    
    
def Hct_sO2_from_T1_T2(coeff_list, T1=1.5, T2=0.1): 
    """Calculate Hct and sO2 given T1 and T2 relaxation times (in seconds)"""
    K0,K1,K2,K3,K4,K5,M0,M1,M2=[i for i in coeff_list]  
    R1=1.0/T1
    R2=1.0/T2
    A1 = K2*M2**2 - K4*M1*M2 - K5*M1**2
    A2 = K5*(R1*M2 - M0*M2)
    B1 = K1*M2**2 - K3*M1*M2 + K4*(R1*M2 - M0*M2) + K5*(M1**2 + 2*R1*M1 - 2*M0*M1)
    B2 = K5*(R1*M1 - M0*M1 - (M0-R1)**2) + K3*(R1*M2 - M0*M2) + K0*M2**2 - R2*M2**2
    C1 = K0*M2**2 + K3*(R1*M2 - M0*M2) + K5*(2*M0*M1 - 2*M1*R1 - (R1-M0)**2) - R2*M2**2
    C2 = K4*(M0 - R1)**2 + K3*(R1*M1 - M0*M1) + K1*(R1*M2 - M0*M2) + 2*K0*M1*M2 - 2*R2*M1*M2
    D1 = K5*(R1 - M0)**2
    D2 = K2*(M0 - R1)**2 + K1*(R1*M1 - M0*M1) + K0*M1**2 - R2*M1**2
   
    Hct_sol=np.roots([A1,B1,C1,D1])
    sO2_sol=np.roots([A2,B2,C2,D2])
    
    return Hct_sol, sO2_sol    
    

    
