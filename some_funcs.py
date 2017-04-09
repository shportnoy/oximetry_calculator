# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 19:16:51 2017

@author: sportnoy
"""
import numpy as np
np.seterr(all=None)
from ipywidgets import *
import pylab as plt
import matplotlib

#Read in parameter text file ('par_values.txt')
f = open('par_values.txt')
par_dict = {}
for line in f:
    line = line.strip()
    if not line: continue

    if line.startswith('Parameters'):
        if line not in par_dict: par_dict[line] = []
        heading = line
        continue
    par_dict[heading].append(line)

f.close()

def display_widgets():
    
    from IPython.display import display
    blood_type=widgets.ToggleButtons(
        options=['Adult', 'Fetal'],
        description='Blood type:',
        disabled=False,
        button_style='', 
        tooltip='Select blood type - adult or fetal',
        )

    field_strength=widgets.ToggleButtons(
        options=['1.5T', '3T'],
        description='Field Strength:',
        disabled=False,
        button_style='', 
        tooltip='Select field strength - 1.5T or 3T',
    ) 
    
    tau_180=widgets.IntSlider(
        value=16,
        min=1,
        max=50,
        step=1,
        description='Refocusing interval (ms)',
        disabled=False,
        continuous_update=False,
        orientation='horizontal',
        readout=True,
        readout_format='i',
        slider_color='white'
        )
        
    display(Box([field_strength, blood_type, tau_180]))  
   
    return(field_strength, blood_type, tau_180)


def display_sO2_T1_boxes():
    
    from IPython.display import display
    T1=widgets.BoundedIntText(
        value=1500,
        min=0,
        max=5000,
        step=1,
        description='T1 (ms):',
        disabled=False,
        color='black'
    )
    
    sO2=widgets.BoundedFloatText(
        value=0.85,
        min=0,
        max=1,
        step=0.01,
        description='sO2:',
        disabled=False,
        color='black'
    )

    display(Box([sO2, T1]))  
    
    return(sO2, T1)
    
def display_Hct_T2_boxes():
    
    from IPython.display import display
    T2=widgets.BoundedFloatText(
        value=150,
        min=0,
        max=500,
        step=1,
        description='T2 (ms):',
        disabled=False,
        color='black'
    )
    
    Hct=widgets.BoundedFloatText(
        value=0.45,
        min=0,
        max=1,
        step=0.01,
        description='Hct:',
        disabled=False,
        color='black'
    )

    display(Box([Hct, T2]))  
    
    return(Hct, T2)    
    
def display_T1_T2_boxes():  
    
    from IPython.display import display
    T2=widgets.BoundedFloatText(
        value=150,
        min=0,
        max=500,
        step=1,
        description='T2 (ms):',
        disabled=False,
        color='black'
    )
    
    T1=widgets.BoundedIntText(
        value=1500,
        min=0,
        max=5000,
        step=1,
        description='T1 (ms):',
        disabled=False,
        color='black'
    )
    
    display(Box([T1, T2]))  
    
    return(T1, T2)  


def fetch_model_pars(param_file_dict, blood_type='Fetal', field_strength='3T'):
    """After reading 'par_values.txt' file into dictionary 
      ('param_file_dict'), make dictionary of two-compartment model 
       parameters appropriate for blood type (Adult or Fetal) and field
       strength (1.5T or 3T)"""
    
    if (blood_type not in ['Fetal', 'Adult']):
        print "input parameter blood_type must be 'Adult' or 'Fetal'"
        return
        
    if (field_strength not in ['1.5T', '3T']):
        print "input parameter field_strength must be '1.5T' or '3T'"
        return    
     
    
    key_str = [k for k in param_file_dict.keys() if blood_type in k and field_strength in k][0]
    string_pars=param_file_dict[key_str]   
    par_dict={}
    for item in string_pars:
        key=item.split('=')[0]
        value=float(item.split('=')[1])
        par_dict[key]=value
    
    if '3T' in field_strength:    
        par_dict['field']=3.0
    if '1.5T' in field_strength:
        par_dict['field']=1.5

    output_str='Blood type = %s'%blood_type 
    print output_str
    
    output_str='Field strength = %s'%field_strength
    print output_str
   
    
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
    w_dia_oxy=par_dict['w_dia_plus_w_oxy']*42.576*2*np.pi*par_dict['field']  #convert from ppm to rad/sec
    w_de_oxy=par_dict['wdeo_minus_w_oxy']*42.576*2*np.pi*par_dict['field']
    
    mu=tau*(1-(2*tau/tau_180)*np.tanh(tau_180/(2*tau)))
   
    output_str='Refocusing interval = %.0f ms'%(1000*tau_180)
    print output_str
    
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
    if not isinstance(Hct,np.ndarray):  
        output_str='solution 1: sO2=%.2f\nsolution 2: sO2=%.2f'%(sO2_alt, sO2)
        print output_str
    sols=[sO2, sO2_alt]
    return sols    

def Hct_from_sO2_T1(coeff_list, sO2=0.85, T1=1.3):
     """Calculate hematocrit from T1 relaxation time (in seconds),
       given a specific sO2 value"""
     K0,K1,K2,K3,K4,K5,M0,M1,M2=[i for i in coeff_list]  
     R1=1.0/T1
     Hct=(R1-M0)/(M1+M2*sO2)
     if not isinstance(sO2,np.ndarray):
         output_str='Hct=%.2f'%Hct
         print output_str
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
    
    for jj in np.arange(3):
        if ~np.isreal(sO2_sol[jj] or Hct_sol[jj]):
            output_str='solution {0}: sO2 = {1:.2f} {2} {3:.2f}i, Hct = {4:.2f} {5} {6:.2f}i'.format(jj+1, sO2_sol[jj].real, '+-'[sO2_sol[jj].imag < 0], abs(sO2_sol[jj].imag), 
                                                                                                           Hct_sol[jj].real, '+-'[Hct_sol[jj].imag < 0], abs(Hct_sol[jj].imag))
        else:
            output_str='solution {0}: sO2 = {1:.2f}, Hct = {2:.2f}'.format(jj+1, sO2_sol[jj].real, Hct_sol[jj].real)                                                                                        
            
        print output_str
        
    return Hct_sol, sO2_sol    

def calc_R2(x, R2plas, R2ery0, R2deo, w0, wdeo, tau):
    mu = tau*(1-(2*tau)/x[2]*np.tanh(x[2]/(2*tau))) 
    R2ery = R2ery0 + (1-x[1])*R2deo
    R20=x[0]*R2ery + (1-x[0])*R2plas
    delta=w0+(1-x[1])*wdeo
    R2ex=x[0]*(1-x[0])*delta**2*mu
    R2=R20+R2ex
    return R2 
    
def calc_R1(x, R1plas, R1ery0, r1dHb):   
    R1=x[0]*(R1ery0+r1dHb*(1-x[1]))+(1-x[0])*R1plas
    return R1
    
def plot_T2_sO2(par_dict, sO2_sols, tau_180=0.016, Hct=0.45, T2=0.15):
    tau=par_dict['tau']
    R2_de_oxy=par_dict['R2_deoxy_minus_R2_oxy']    
    R2_dia_oxy=par_dict['R2_dia_plus_R2_oxy']
    R2plas=par_dict['R2plas']
    w_dia_oxy=par_dict['w_dia_plus_w_oxy']*42.576*2*np.pi*par_dict['field']  #convert from ppm to rad/sec
    w_de_oxy=par_dict['wdeo_minus_w_oxy']*42.576*2*np.pi*par_dict['field']
    
    R2ery0=R2_dia_oxy+R2plas
    w0=w_dia_oxy
    
    sO2_axis=np.arange(0,1.001,0.001)
    xtemp=np.zeros([3,len(sO2_axis)])
    xtemp[0,:]=Hct
    xtemp[1,:]=sO2_axis
    xtemp[2,:]=tau_180
    
    T2s=1000/calc_R2(xtemp, R2plas, R2ery0, R2_de_oxy, w0, w_de_oxy, tau)
    
    plt.xlabel('oxygen-saturation, $sO_2$', fontsize=16)
    plt.ylabel('$T_2$ relaxation time (ms)', fontsize=16)
    
    plt.plot(sO2_axis,T2s,lw=2)
    plt.ylim([0.95*np.min(T2s),1.05*np.max(T2s)])
    plt.xlim([0,1])

    plt.hlines(1000*T2,0,1,linestyle='dashed')
    plt.vlines(sO2_sols[0],0.95*np.min(T2s),1.05*np.max(T2s),linestyle='dashed')
    plt.vlines(sO2_sols[1],0.95*np.min(T2s),1.05*np.max(T2s),linestyle='dashed')
    
    return
    

    
def plot_T1_Hct(par_dict, Hct_sol, sO2=0.5, T1=1.5):
    
    R1ery0=par_dict['R1ery0']
    R1plas=par_dict['R1plas']
    r1dHb=par_dict['r1_prime_dHb']
    
    Hct_axis=np.arange(0,1.001,0.001)
    xtemp=np.zeros([2,len(Hct_axis)])
    xtemp[0,:]=Hct_axis
    xtemp[1,:]=sO2
    
    T1s=1000/calc_R1(xtemp, R1plas, R1ery0, r1dHb)
    
    plt.xlabel('hematocrit (Hct)', fontsize=16)
    plt.ylabel('$T_1$ relaxation time (ms)', fontsize=16)
    
    plt.plot(Hct_axis,T1s,lw=2)
    
    plt.ylim([0.95*np.min(T1s),1.05*np.max(T1s)])
    plt.xlim([0,1])
    
    plt.hlines(1000*T1,0,1,linestyle='dashed')
    plt.vlines(Hct_sol,0.95*np.min(T1s),1.05*np.max(T1s),linestyle='dashed')
    
    return
    
def plot_T1_T2_Hct_sO2(coeff_list, Hct_sol, sO2_sol, T1=1.5, T2=0.15):   
    
    Hct_axis=np.arange(0,1.001,0.001)
    sO2_axis=np.arange(0,1.001,0.001)
    
    T1_curve=Hct_from_sO2_T1(coeff_list,sO2=sO2_axis,T1=T1)
    T2_curves=sO2_from_Hct_T2(coeff_list,Hct=Hct_axis,T2=T2)
    
    plt.plot(T1_curve,sO2_axis,color='r',lw=2, label='$T_1=%.0f$'%(1000*T1))
    plt.plot(Hct_axis,T2_curves[0],color='b',lw=2,label='$T_2=%.1f$'%(1000*T2))
    plt.plot(Hct_axis,T2_curves[1],color='b',lw=2)
    plt.xlim([0,1])
    plt.ylim([0,1])
    
    plt.xlabel('hematocrit (Hct)', fontsize=16)
    plt.ylabel('oxygen-saturation, $sO_2$', fontsize=16)
    plt.legend(loc='upper right', bbox_to_anchor=(1.45,1.04), fontsize=16)
    
    sols=zip(Hct_sol, sO2_sol)
    
    for sol in sols:
        if np.any(np.array(sol)>1) or np.any(np.array(sol)<0) or np.any(np.abs(np.imag(np.array(sol)))>1e-5):
            sol=[np.nan, np.nan]
        plt.plot(sol[0], sol[1], 'o', ms=12, mfc='k', alpha=0.5)
    
    return
    
