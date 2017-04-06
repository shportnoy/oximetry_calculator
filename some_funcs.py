# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 19:16:51 2017

@author: sportnoy
"""

def make_par_dict(param_file_dict, blood_type='Fetal', field_strength='3T'):
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