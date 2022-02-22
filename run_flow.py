#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 12:04:03 2020

@author: rosie
"""

import os, sys
from time import time
from pydfnworks import * 
from numpy import *

def rewrite_file(filename,b,three=False):
    with open(filename,'r') as f:
        lines = f.readlines()
    
    j = 0
    header = True
    s = ' '
    for line in lines:
        if header:
            header = False
        else:
            values = line.strip().split()
            if three:
                values[-3] = str(b[j])
                values[-2] = str(b[j])
                
            values[-1] = str(b[j]) + '\n'
    
            j = j + 1
            lines[j] = s.join(values)
    
    with open(filename,'w') as f:
        for line in lines:
            f.writelines(line)
        
        
define_paths()
main_time = time()
DFN = create_dfn()

DFN.set_flow_solver("PFLOTRAN")
DFN.dfn_gen(output=True)

b = zeros(DFN.num_frac)
b[:3]= 1e-3
b[3] = 5e-4 #change the aperture of the fourth fracture

cwd = os.getcwd()
##############
#perm

perm = b**2/12 # calculate permeability as function of aperture

############
#dump into file
filename_aper = cwd + '/aperture.dat'
filename_perm = cwd + '/perm.dat'

rewrite_file(filename_aper,b)
rewrite_file(filename_perm,perm,three=True)

DFN.dfn_flow()
