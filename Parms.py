# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 15:07:01 2017

@author: cathy
"""
from random import Random

class Parms:
    Rnd=Random(12345)
    initN={'Sh':499,'Sl':9500,'I1h':0, 'I1l':0,'I2h':1, 'I2l':0} 
    initN_homo={'S':9999, 'I1':0, 'I2':1}
    N=1e4
    y0=[initN['Sh'], initN['Sl'], initN['I1h'], initN['I2h'], initN['I2h'], initN['I2l']]
    #y0_homo=initN_homo
    m=1/40.		           #Birth and death rates per year in both children and older
    gammas=[4., 0.1]
    
    w=0.0
    fH=0.1
    fL=1-fH
    g1=4.
    g2=0.1
    r=26
    B1=0.05
    B2=B1/r
    mH=0.0
    mL=0.0
    MaxSimTime=10
    C=24.
    rho=8
    Ch=C*rho/(fL+rho*fH)
    Cl=C/(fL+rho*fH)