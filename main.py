# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 17:13:24 2019

@author: 文明研究员
"""
import numpy as np 
import matplotlib.pyplot as plt
from calculator import Methods,Calculator
from const import Const
from scipy.optimize import leastsq

def func(p,x):
    k,b=p
    return k*x+b 
def error(p,x,y):
    return func(p,x)-y/1e+13 

if __name__=='__main__':
    cal=Calculator()
    cal.me.change_const(Const(0.3))
    M_p=np.linspace(5.15,12,200)
    cal.creat_data(M_p,'Data2.xlsx','test.png')