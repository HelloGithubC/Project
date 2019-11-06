# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 17:13:24 2019

@author: 文明研究员
"""
import numpy as np 
import matplotlib.pyplot as plt
from calculator import Methods,Calculator
from scipy.optimize import leastsq

def func(p,x):
    k,b=p
    return k*x+b 
def error(p,x,y):
    return func(p,x)-y/1e+13 

if __name__=='__main__':
    cal=Calculator()
    cal.read_excel('Data.xlsx')
    order=1
    m=cal.M_p[order]
    cal.me.con.set_M_p(m)
    ST_max=cal.ST[order]
    L_s=cal.L_init[order]
    num=200

    ST=ST_max*9/10
    dG=np.linspace(0,1e-9,num)
    g_in=np.ones(num)
    for i in range(num):
        g_in[i]=cal.me.cal_ML_simple_B(ST,L_s,True,1e-10,dG[i])[2]
    p0=[-1,2.5]
    para=leastsq(error,p0,args=(dG,g_in))
    k,b=para[0] 
    print(k,b)
    print(-b/k)
    