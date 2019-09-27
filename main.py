# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 17:13:24 2019

@author: 文明研究员
"""
from calculator import Methods,Calculator
import numpy as np
import matplotlib.pyplot as plt

if __name__=='__main__':
    #arange=np.linspace(5.07,12,200)
    #cal=Calculator()
    #cal.creat_data(arange)
    cal=Calculator()
    cal.read_excel('Data.xlsx')
    m1,ST1,L_s1=cal.M_p[0],cal.ST[0],cal.L_init[0]
    cal.me.con.set_M_p(m1)
    g=1e-8
    ST_test=np.linspace(0,ST1,100)
    L_in=np.ones(100)
    g_in=np.ones(100)
    index=0
    for i in iter(ST_test):
        L_in[index],g_in[index]=cal.me.cal_ML_simple_B(i,L_s1,True,g)[1:3]
        index+=1
    fig=plt.figure(dpi=100)
    ax1=fig.add_subplot(211)
    ax1.plot(ST_test,L_in,'r--')
    ax2=fig.add_subplot(212)
    ax2.plot(ST_test,g_in,'g--')
    plt.savefig('test.png')