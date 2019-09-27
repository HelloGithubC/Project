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
    g_test=np.linspace(1e-7,1e-4,100)
    g_in=np.ones(100)
    L_in=np.ones(100)
    index=0
    for i in iter(g_test):
        L_in[index],g_in[index]=cal.me.cal_ML_simple_B(ST1,L_s1,True,i)[1:3]
        index+=1
    fig=plt.figure(dpi=100)
    ax1=fig.add_subplot(211)
    ax1.semilogx(g_test,g_in,'r--')
    ax2=fig.add_subplot(212)
    ax2.semilogx(g_test,L_in,'g--')
    plt.savefig('test.png')