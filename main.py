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
    for j in range(7):
        m,ST,L_s=cal.M_p[j],cal.ST[j],cal.L_init[j]
        cal.me.con.set_M_p(m)
        g_test=np.linspace(1e-15,1e-9,100)
        ST/=10
        L_in=np.ones(100)
        g_in=np.ones(100)
        index=0
        for i in iter(g_test):
            L_in[index],g_in[index]=cal.me.cal_ML_simple_B(ST,L_s,True,i)[1:3]
            index+=1
        fig=plt.figure(dpi=100)
        ax1=fig.add_subplot(211)
        ax1.plot(g_test,L_in,'r--')
        ax2=fig.add_subplot(212)
        ax2.plot(g_test,g_in,'g--')
        plt.savefig('test_{0}.png'.format(j+1))
        plt.close()