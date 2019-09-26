# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 17:59:31 2019

@author: 文明研究员
"""
import numpy as np
import matplotlib.pyplot as plt
from calculator import Calculator


if __name__=='__main__':
    cal=Calculator()
    arange=np.linspace(5.07,12,200) #5.07 is the min value could calculate
    cal.cal_M_range(arange)
    cal.creat_t(arange,False)
    plt.figure(dpi=100)
    plt.plot(cal.t,arange,label=r'$M_{p}-t$')
    plt.xlabel('t/Myr')
    plt.ylabel(r'$M_{p}$')
    plt.show()
    cal.write_excel(arange)
    