# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 17:13:24 2019

@author: 文明研究员
"""
from calculator import Methods,Calculator
import numpy as np
import matplotlib.pyplot as plt

if __name__=='__main__':
    arange=np.linspace(5.07,12,500)
    cal=Calculator()
    cal.creat_data(arange)
    fig=plt.figure(dpi=100)
    plt.plot(cal.t,cal.M_p)
    plt.savefig('result.png')
    cal.write_excel(arange,'Data1.xlsx')
    