# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 17:13:24 2019

@author: 文明研究员
"""
from calculator import Methods,Calculator
from test_total import Test
import numpy as np
import matplotlib.pyplot as plt

if __name__=='__main__':
    cal=Calculator()
    cal.read_excel('Data.xlsx')
    order=1
    m=cal.M_p[order]
    ST_max=cal.ST[order]
    L_s=cal.L_init[order]
    num=200
    test=Test()
    test.con.set_M_p(m)
    test.cal_ML_simple_B(ST_max*9/10,L_s,True,1e-9,2e-9) 
    test.creat_test_data()
    