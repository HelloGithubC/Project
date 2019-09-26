'''
@Author: Xiao Liang
@Date: 2019-07-26 11:00:27
@LastEditors: Xiao Liang
@LastEditTime: 2019-09-22 21:08:00
@Description: file content
'''
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 10:59:08 2019

@author: 文明研究员
"""
import numpy as np
from calculator import Methods,ToAdjust
import matplotlib.pyplot as plt
from test import Test
import sys

if __name__=='__main__':
    me=Methods()
    ST,L=2.46786405544956E-06,2.81916528606318
    me.cal_ML_simple(ST,L)
        
    #fig=plt.figure(dpi=100)
    #plt.loglog(me.r,me.L,'g-')
    store1=[]
    for i in range(len(me.r)):
        store1.append([me.r[i],me.P[i],me.T[i],me.M[i],me.L[i],me.g[i],me.dg[i]])

    while True:
        try:
            me.cal_ML_simple(ST,L,True)
        except ToAdjust:
            ST-=ST/10
            print(ST,me.L_test)
        else:
            break
    #me.cal_ML_simple(0.0,L,True)
    print(ST)
    store2=[]
    for i in range(len(me.r)):
        store2.append([me.r[i],me.P[i],me.T[i],me.M[i],me.L[i],me.g[i],me.dg[i]])
    test=Test()
    test.write_excel('Store.xlsx',store1,store2)
   # plt.loglog(me.r,me.L,'r--')
   # plt.savefig('test.png')
    
