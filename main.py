# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 17:13:24 2019

@author: 文明研究员
"""
from calculator import Methods
from test_total import Test
import matplotlib.pyplot as plt

if __name__=='__main__':
    M_p=6.4
    me=Methods()
    me.con.set_M_p(M_p)
    m1,m2=me.find_L_M([5.0,0.0],[0.0,2.0],[1e-4,1e-6],[1e-4,3e-2])
    test=Test()
    test.cal_ML_simple(m1[0],m1[1],0)
    P,T,L,M=test.P,test.T,test.L,test.M
    test.cal_ML_simple(m1[0]*0.9,m1[1],1.8e-3)
    print("L_c={0}".format(test.L[-1]))
    
    plt.figure(dpi=120)
    """plt.subplot(221)
    plt.loglog(test.me.r,test.me.P)
    plt.loglog(test.me.r,P,'r--')
    plt.subplot(222)
    plt.loglog(test.me.r,test.me.T)
    plt.loglog(test.me.r,T,'r--')
    plt.subplot(223)
    plt.loglog(test.me.r,test.me.L)
    plt.loglog(test.me.r,L,'r--')
    plt.subplot(224)
    plt.semilogx(test.me.r,test.me.M)
    plt.semilogx(test.me.r,M,'r--')"""
    plt.loglog(test.P,test.T)
    plt.loglog(P,T,'r--')
    
    fig=plt.figure(dpi=120)
    plt.loglog(test.P,test.L)
    plt.loglog(P,L,'r--')
    plt.show()