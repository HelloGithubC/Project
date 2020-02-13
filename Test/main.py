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
from calculator import Methods,Calculator
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

if __name__=='__main__':
    me=Methods()
    #print(me.con.P_0)

    p=np.linspace(1,1e+4,1000)
    t=np.linspace(1,100,1000)
    P,T=np.meshgrid(p,t)
    stepP=1
    stepT=1e-3
    sigma_0=me.sigma(P,T)
    dP=(me.sigma(P+stepP,T)-sigma_0)/stepP/sigma_0
    dT=(me.sigma(P,T+stepT)-sigma_0)/stepT/sigma_0

    fig=plt.figure(dpi=300)
    ax1=fig.add_subplot(221,projection='3d')
    ax1.plot_surface(P,T,dP,cmap='rainbow')
    ax2=fig.add_subplot(222,projection='3d')
    ax2.plot_surface(P,T,dT,cmap='rainbow')
    
    ax3=fig.add_subplot(223)
    ax3.plot(T[:,300],dP[:,300])
    ax4=fig.add_subplot(224)
    ax4.plot(P[300],dT[300])
    plt.show()