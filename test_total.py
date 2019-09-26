# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 20:21:57 2019

@author: 文明研究员
"""
import numpy as np
import matplotlib.pyplot as plt
from calculator import Methods
from scipy.integrate import odeint
import warnings
import sys

class Test(object):
    me=Methods()
    con=me.con
    @classmethod
    def test_range(cls,M_p,L_init,ST_range,plot=False):
        """Use to test the range of ST.
        Can ues plot to decide whether to plot.
        Recommend to set a corrent L."""
        cls.me.con.set_M_p(M_p)
        M_c=cls.me.con.M_c
        L=L_init
        M=cls.me.cal_ML_simple(0,L)[0]
        while abs(M-M_c)>0.001:
            L_next=L+0.01
            point1=[L,M]
            point2=[L_next,cls.me.cal_ML_simple(0,L_next)[0]]
            k=cls.me.grad(point1,point2)
            if k==0:
                raise ZeroDivisionError('k==0')
            else:
                L=(M_c-M)/k+L
                M=cls.me.cal_ML_simple(0,L)[0]
        L_result=[0]*len(ST_range)
        for i in range(len(ST_range)):
            try:
                warnings.simplefilter('error')
                L_result[i]=cls.me.cal_ML_simple(ST_range[i],L)[1]
            except:
                L_result[i]=-10
                continue
        if plot==True:
            fig=plt.figure(dpi=100)
            plt.plot(ST_range,L_result)
            plt.show()
        return L_result
    
    @classmethod
    def ode_fun_test(cls,Z,r,c_P,c_T,c_M,others,ST,J):
        P,T,M,L=Z
        g,alpha,beta=others
        if M<4.08:
            print("M={0}".format(M))
            pass

        dP=c_P*M*P/(T*r**2)
        dT=min(c_T*1e24*L*P**(alpha+1)*T**(beta-4)/M,g)*(T*dP/P)

        dM=c_M*r**2*P/T
        
        if T<0 or P<0:
            print("T,P,L={0},{1},{2}".format(T,P,L))
            print("dT,dP={0},{1}".format(dT,dP))
            breakpoint()
            raise ValueError('Error')
        sigma=cls.con.sigma_1*T**(3./4.)*P**(-1./2.)*np.exp(-cls.con.sigma_2/T)
        if dT==g*(T*dP/P):
            dL=1e-24*cls.me.con.M_e*cls.me.con.T_0*dM*ST
        else:
            dL=0

        dL+=cls.me.con.R_B**3*r**2*J**2/sigma*1e-24
        return [dP,dT,dM,dL]
    
    @classmethod
    def cal_ML_simple(cls,ST,L_s,J):
        c_P,c_T,c_M=cls.con.cal_const()
        g,alpha,beta=cls.con.g_ad,cls.con.alpha,cls.con.beta
        R_in,R_out=cls.con.R_p,cls.con.R_out
        num=800

        r=np.linspace(1,R_in/R_out,num)
        initial=(1.,1.,cls.con.M_p,L_s)
        others=(g,alpha,beta)
        result=odeint(cls.ode_fun_test,initial,r,args=(c_P,c_T,c_M,others,ST,J))
        P,T,M,L=result[:,0],result[:,1],result[:,2],result[:,3]
        cls.P,cls.T,cls.M,cls.L=P,T,M,L
        cls.r=r
        sigma=cls.con.sigma_1*T**(3./4.)*P**(-1./2.)*np.exp(-cls.con.sigma_2/T)
        cls.sigma=sigma

        return M[-1],L[-1]