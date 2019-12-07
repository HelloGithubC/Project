'''
@Author: Xiao Liang
@Date: 2019-09-08 23:02:20
@LastEditors: Xiao Liang
@LastEditTime: 2019-09-22 20:57:46
@Description: file content
'''

import numpy as np 
import matplotlib.pyplot as plt
from const import Const
from scipy.integrate import odeint

class Test(object):
    def __init__(self):
        self.con=Const()
        self.con.depth=1.0
    
    def ode_fun(self,Z,r,c_P,c_T,c_M,others,ST):
        P,T,M,G,DG,L=Z 
        g,alpha,beta,Lambda,R_B=others 
        dP=c_P*M*P/(T*r**2)
        dT=min(c_T*1e24*abs(L)*P**(alpha+1)*T**(beta-4)/M,g)*(T*dP/P)
        dM=c_M*r**2*P/T
        dG=DG
        ddG=6*G/r**2+((-0.5*dP/P)+(0.75/T+self.con.sigma_2/T**2)*dT)*dG-(r>self.con.depth)*2*self.con.M_v*(1/r-self.con.depth**2/r**3)/self.con.c**2/self.con.R_B
        
        if dT==g*(T*dP/P):
            dL=1e-24*self.con.M_e*self.con.T_0*dM*ST 
        else:
            dL=0
        
        return [dP,dT,dM,dG,ddG,dL]
    
    def cal_ML_simple(self,ST,L_s,G=0.0,dG=1e-11):
        c_P,c_T,c_M=self.con.cal_const()
        g,alpha,beta=self.con.g_ad,self.con.alpha,self.con.beta
        R_in,R_out=self.con.R_p,self.con.R_out
        num=5000

        r=np.linspace(1,R_in/R_out,num)
        initial=(1.,1.,self.con.M_p,G,dG,L_s)
        others=(g,alpha,beta,self.con.Lambda,self.con.R_B)
        result=odeint(self.ode_fun,initial,r,args=(c_P,c_T,c_M,others,ST))
        P,T,M,G,dg,L=result[:,0],result[:,1],result[:,2],result[:,3],result[:,4],result[:,5]
        self.P,self.T,self.M,self.L=P,T,M,L
        self.r,self.g=r,G
        self.dg=dg
        return M[-1],L[-1],G[-1]