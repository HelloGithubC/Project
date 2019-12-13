# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 20:21:57 2019

@author: 文明研究员
"""
import numpy as np
import matplotlib.pyplot as plt
from calculator import Methods
from scipy.integrate import odeint
from openpyxl import Workbook,load_workbook
from valueDataBase import DataBase

class Test(object):
    me=Methods()
    con=me.con
    values=[]
    dValues=[]
    dr=0.0
    num=800
    f=[]
    @classmethod
    def ode_fun_B(cls,Z,r,c_P,c_T,c_M,others,ST):
        P,T,M,G,DG,L=Z
        cls.values.append([r,P,T,M,G,DG,L])
        #if L<-1000:
            #raise LOutOfRangeError('L<-1000')
        g,alpha,beta,Lambda,R_B=others

        dP=c_P*M*P/(T*r**2)
        dT=min(c_T*1e24*abs(L)*P**(alpha+1)*T**(beta-4)/M,g)*(T*dP/P)
        dM=c_M*r**2*P/T
        dG=DG
        dL=Lambda*7.15e-5*(dG**2+G**2/r**2)/(cls.sigma(P,T)*R_B)
        if dT==g*(T*dP/P):
            dL+=1e-24*cls.con.M_e*cls.con.T_0*dM*ST
        else:
            dL+=0
        P_0=cls.con.P_0
        ddG=6*G/r**2+(0.5*1e3*(cls.f(np.log10(P*P_0)-6+1e-3)-cls.f(np.log10(P*P_0)-6-1e-3))*dP/P)*dG+(r>cls.con.depth)*2*cls.con.M_v*(1/r-cls.con.depth**2/r**3)/cls.con.c**2/cls.con.R_B

        return [dP,dT,dM,dG,ddG,dL]

    @classmethod
    def find_RCB_index(cls):
        P,T,M,L=cls.P,cls.T,cls.M,cls.L
        judge=(cls.con.c_T*1e24*abs(L)*P**(cls.con.alpha+1)*T**(cls.con.beta-4)/M)<(np.ones(len(L))*cls.con.g_ad)
        for i in range(len(judge)):
            if judge[i]!=judge[0]:
                cls.RCB_index=i
                return i 
        return 0

    @classmethod
    def cal_ML_simple_B(cls,ST,L_s,B=False,G=1e-12,dG=1e-15):
        c_P,c_T,c_M=cls.con.cal_const()
        g,alpha,beta=cls.con.g_ad,cls.con.alpha,cls.con.beta
        R_in,R_out=cls.con.R_p,cls.con.R_out
        num=cls.num
        base=DataBase()
        base.creat_sigma()
        cls.f=base.f_sigma

        r=np.linspace(1,R_in/R_out,num)
        cls.dr=r[0]-r[1]
        if B:
            initial=(1.,1.,cls.con.M_p,G,dG,L_s)
        else:
            initial=(1.,1.,cls.con.M_p,0.0,0.0,L_s)
        others=(g,alpha,beta,cls.con.Lambda,cls.con.R_B)
        result=odeint(cls.ode_fun_B,initial,r,args=(c_P,c_T,c_M,others,ST))
        P,T,M,G,dg,L=result[:,0],result[:,1],result[:,2],result[:,3],result[:,4],result[:,5]
        cls.P,cls.T,cls.M,cls.L=P,T,M,L
        cls.r,cls.g=r,G
        cls.dg=dg

        cls.find_RCB_index()
        return M[-1],L[-1],G[-1]

    @classmethod
    def sigma(cls,P,T):
        return cls.con.sigma_1*T**(3./4)*P**(-1./2)*np.exp(-cls.con.sigma_2/T)

    @classmethod
    def creat_test_data(cls,fileName='TestData.xlsx'):
        wb=Workbook()
        wb.create_sheet('result',index=0)
        sheet=wb[wb.sheetnames[0]]
        length=len(cls.values)
        sheet.append(['Information: ','r_range=linspace({0},{1},{2})'.format(cls.values[0][0],cls.values[length-1][0], length),'M_c={0}M_e'.format(cls.con.M_e),'M_p={0}M_e'.format(cls.con.M_p),'a={0}AU'.format(cls.con.a)])
        sheet.append(['r','P','T','M','G','DG','L'])
        for i in range(len(cls.values)):
            sheet.append(cls.values[i])
        sheet.append(['dr','dP','dT','dM','dG','dDG','dL'])
        for i in range(len(cls.values)):
            sheet.append(cls.dValues[i])
        wb.save(fileName)