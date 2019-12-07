import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_bvp
from const import Const

class Solve:
    def __init__(self):
        self.con=Const()
    def fun(self,Y,r,c_P,c_T,c_M,others,ST):
        P,T,M,L=Y
        g,alpha,beta=others
        dP=c_P*M*P/(T*r**2)
        dT=min(c_T*1e24*abs(L)*P**(alpha+1)*T**(beta-4)/M,g)*(T*dP/P)

        dM=c_M*r**2*P/T

        if dT==g*(T*dP/P):
            dL=1e-24*self.con.M_e*self.con.T_0*dM*ST

        else:
            dL=0
        return np.vstack([dP,dT,dM,dL])      