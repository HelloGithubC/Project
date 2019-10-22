import numpy as np 
import matplotlib.pyplot as plt 
from calculator import Methods

if __name__=='__main__':
    me=Methods()
    me.con.set_M_p(5.04)
    ST=0.0
    L_s_test=np.linspace(1000,20000,1000)
    m=np.ones(1000)
    index=0
    for i in iter(L_s_test):
        m[index]=me.cal_ML_simple(ST,i)[0]
        index+=1
    
    plt.figure(dpi=100)
    plt.plot(L_s_test,m)
    plt.show()
