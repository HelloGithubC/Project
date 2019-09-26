import numpy as np 
import matplotlib.pyplot as plt 
from calculator import Methods

if __name__=='__main__':
    me=Methods()
    me.con.set_M_p(6.00150753768844)
    L_s=1.04275883779695
    ST_test=np.linspace(1e-7,3e-6,100)
    L_c=np.ones(100)
    index=0
    for i in iter(ST_test):
        L_c[index]=me.cal_ML_simple(i,L_s)[1]
        index+=1
    
    plt.figure(dpi=100)
    plt.plot(ST_test,L_c)
    plt.show()
