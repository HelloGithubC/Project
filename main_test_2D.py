import numpy as np 
import matplotlib.pyplot as plt 
from numba import jit
from calculator import Methods,Calculator

def creat_plot_2D(x,y,ref,plot_style='r--',plot_style_ref='g--',save=False):
    plt.figure(dpi=100)
    plt.plot(x,y,plot_style)
    plt.plot(x,np.ones(y.shape)*ref,plot_style_ref)
    if save:
        if type(save) is str:
            plt.savefig(save)
        else:
            plt.savefig('result.png')
    else:
        plt.show()

if __name__ == "__main__":
    cal=Calculator()
    cal.read_excel('Data.xlsx')
    order=1
    m=cal.M_p[order]
    cal.me.con.set_M_p(m)
    ST_max=cal.ST[order]
    L_s=cal.L_init[order]
    ST_test=np.linspace(0,ST_max,100)
    g_out_test=np.linspace(0,1e-9,100)
    M_c=np.ones(ST_test.shape)
    for i in range(ST_test.shape[0]):
        M_c[i]=cal.me.cal_ML_simple_B(ST_test[i],L_s,True,1e-11)[0]
    creat_plot_2D(ST_test,M_c,5.0)
    for i in range(g_out_test.shape[0]):
        M_c[i]=cal.me.cal_ML_simple_B(ST_max/2,L_s,True,g_out_test[i])[0]
    creat_plot_2D(g_out_test,M_c,5.0)