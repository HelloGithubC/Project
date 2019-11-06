import numpy as np 
import matplotlib.pyplot as plt 
from numba import jit
from calculator import Methods,Calculator
from methods.NA import grad
from scipy import optimize

def creat_plot_2D(x,y,ref,plot_style='r--',plot_style_ref='g--',save=False,scatter=False):
    plt.figure(dpi=100)
    y=abs(y)
    plt.loglog(x,y,plot_style)
    plt.loglog(x,np.ones(y.shape)*ref,plot_style_ref)
    if scatter:
        plt.scatter(scatter[0],scatter[1],marker='o')
    if save:
        if type(save) is str:
            plt.savefig(save)
        else:
            plt.savefig('result.png')
    else:
        plt.show()

def solve_2D(func,x_init,y_target,x_step,error):
    x=x_init
    y=func(ST_max/2,L_s,True,1e-9,x)[2] #Must modify according the func you need
    while abs(y-y_target)>error:
        x_next=x+x_step
        y_next=func(ST_max/2,L_s,True,1e-9,x)[2]
        k=grad([x,y],[x_next,y_next])
        if k==0:
            raise ZeroDivisionError('k=0')
        x=(y_target-y)/k+x 
        y=func(ST_max/2,L_s,True,1e-9,x)[2]
    return (x,y)


if __name__ == "__main__":
    cal=Calculator()
    cal.read_excel('Data.xlsx')
    order=1
    m=cal.M_p[order]
    cal.me.con.set_M_p(m)
    ST_max=cal.ST[order]
    L_s=cal.L_init[order]
    num=200
    




