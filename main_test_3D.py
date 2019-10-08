import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from numba import jit
from calculator import Methods,Calculator

def creat_plot_3D(X,Y,Z,ref,save=False):
    plt.figure(dpi=100)
    ax1=plt.axes(projection='3d')
    ax1.plot_surface(X,Y,Z,cmap='rainbow')
    ax1.plot_surface(X,Y,np.ones(Z.shape)*ref)
    if save:
        if type(save) is str:
            plt.savefig(save)
        else:
            plt.savefig('result.png')
    else:
        plt.show()


@jit(nopython=True)
def check_values(numbers1,numbers2,ref1,ref2,error1,error2):
    result=[]
    judge=False
    for i in range(numbers1.shape[0]):
        for j in range(numbers1.shape[1]):
            if abs(numbers1[i,j]-ref1)<error1:
                if abs(numbers2[i,j]-ref2)<error2:
                    result.append([i,j])
                    judge=True
    if judge:
        print('True')
        return result
    else:
        print('False')



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
    result=np.ones((100,100,2))
    X,Y=np.meshgrid(ST_test,g_out_test)
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            result[i,j]=cal.me.cal_ML_simple_B(X[i,j],L_s,True,Y[i,j])[1:3]
    creat_plot_3D(X,Y,result[:,:,0],0.0)
    creat_plot_3D(X,Y,result[:,:,1],cal.me.con.g_in)