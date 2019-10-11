import numpy as np
import matplotlib.pyplot as plt
from calculator import Methods,Calculator
import matplotlib.animation as animation

if __name__ == "__main__":
    cal=Calculator()
    cal.read_excel('Data.xlsx')
    r=np.ones((500,800))
    P=np.ones((500,800))
    T=np.ones((500,800))
    M=np.ones((500,800))
    RCB=np.ones(500)
    RCB_index=np.ones(500)
    for i in range(500):
        M_p=cal.M_p[i]
        cal.me.con.set_M_p(M_p)
        cal.me.cal_ML_simple(cal.ST[i],cal.L_init[i])
        r[i]=cal.me.r
        P[i]=cal.me.P
        T[i]=cal.me.T
        M[i]=cal.me.M
        RCB[i]=cal.me.RCB
        for j in range(800):
            if abs(r[i,j]-RCB[i])<1e-3:
                RCB_index[i]=j
                break
    fig=plt.figure(figsize=(10,20))
    ax1=fig.add_subplot(211)
    line1=ax1.loglog(r[0],P[0],label='r-P')[0]
    line2=ax1.loglog(r[0],T[0],label='r-T')[0]
    point1=ax1.scatter([r[RCB_index]],[P[0,RCB_index]],marker='o')[0]
    point2=ax1.scatter([r[RCB_index]],[T[0,RCB_index]],marker='o')[0]

    ax2=fig.add_subplot(212)
    line3=ax2.semilogx(r[0],M[0],label='r-M')[0]
    point3=ax2.scatter([r[RCB_index]],[M[0,RCB_index]],marker='o')[0]

    def update(num):
        line1.set_data(r[num],P[num])
        line2.set_data(r[])
