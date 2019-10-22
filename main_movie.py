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
        RCB_index[i]=cal.me.RCB_index
        RCB[i]=r[i,int(RCB_index[i])]
    
    fig=plt.figure(figsize=(20,10))
    ax1=fig.add_subplot(121)
    line1=ax1.loglog(r[0],P[0],label='r-P')[0]
    line2=ax1.loglog(r[0],T[0],label='r-T')[0]
    RCB_index_need=int(RCB_index[0])
    point1=ax1.plot([r[0,RCB_index_need]],[P[0,RCB_index_need]],'ro')[0]
    point2=ax1.plot([r[0,RCB_index_need]],[T[0,RCB_index_need]],'go')[0]
    ax1.grid(ls='--')

    ax2=fig.add_subplot(122)
    ax2.set_ylim([4.9,12.3])
    line3=ax2.semilogx(r[0],M[0],label='r-M')[0]
    point3=ax2.plot([r[0,RCB_index_need]],[M[0,RCB_index_need]],'wo')[0]
    ax2.grid(ls='--')

    def update(num):
        line1.set_data(r[num],P[num])
        line2.set_data(r[num],T[num])
        RCB_index_need=int(RCB_index[num])
        point1.set_data([r[num,RCB_index_need]],[P[num,RCB_index_need]])
        point2.set_data([r[num,RCB_index_need]],[T[num,RCB_index_need]])
        line3.set_data(r[num],M[num])
        point3.set_data([r[num,RCB_index_need]],[M[num,RCB_index_need]])
        return [line1,line2,point1,point2,line3,point3]
    
    ani=animation.FuncAnimation(fig,update,np.arange(500),interval=40, blit=True)

    ani.save('test.mp4',writer='ffmpeg',fps=25)
