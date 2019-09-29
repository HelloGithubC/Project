import numpy as np 
import matplotlib.pyplot as plt 
from calculator import Methods,Calculator

if __name__ == "__main__":
    cal=Calculator()
    cal.read_excel('Data.xlsx')
    m,ST_init,L_s=cal.M_p[0],cal.ST[0],cal.L_init[0]
    ST,g_in=cal.cal_g_in(m,ST_init,L_s,1e-9,0.0)
    print(cal.me.cal_ML_simple_B(ST,L_s,True))
    print(g_in)