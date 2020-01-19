import numpy as np 
import matplotlib.pyplot as plt 
from calculator import Methods,Calculator
from const import Const
from test1 import Test

cal=Calculator()
cal.read_excel('Data2.xlsx')

te=Test()

start,end=198,199
result=np.ones((end+1-start,5))
for order in range(start,end+1,1):
    m=cal.M_p[order]
    ST=cal.ST[order]
    L_s=cal.L_init[order]
    te.con.set_M_p(m,1e+6)
    ST,L_c=te.find_L_two(0.0,L_s,0.0,ST,1e-4,True)
    result[order-start]=[int(order),te.P[-1],te.T[-1],L_s,ST]
    print("{0} is done.".format(order))