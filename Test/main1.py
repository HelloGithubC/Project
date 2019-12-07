import numpy as np 
import matplotlib.pyplot as plt
from calculator import Methods,Calculator
from scipy.optimize import leastsq

def func(p,x):
    k,b=p
    return k*x+b 
def error(p,x,y):
    return func(p,x)-y/1e+13 

if __name__=='__main__':
    cal=Calculator()
    M_p=np.linspace(5.15,12,200)
    cal.creat_data(M_p,'Data2.xlsx','test.png')