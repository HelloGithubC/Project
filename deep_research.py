import numpy as np 
import matplotlib.pyplot as plt 
from calculator import Methods

me=Methods()
me.con.set_M_p(10.7)
result=me.find_L_M([5.0,0.0],[0.0,1.0],[1e-8,1e-3],[1e-5,3e-3])
print(result)