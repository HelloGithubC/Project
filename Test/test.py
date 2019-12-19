from valueDataBase import DataBase
from calculator import Calculator
from const import Const
from test1 import Test 
from openpyxl import Workbook,load_workbook

if __name__ == '__main__':
    """cal=Calculator()
    cal.read_excel('Data2.xlsx')
    M_p=cal.M_p
    L_s=cal.L_init
    ST=cal.ST
    te=Test()
    te.creat_t([M_p,ST,L_s])"""

    cal=Calculator()
    cal.me.change_const(Const(0.3))
    cal.read_excel('Data2.xlsx')
    order=20
    m=cal.M_p[order]
    cal.me.con.set_M_p(m)
    L_s=cal.L_init[order]
    ST=cal.ST[order]

    te=Test()
    te.con.set_M_p(m)
    print(te.cal_ML_simple_B(ST,L_s,True,0.0,-1))
    result=te.find_L_simple(0.0,[ST,L_s],1e-8,1e-3)
    print(result)
    print(te.cal_ML_simple_B(result[0],L_s,True,0.0,-1))