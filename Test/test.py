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
    te.creat_t([M_p,ST,L_s],test_start=0)"""

    te=Test()
    te.find_L(39)
