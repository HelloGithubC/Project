'''
@Author: Xiao Liang
@Date: 2019-09-08 23:02:20
@LastEditors: Xiao Liang
@LastEditTime: 2019-09-22 20:57:46
@Description: file content
'''
from openpyxl import Workbook,load_workbook
import datetime
class Test:
    @classmethod
    def write_excel(cls,filename,value1,value2):
        wb=Workbook()
        index_value=0
        while True:
            if 'result:{0}'.format(index_value) in wb.sheetnames:
                index_value+=1
            else:
                break
        file='result_{0}'.format(index_value)
        wb.create_sheet(file,index=index_value)
        sheet=wb[wb.sheetnames[index_value]]
        sheet.append([str(datetime.datetime.now())])
        sheet.append(['r1','P1','T1','M1','L1','g1','dg1'])
        for i in value1:
            sheet.append(i)
        sheet.append([])
        sheet.append(['r2','P2','T2','M2','L2','g2','dg2'])
        for i in value2:
            sheet.append(i)
        wb.save(filename)