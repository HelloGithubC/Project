import numpy as np 
import matplotlib.pyplot as plt 
from scipy import interpolate 

class DataBase(object):
    f_sigma=[]
    @classmethod
    def creat_sigma(cls):
        P=np.array([280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,447,464,470,474,477,480,483,488,494,500,505,507,512,520,528,534,540,546,556,566,575,587,600,612,624,634,655,678,700,720,740,756,763,769,776,787,800,805,810,815,820,825,830,840,846,850,860,870,875,880,890,900,910,918,930,945,960,970,975,980,985,990,995,1000])
        ln_sigma=np.array([1321,1327,1334,1340,1348,1354,1360,1366,1372,1378,1385,1390,1396,1402,1409,1414,1407,1401,1396,1393,1388,1384,1378,1367,1352,1338,1332,1317,1290,1262,1243,1220,1198,1171,1154,1149,1156,1168,1174,1176,1179,1193,1208,1219,1232,1243,1247,1246,1245,1244,1237,1220,1204,1194,1180,1172,1150,1127,1097,1075,1064,1040,1012,1000,990,970,952,935,925,909,892,875,860,853,850,840,828,800,760])
        
        P=np.array([-6+2*(P[i]-145)/149 for i in range(len(P))])
        ln_sigma=np.array([-6+(1414-ln_sigma[j])/100 for j in range(len(ln_sigma))])
        cls.f_sigma=interpolate.interp1d(P,ln_sigma,kind='cubic')
