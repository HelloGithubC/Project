import numpy as np
import matplotlib.pyplot as plt
from xlrd import open_workbook
import re

class Plot(object):
	"""Used to plot."""
	def __init__(self,dpi=100):
		self.dpi=dpi

	def plot_from_memorry(self, vars, labels=None, xylabel=None):
		"""The first element of vars is independent variable."""
		length=len(vars)-1
		vars=np.array(vars)
		x=vars[0]

		fig=plt.figure(dpi=self.dpi)
		if labels!=None:
			for i in range(length):
				plt.loglog(x,vars[i+1],label=labels[i],lw=2)
		else:
			for i in range(length):
				plt.loglog(x,vars[i+1],lw=2)
		
		if xylabel!=None:
			plt.xlabel=xylabel[0]
			plt.ylabel=xylabel[1]
		plt.show()

	def plot_from_file(self, filename, dtype, vars_number,labels=None,xylabels=None):
		"""The order is the same as it in vars."""
		if dtype=='excel':
			workshop=open_workbook(filename)
			sheet=workshop.sheet_by_index(0)
			vars=[]
			for i in range(vars_number):
				store=sheet.col(i)
				for j in range(len(store)):
					store[j]=store[j].value
				vars.append(store)
			self.plot_from_memorry(vars,labels,xylabels)
		elif dtype=='txt':
			file=open(filename)
			store=file.readline()
			store=re.sub(r'\D',' ',store).split()
			length=vars_number
			vars=[[] for i in range(length)]
			for i in range(length):
				vars[i].append(store[i])
			while True:
				store=file.readline()
				if not store:
					break
				store=re.sub(r'\D',' ',store).split()
				for i in range(length):
					vars[i].append(store[i])
			self.plot_from_memorry(vars,labels,xylabels)
		else:
			raise ValueError('The value of dtype is invalid.')