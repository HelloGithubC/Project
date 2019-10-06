class Const:
	def __init__(self):
		self.M_e=5.965e27 #g
		self.M_p=5.5  #/M_e
		self.M_c=5      #/M_e
		self.G=6.674e-8   #erg*cm**2/g**2
		self.Myr=3.156e13

		self.a=0.5 #/AU
		self.c_0=1e5*self.a**(-1/4)
		self.sigma=5.6704e-5 #erg/cm**2/s**1/K**4
		self.pi=3.14159
		self.alpha=1.
		self.beta=1.
		self.k_0=1e-3    #Opacity
		self.R=8.314e7   #Universal gas constant
		self.T_0=300*self.a**(-1/2)
		self.mu=self.R*self.T_0/self.c_0**2       #Relative molecular mass of gases
		self.gamma=7./5
		self.g_ad=(self.gamma-1)/self.gamma #gradient when ad

		self.rho_in=7.
		self.rho_0=2.4e-9*self.a**(-11/4)
		self.R_p=(3*self.M_c*self.M_e/(4*self.pi*self.rho_in))**(1/3)
		self.R_H=self.a*2e11*self.M_p**(1/3)
		#R_B=4e10*a**(1/2)*M_p
		self.R_B=self.G*self.M_p*self.M_e/self.c_0**2
		#R_out=min(R_H,R_B)
		self.R_out=self.R_B
		self.P_0=self.rho_0*self.R*self.T_0/self.mu
		#calculate sigma
		self.sigma_1=9e9*906.4*(self.T_0/1002)**(5./4.)*(self.T_0/1500)**(
			-1./2.)*(self.P_0/1e6)**(-1./2.)
		self.sigma_2=(4.35*0.5)/(self.T_0*8.617e-5)
		self.k_P=self.P_0/1e6

		self.Lambda=1

		self.g_in=183870671770695.2 # When ST = 0.20143836115164995, L_s = 13888.2056137455 and g=1e-9. Use calculator.cal_g_in to get it.

	def cal_const(self):
		self.c_M=4*self.pi*self.rho_0*self.R_out**3/self.M_e
		self.c_T=3*self.k_0*self.P_0/(
                64*self.pi*self.G*self.sigma*self.M_e*self.T_0**4)
		self.c_P=-(self.rho_0/(self.P_0*self.R_B))*self.G*self.M_e
		return [self.c_P,self.c_T,self.c_M]

	def set_M_p(self,M_p):
		self.M_p=M_p
		self.R_H=self.a*2e11*self.M_p**(1/3)
		self.R_B=self.G*self.M_p*self.M_e/self.c_0**2
		self.R_out=self.R_B
