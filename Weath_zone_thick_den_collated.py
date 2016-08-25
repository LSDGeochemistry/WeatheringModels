##Weath_zone_thick_den_collated.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## Outputs the weatrerhing zone thickness as a function of the denudation rate for various authors 
##(Gabet and Mudd,2009; West,2012; Maher and Chamberlain, 2014; Lebedeva et al 2010)
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## LK 19/08/2016
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

import numpy as np
import matplotlib.pyplot as pp
## Erosion rate (in m/yr)
E = np.arange(0.000001,0.01,0.00000001)
#convert to mm/yr
e = E*1000
#Density (t km-3) (from Gabet and Mudd 2009)
rho = 1.3*np.power(10.0,9)
## Depth calculation from West, 2012 and Montgomery and Brandon 2002
#Empirical constants (see Montgomery and Brandon, 2002 figure 5)
a = 1.4*np.power(10.0,-6)
b = 1.8

#thickness in m
z = (E/1000*rho/a)**(1/b)*100000/rho

##Depth calculation from Maher and Chamberlain, 2014 and Gabet and Mudd, 2009
#thickness in m
h = np.log((E/1000*rho/10**4.01))/-2300*1000

## Depth Calculatiom from Lebedeva et al 2010
###Variables - values taken from Lebedeva et al (2010)
#Erosion rate converted to m s-1
W = E/60/60/24/365
#Darcy velocity of pore fluid through regolith (in m a-1)
u = 0.2
#Darcy velocity converted to m s-1
ums = u/365/24/60/60
#Reduced Diffusivity in pore fluid (Diffusivity/tortuosity)
D = (0.8*np.power(10.0,-9))
#Porosity
phi = 0.2
#Porosity of ALbite
phi_ab = 0.4
#Specific volume of albite (m3 mol-1)
V_ab = 1*np.power(10.0,-4)
#Dissolution Rate of albite (mol m-2 s-1)
k_ab = 3.87*np.power(10.0,-10)
#Initial volume fraction of albite
phi_ab = 0.4
#Specific surface area of albite
s_ab = 3.5*np.power(10.0,4)
#Fitting factor for the allbite dissolution rate
Psi_ab = 1.2*np.power(10.0,-2)
#Concentration of aqueous compontent in equilibrium with bedrock minerals
Ce = 0.2
#Concentration of aqueous compnonet at inlet
Cr = 0
# Nucleation barrier for Kaolinite precipitation
Cl = 0.9*Ce 
#Weathering zone thickness with erosion (equation 19)
def weath_zone_thick_eq(ums, phi, phi_ab, D , W,Cr,Cl, V_ab):

  leading_term_eq = -((D*phi)/(ums-phi*W))
  log_term_eq = 1-(((ums-phi*W)*(Cl-Cr))/(phi_ab/V_ab*W))
  #print leading_term_eq
  #print log_term_eq
  delta_eq = leading_term_eq*np.log(log_term_eq)
  return delta_eq  
#Equation 19
delta_eq = weath_zone_thick_eq(ums, phi, phi_ab, D ,W ,Cr,Cl, V_ab)

### Estimating values used for tau_kj in Li et al 2014
#Method one
tau_kj_mb = z/E
#Method two
tau_kj_gm = h/E
#Method three
tau_kj_l = delta_eq/E

###Finding the erosion rate and depth from a given value of tau_kj
tau_kj_value = 10000
##Method one





fig = pp.figure()
ax = fig.add_subplot(1,2,1)
pp.tick_params(axis='both', which='major', labelsize=8)
ax.plot(e,z,'k')
ax.plot(e,h,'r')
ax.plot(e,delta_eq,'g')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0,max(e))
ax.set_ylim(0.00001,max(z))
ax.set_xlabel('Erosion Rate ($m\ a^{-1}$)',fontsize = '10')
ax.set_ylabel('Weathering Zone Thickness ($m$)',fontsize = '10')
ax2 = fig.add_subplot(1,2,2)
pp.tick_params(axis='both', which='major', labelsize=8)
ax2.plot(e, tau_kj_mb, 'k--')
ax2.plot(e, tau_kj_gm, 'r--')
ax2.plot(e, tau_kj_l, 'g--')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_ylim(0.01,max(tau_kj_gm))
ax2.set_xlabel('Erosion Rate ($m\ a^{-1}$)',fontsize = '10')
ax2.set_ylabel('Kinetic Residence Time for Mineral j ($a^{-1}$)',fontsize = '10')
pp.tight_layout()
pp.show()