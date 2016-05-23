## Weathering_zone_thicknes_erosions.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## Makes plots of weathering zone thickness as a function of erosion rate following Lebedeva et al 2010 ESPL
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## LK 23/05/2016
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
import matplotlib.pyplot as pp
import numpy as np
import matplotlib.lines as mpllines
from matplotlib.ticker import FormatStrFormatter
  
#Weathering zone thickness with erosion (equation 19)
def weath_zone_thick_eq(ums, phi, phi_ab, D , W,Cr,Cl, V_ab):

  leading_term_eq = -((D*phi)/(ums-phi*W))
  log_term_eq = 1-(((ums-phi*W)*(Cl-Cr))/(phi_ab/V_ab*W))
  #print leading_term_eq
  #print log_term_eq
  delta_eq = leading_term_eq*np.log(log_term_eq)
  return delta_eq  

#Erosion Rate (m ka -1)   
w = np.arange(0.01,2,0.0005) 
#Erosion rate converted to m s-1 
W = w/60/60/24/365/1000


###Variables - values taken from Lebedeva et al (2010)
#Darcy velocity of pore fluid through regolith (in m s-1)
u = 0.2
#Darcy velocity converted to m s-1
ums = u/365/24/60/60
#Reduced Diffusivity in pore fluid (Diffusivity/tortuosity)
D = 0.8*np.power(10.0,-9)
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
#Equation 19
delta_eq = weath_zone_thick_eq(ums, phi, phi_ab, D , W,Cr,Cl, V_ab)
#PLotting the results  
fig = pp.figure(1, facecolor='white',figsize=(10,7.5))
ax1 = fig.add_subplot(1,1,1)  
  
pp.plot(w,delta_eq,'k',linewidth=3)


pp.rcParams['xtick.direction'] = 'out'
pp.rcParams['ytick.direction'] = 'out'
 
ax1.spines['top'].set_linewidth(2.5)
ax1.spines['left'].set_linewidth(2.5)
ax1.spines['right'].set_linewidth(2.5)
ax1.spines['bottom'].set_linewidth(2.5) 
ax1.tick_params(axis='both', width=2.5)    
ax1.xaxis.set_major_formatter(FormatStrFormatter('%.01f'))
ax1.set_yscale('log')

for line in ax1.get_xticklines():
   line.set_marker(mpllines.TICKDOWN)

for line in ax1.get_yticklines():
    line.set_marker(mpllines.TICKLEFT)

pp.ylabel('Weathering zone thickness ($m$)',fontsize = 20)
pp.xlabel('Erosion Rate, $m$ $ka^{-1}$',fontsize = 20) 

pp.tight_layout()
       
pp.show()
