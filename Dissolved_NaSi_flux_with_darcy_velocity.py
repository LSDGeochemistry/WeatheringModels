## Dissolved_Na_Si_flux_with_darcy_velocity.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## Makes plots of dissolved NaSi as a function of Darcy velocity following Lebedeva et al 2010 ESPL
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## LK 23/05/2016
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
import matplotlib.pyplot as pp
import numpy as np
import matplotlib.lines as mpllines
from matplotlib.ticker import FormatStrFormatter


###Variables - values taken from Lebedeva et al (2010)
#Reduced Diffusivity in pore fluid (Diffusivity/tortuosity)
D = 0.8*np.power(10.0,-9)
#Porosity
phi = 0.2
#Porosity of ALbite
phi_ab = 0.4
#Specific volume of albite (m3 mol-1)
V_ab = 1*np.power(10.0,-4)
#Darcy velocity of pore fluid through regolith (in m/a)
ua= np.arange(0,1,0.001)
#Darcy Velocity (m/s)
u = ua/365/60/24/24
#Dissolution Rate of albite (mol m-3 s-1)
k_ab = 2.51*np.power(10.0,-10)
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
##Initial concentration of Albite mol/m-3 (see White et al 1998?)
Q_ab = phi_ab/V_ab
#Initial concentration of Inert mol/m-3
Q_in = 0.7/22.7*np.power(10.0,-6)
###Changing Variable
#Erosion Rate (m/ka)
w = 0.5
W = w/1000/365/24/60/60
#Defining k, equation 4 Lebedeva et al (2010)
def calculate_k(k_ab, phi_ab, s_ab, Psi_ab):
  k = 2*k_ab*phi_ab*s_ab*Psi_ab
  return k 
  
#Defining beta, equation 24 Lebedeva et al, (2010)
def calculate_beta(k, D, phi, u):
  beta = np.sqrt(1+4*k*D*phi/(u*u))
  return beta
 
k = calculate_k(k_ab, phi_ab, s_ab, Psi_ab)
beta = calculate_beta(k, D, phi, u) 
#Calculating the NaSi2 flux for kinetic limited (equation 25)
omega = (u*Cr+(((u-phi*W)*(1-beta)*(Ce-Cr))/2))*262*365*24*60*60

 # now plot the results    
fig = pp.figure(1, facecolor='white',figsize=(10,7.5))
ax1 = fig.add_subplot(1,1,1)  

pp.plot(ua,omega,'k', linewidth=3)

pp.rcParams['xtick.direction'] = 'out'
pp.rcParams['ytick.direction'] = 'out'

  
ax1.spines['top'].set_linewidth(2.5)
ax1.spines['left'].set_linewidth(2.5)
ax1.spines['right'].set_linewidth(2.5)
ax1.spines['bottom'].set_linewidth(2.5) 
ax1.tick_params(axis='both', width=2.5)    
ax1.xaxis.set_major_formatter(FormatStrFormatter('%.01f'))

for line in ax1.get_xticklines():
    line.set_marker(mpllines.TICKDOWN)
for line in ax1.get_yticklines():
    line.set_marker(mpllines.TICKLEFT)

pp.ylabel('Dissolved $NaSi_2$ flux ($tonnes\ km^{-2} yr^{-1}$)',fontsize = 20)
pp.xlabel('Darcy Velocity ($m\ yr^{-1}$)',fontsize = 20) 

pp.tight_layout()
      
  #pp.savefig("Particle_size_thick.svg", format='svg')  
pp.show()    
  
    
    
    
    