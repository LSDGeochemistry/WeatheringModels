## Dissolved_Na_Si_flux_with_erosion.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## Makes plots of dissolved NaSi as a function of erosion rate following Lebedeva et al 2010 ESPL
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## LK 23/05/2016
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

import matplotlib.pyplot as pp
import numpy as np
import matplotlib.lines as mpllines


### How to calulate the NaSi2 flux in equilibrium conditions, using the diffusion-advective model

###Variables - values taken from Lebedeva et al (2010)
#Specific volume of albite (m3 mol-1)
V_ab = 1*np.power(10.0,-4)
#Darcy velocity of pore fluid through regolith (in m/a)
u = 0.15
#Initial volume fraction of albite
phi_ab = 0.1
#Concentration of aqueous compontent in equilibrium with bedrock minerals
Ce = 0.2
# Nucleation barrier for Kaolinite precipitation
Cl = 0.9*Ce 
##Initial concentration of Albite mol/m-3 
Q_ab = phi_ab/V_ab
#Initial concentration of Inert mol/m-3
Q_in = 0.7/(22.7*np.power(10.0,-6))
#Erosion Rate (m/ka)
w = np.arange(0,0.5,0.0001)
#Erosion Rate (m/a)
W = w/1000
#Denudation Rate (tonnes/km2/yr) (equation 35)
Den = W*(Q_ab*262+Q_in*60)
#Calculating the flux for local equilibrium (equation 21)
omega_eq = ((u*Cl)-(W*Q_ab))*262
#Plot the results    
fig = pp.figure(1, facecolor='white',figsize=(10,7.5))
ax1 = fig.add_subplot(1,1,1)  
pp.plot(Den,omega_eq,'k', linewidth=3)
pp.rcParams['xtick.direction'] = 'out'
pp.rcParams['ytick.direction'] = 'out'

ax1.spines['top'].set_linewidth(2.5)
ax1.spines['left'].set_linewidth(2.5)
ax1.spines['right'].set_linewidth(2.5)
ax1.spines['bottom'].set_linewidth(2.5) 
ax1.tick_params(axis='both', width=2.5)    
ax1.set_xlim(0,1000)

for line in ax1.get_xticklines():
    line.set_marker(mpllines.TICKDOWN)
for line in ax1.get_yticklines():
    line.set_marker(mpllines.TICKLEFT)

pp.ylabel('Dissolved $NaSi_2$ flux ($tonnes\ km^2 yr^{-1}$)',fontsize = 20)
pp.xlabel('Denudation Rate  ($tonnes\ km^2 yr^{-1}$) ',fontsize = 20) 

pp.tight_layout()
      
pp.show()    
  
    
    
    
    