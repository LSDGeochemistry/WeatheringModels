## Chemical_depletion_factor_with_erosion_na_si.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## Makes plots of the CDF for Na and Si as a function of denudation rate following Lebedeva et al 2010 ESPL
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## LK 23/05/2016
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
import matplotlib.pyplot as pp
import numpy as np
import matplotlib.lines as mpllines

#Variables - values taken from Lebedeva et al (2010)
#Reduced Diffusivity in pore fluid
D = 0.8*np.power(10.0,-9)
#Porosity
phi = 0.2
#Darcy velocity of pore fluid through regolith (in m a-1)
u = 0.15
#Darcy velocity of pore fluid through regolith (in m s-1)
ums = u/365/24/60/60
#Initial volume fraction of albite
phi_ab = 0.4
#Specific volume of albite (m3 mol-1)
V_ab = 1*np.power(10.0,-4)
#Specific surface area of albite
s_ab = 3.5*np.power(10.0,4)
#Fitting factor for the allbite dissolution rate
Psi_ab = 1.2*np.power(10.0,-2)
#Concentration of aqueous compontent in equilibrium with bedrock minerals
Ce = 0.2
#Concentration of aqueous compnonet at inlet
Cr = 0
#Initial concentration of Albite mol m-3 
Q_ab = phi_ab/V_ab
#Initial concentration of Inert mineral mol m-3
Q_in = 0.7/22.7*np.power(10.0,-6)
#Kinetic Constant (s)
k = 2.5*np.power(10.0,-10)
#Erosion Rate (m/a)
W = np.arange(0,0.001,0.00001)
#Erosion Rate (m/s)
w = W/365/24/60/60

#Defining beta, equation 24 Lebedeva et al, (2010)
def calculate_beta(k, D, phi, ums):
  beta = np.sqrt(1+((4*k*D*phi)/(ums*ums)))
  return beta 
#Defining eta, equation 39 Lebedeva et al (2010)
def calculate_eta(ums,Ce,Cr,beta,Q_ab,w): 
  beta = calculate_beta(k, D, phi, ums)
  eta = (ums*((Ce-Cr)*(beta+1)))/(2*Q_ab*w)
  return eta

#Chemical Depletion Factor, equation 40 Lebedeva et al (2010)
beta = calculate_beta(k, D, phi, u) 

eta_Na = calculate_eta(u,Ce,Cr,beta,Q_ab,W)
eta_Si = (2.0/3.0)*eta_Na
#Convert the erosion rate to denudation rate, equation 35 Lebedeva et al (2010)
Den = W*(Q_ab*262+Q_in*60)
#Plotting the Results     
fig = pp.figure(1, facecolor='white',figsize=(10,7.5))
ax1 = fig.add_subplot(1,1,1)  


pp.plot(Den,eta_Na,'k', linewidth=3, label='Na')
pp.plot(Den,eta_Si,'r', linewidth=3, label='Si')
  
pp.rcParams['xtick.direction'] = 'out'
pp.rcParams['ytick.direction'] = 'out'

  
  
ax1.spines['top'].set_linewidth(2.5)
ax1.spines['left'].set_linewidth(2.5)
ax1.spines['right'].set_linewidth(2.5)
ax1.spines['bottom'].set_linewidth(2.5) 
ax1.tick_params(axis='both', width=2.5)    
ax1.set_xlim(0,1000)
ax1.set_ylim(0,1)
ax1.legend()

for line in ax1.get_xticklines():
    line.set_marker(mpllines.TICKDOWN)
for line in ax1.get_yticklines():
    line.set_marker(mpllines.TICKLEFT)

pp.ylabel('Chemical Depletion Factor',fontsize = 20)
pp.xlabel('Denudation Rate ($tonnes/km^2/yr$)',fontsize = 20) 


pp.tight_layout()
      
pp.show()    
  
    
    
    
    