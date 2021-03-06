##Maher_clay_flux_erosion_rate.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## Makes plots of clay flux as a function of erosion rate for different runoffs following Maher and Chamberlain 2014 Science, 
## weathering zone thickness is calculated as a funciton of erosion rate and Flowpath length is either asumed to be equal
## to one or to weathering zone thickness.
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## LK 24/08/2016
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


import matplotlib.pyplot as pp
import numpy as np
from matplotlib import rcParams
import matplotlib.lines as mpllines
import MaherVariables as MV




# These set the font size
label_size = 20
#title_size = 30
axis_size = 20

rcParams['font.size'] = label_size 
    


#Specific surface area (m2 g-1)
A = 0.1
# Erosion rate m/year
E = np.arange(0.000001,0.01,0.00000001)
#Soil Density tonne km-3
rho = 1.3*np.power(10.0,9)
#Denudation rate tonnes km-2 yr
Den = E/1000*rho
#Weathering Zone thickness (m) as a function of erosion rate
H = np.log((E/1000*rho)/10**4.01)/-2300*1000
H[H < 0] = 0        
#Flowpath length, assumed to be either equal to weathering zone thickness (h) or one, see Maher and Chamberlain 2014 suppl. mat. fig. 5  
L = H
#Different values of flow rate (m a-1) used, see Maher and Chamberlain 2014 suppl. mat. fig. 5
q = np.zeros(5)
q[0] = 0.1
q[1] = 0.3
q[2] = 1
q[3] = 3
q[4] = 10
  

#Calculate the Damkohler coefficent using the formula in MaherVariables.py
Dw = MV.calculate_Dw(H,E,A,L)
nDw = Dw.shape[0]
nq = q.shape[0]
Q = np.zeros((nq,nDw))   

#Find Q for each value of q defined above
i = 0
for thisDw in Dw:
      j = 0
      for thisq in q:
          #print "H: " + str(thisH)+ " E: " + str(thisE)
          #print "i: "+str(i)+ " j: " +str(j)
          Q[j][i]=MV.calculate_Q(thisDw,thisq)
          j = j+1
      i = i+1

  
# now plot the results    
fig = pp.figure(1, facecolor='white',figsize=(10,7.5))
ax1 = fig.add_subplot(1,1,1)  

#Removes negative flux values for plotting  
E = E[Q[0] > 0]
Den = E/1000*rho    
Q_0 = Q[0][Q[0] > 0]
Q_1 = Q[1][Q[1] > 0]  
Q_2 = Q[2][Q[2] > 0]  
Q_3 = Q[3][Q[3] > 0]  
Q_4 = Q[4][Q[4] > 0]  


molar_SiO2 = 2*15.9994+28.0855
clay_Q_0 = Q_0/molar_SiO2/2.8*256.8
clay_per_0 = clay_Q_0/Den*100 
clay_Q_1 = Q_1/molar_SiO2/2.8*256.8
clay_per_1 = clay_Q_1/Den*100 
clay_Q_2 = Q_2/molar_SiO2/2.8*256.8
clay_per_2 = clay_Q_2/Den*100 
clay_Q_3 = Q_3/molar_SiO2/2.8*256.8
clay_per_3 = clay_Q_3/Den*100
clay_Q_4 = Q_4/molar_SiO2/2.8*256.8
clay_per_4 = clay_Q_4/Den*100



pp.plot(Den,clay_Q_0,linewidth=3,label = ("q = "+str(q[0])))
pp.plot(Den,clay_Q_1,linewidth=3,label = ("q = "+str(q[1])))
pp.plot(Den,clay_Q_2,linewidth=3,label = ("q = "+str(q[2])))
pp.plot(Den,clay_Q_3,linewidth=3,label = ("q = "+str(q[3])))
pp.plot(Den,clay_Q_4,linewidth=3,label = ("q = "+str(q[4]))) 

pp.legend(loc=1,fontsize='20')
pp.rcParams['xtick.direction'] = 'in'
pp.rcParams['ytick.direction'] = 'in'
ax1.set_xlim(min(Den),10000)
ax1.set_ylim(0,100)
#ax1.set_xscale('log')
  

  
ax1.spines['top'].set_linewidth(2.5)
ax1.spines['left'].set_linewidth(2.5)
ax1.spines['right'].set_linewidth(2.5)
ax1.spines['bottom'].set_linewidth(2.5) 
ax1.tick_params(axis='both', width=2.5)    

for line in ax1.get_xticklines():
    line.set_marker(mpllines.TICKUP)

for line in ax1.get_yticklines():
    line.set_marker(mpllines.TICKRIGHT)

pp.xlabel('Denudation Rate  (($tonnes\ km^2\ yr^{-1}$))',fontsize = axis_size)
pp.ylabel('Halloysite (clay) flux ($tonnes\ km^2\ yr^{-1}$)',fontsize = axis_size) 

pp.tight_layout()
  
pp.show()



 

