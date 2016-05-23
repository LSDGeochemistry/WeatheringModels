##Maher_flux_erosion_rate.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## Makes plots of SiO2 as a function of erosion rate for different runoffs following Maher and Chamberlain 2014 Science
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## LK 23/05/2016
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


import matplotlib.pyplot as pp
import numpy as np
from matplotlib import rcParams
import matplotlib.lines as mpllines
from matplotlib.ticker import FormatStrFormatter
import MaherVariables as MV




# These set the font size
label_size = 20
#title_size = 30
axis_size = 28

rcParams['font.size'] = label_size 
    
#get the fw as a function of the erosion rate
# depth  of weathering zone (m)
H = 5
#Specific surface area ()
A = 0.1
log_E = np.linspace(-6,-3,100)
E = np.power(10,log_E)*1000
L = 50
q = np.zeros(5)
q[0] = 1
q[1] = 2
q[2] = 4
q[3] = 8
q[4] = 10
  
#print E

Dw = MV.calculate_Dw(H,E/1000,A,L)
nDw = Dw.shape[0]
nq = q.shape[0]
Q = np.zeros((nq,nDw))   
#print Dw
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
  
  
  
pp.plot(E,Q[0],linewidth=3,label = ("q = "+str(q[0])))
pp.plot(E,Q[1],linewidth=3,label = ("q = "+str(q[1])))
pp.plot(E,Q[2],linewidth=3,label = ("q = "+str(q[2])))
pp.plot(E,Q[3],linewidth=3,label = ("q = "+str(q[3])))
pp.plot(E,Q[4],linewidth=3,label = ("q = "+str(q[4]))) 
#pp.plot(H,fw[:,3],linewidth=3,label = ("E = "+str(E[3]*1000)+ "$mm/yr$"))
pp.legend(loc=2,fontsize='22')
pp.rcParams['xtick.direction'] = 'in'
pp.rcParams['ytick.direction'] = 'in'
ax1.set_xlim(min(E),max(E))
ax1.set_xscale('log')
ax1.set_yscale('log')
  

  
ax1.spines['top'].set_linewidth(2.5)
ax1.spines['left'].set_linewidth(2.5)
ax1.spines['right'].set_linewidth(2.5)
ax1.spines['bottom'].set_linewidth(2.5) 
ax1.tick_params(axis='both', width=2.5)    
ax1.xaxis.set_major_formatter(FormatStrFormatter('%.03f'))
for line in ax1.get_xticklines():
    line.set_marker(mpllines.TICKUP)

for line in ax1.get_yticklines():
    line.set_marker(mpllines.TICKRIGHT)

pp.xlabel('Erosion Rate ($mm\ a^{-1}$)',fontsize = axis_size)
pp.ylabel('SiO$_2$ (aq) flux ($tonnes\ km^2\ yr^{-1}$)',fontsize = axis_size) 
#pp.title("$A$ is "+str(A)+ " $m^2/g$")
#pp.ylim(0,1)
pp.tight_layout()
pp.show()



 