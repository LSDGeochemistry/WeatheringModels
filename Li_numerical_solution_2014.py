##Li_numerical_solution_2014.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## Outputs dimensionless concentration of mineral j and solute i as a function of depth, time, and residence times using a finite difference model
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## LK 24/05/2016
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

import numpy as np

###set the step sizes, taken from Li et al (2014), solution should be grid-independent
#dimensionless depth steps
d_epsilon = 1*np.power(10.0,-2)
#dimensionless time steps
d_theta = 1*np.power(10.0,-2)

###Set the residence times used to calculate Damkohler numbers
#bulk rock residence time
tau_r = 1*np.power(10.0,8)
#bulk water residence time
tau_w = 1*np.power(10.0,1)
#kinetic residence time for mineral j
tau_k = 1*np.power(10.0,4)

### Set other variables
#saturation capacity index for mineral j and solute i
alpha_ij =1*np.power(10.0,-4)

###Calulate Damkohler numbers
Da_r = tau_r/tau_k
Da_w = tau_w/tau_k

#Calculate the Breakthrough time for mineral j (using the dimensionless reaction front velocity)
lambda_j = 1/((alpha_ij-Da_w/Da_r)/(alpha_ij+1))
#Calulate the time interval used as a function of the time steps and the breakthrough time, assuming that lambda_j > 0 (see Li et al, 2014 p 35)
d_thetan = d_theta*lambda_j

###create the matrices
space_nodes = int(1/d_epsilon)
time_nodes = int(1/d_theta)


###populate with zeroes
eta_matrix_i = np.zeros((time_nodes,space_nodes))
psi_matrix_i = np.zeros((time_nodes,space_nodes))



###Set initial and boundary conditions 
for row in range (0,time_nodes):
    for col in range(0,space_nodes):
        ###Initial Conditions
        eta_matrix_i[0][col] = 1
        psi_matrix_i[0][col] = 0
        ###Boundary Conditions
        eta_matrix_i[row][space_nodes-1] = 1
        psi_matrix_i[row][0] = 0

### Create the matrices to be used in calculation, eta matrix is reversed to allow pseudo-backwards differencing               
eta_matrix_n = eta_matrix_i[::1,::-1] 
psi_matrix_n = psi_matrix_i        
eta_matrix = eta_matrix_i[::1,::-1]
psi_matrix = psi_matrix_i



### Construct the finite difference model, reversing matrices where appropriate to allow for pseudo-backwards differencing for eta calculating whilst maintaining boundary conditions
#Set number of iterations, method used here converges after just one iteration
Iteration_limit = 10
It_num = 0
for it_count in range (Iteration_limit):
    for row in range (1 ,time_nodes-1):
        for col in range (1,space_nodes-1):                        
            psi_matrix = psi_matrix[::1,::-1]
            eta_matrix_n[row][col] = (eta_matrix[row-1][col]*(Da_r/Da_w)*d_epsilon+(eta_matrix[row][col-1]*d_thetan))/((Da_r/Da_w)*d_epsilon+d_thetan+d_thetan*d_epsilon*Da_r*(1-psi_matrix[row-1][col]))
            psi_matrix = psi_matrix[::1,::-1]
        for col in range(1,space_nodes-1):
            eta_matrix_n = eta_matrix_n[::1,::-1]
            psi_matrix_n[row][col] = (psi_matrix[row-1][col]*d_epsilon+psi_matrix[row][col-1]*d_thetan+d_thetan*d_epsilon*(Da_w/alpha_ij)*eta_matrix_n[row][col])/(d_epsilon+d_thetan+d_thetan*d_epsilon*(Da_w/alpha_ij)*eta_matrix_n[row][col])
            eta_matrix_n = eta_matrix_n[::1,::-1]
                      
    It_num = It_num+1    
    print It_num
###Gauss-Seidl iteration method
    if np.allclose(eta_matrix,eta_matrix_n, rtol=1e-15) and np.allclose(psi_matrix, psi_matrix_n, rtol = 1e-15):
        print "Convergence"
        break

    eta_matrix = eta_matrix_n
    psi_matrix = psi_matrix_n
###Switches the matrix back from being reversed
eta_matrix = eta_matrix[::1,::-1]
###OUtputs the data in text files for further anlysis
np.savetxt('psimatrix.txt',psi_matrix)     
np.savetxt('etamatrix.txt',eta_matrix)      


