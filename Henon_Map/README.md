# Coupled Henon Map

    Xi (t+1) = henn(Xi(t)) +  eps/sum (aij) * Sum [ Gij * {henn(Xj(t)) - henn(Xi(t))} ]
    Yi(t+1) = beta * Xi(t)
    henn(x,y)=1- alpha *(x^2) + y
  
Code: Henon_Map_Dynamics_on_Single_Layer_Regular_Network.cpp
------------------------------------------------------------------------------------------------
Description of variables

Input:

            henn        : Hennon map function 
            alpha         : Bifurcation Parameter
             Nl           : Size of the regular Network
             Kl           : Node degree or average degree of the Network       
             G           : NxN Adjacency Matrix
             eps        : overall coupling strength     
             beta       : Systems Paramter
              x          :  state variable  represented as L[1:N]
              y          :  state variable  represented as L[N:2N]

              Same special initial condition (supplied) is used in x and y state variable

Output: 

             ---  time x 2N matrix of the x state variable [1:N] and y state variable [N:2N] data
 
-----------------------------------------------------------------------------------------------------------------
# Fig 4a : Incoherent state for Henon Map

Code : Henon_Map_Dynamics_on_Single_Layer_Regular_Network.cpp

Data File: Fig_4a.dat

System State: Fig_4a.png

Parameters:

alpha          : 1.4
Beta       : 0.3
             	Nl              : 100
             	Kl              : 64     
             	G              : NxN Adjacency Matrix of Regular Network with Nl,Kl
             	eps           : 0.1 
             	              
-----------------------------------------------------------------------------------------------------------------
# Fig 4b : Chimera for Henon Map

Code : Henon_Map_Dynamics_on_Single_Layer_Regular_Network.cpp

Data File: Fig_4b.dat

System State: Fig_4b.png

Parameters:

alpha          : 1.4
Beta       : 0.3
             	Nl              : 100
             	Kl              : 64     
             	G              : NxN Adjacency Matrix of Regular Network with Nl,Kl
             	eps           : 0.26

-----------------------------------------------------------------------------------------------------------------
# Fig 4c : Chimera state for Henon Map

Code : Henon_Map_Dynamics_on_Single_Layer_Regular_Network.cpp

Data File: Fig_4c.dat

System State: Fig_4c.png

Parameters:

alpha          : 1.4
Beta       : 0.3
             	Nl              : 100
             	Kl              : 64     
             	G              : NxN Adjacency Matrix of Regular Network with Nl,Kl
             	eps           : 0.3

-----------------------------------------------------------------------------------------------------------------
# Fig 4d : Coherent state for Henon Map

Code : Henon_Map_Dynamics_on_Single_Layer_Regular_Network.cpp

Data File: Fig_4d.dat

System State: Fig_4d.png

Parameters:

alpha          : 1.4
Beta       : 0.3
             	Nl              : 100
             	Kl              : 64     
             	G              : NxN Adjacency Matrix of Regular Network with Nl,Kl
             	eps           : 0.7 
             	              

             	              
 
             	              
