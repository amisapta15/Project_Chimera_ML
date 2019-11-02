# Coupled Logistic Map

    Li (t+1) = logs(Li(t)) +  eps/sum (aij) * Sum [ Gij * {logs(Lj(t)) - logs(Li(t))} ]
    logs(x)=mu * x * (1-x)
  
Code: Logistic_Map_Dynamics_on_Single_Layer_Regular_Network.cpp
------------------------------------------------------------------------------------------------
Description of variables

Input:

             mu         : Bifurcation Parameter
             Nl           : Size of the regular Network
             Kl           : Node degree or average degree of the Network       
             G           : NxN Adjacency Matrix
             eps        : overall coupling strength      
             L          :  state variable 
             
             A special semi-ranom initial condition is attached. 
Output: 

                ---  time x N matrix of the state variable data   
-----------------------------------------------------------------------------------------------------------------
# Fig 3a : Incoherent state for Logistic Map

Code : Logistic_Map_Dynamics_on_Single_Layer_Regular_Network.cpp

Data File: Fig_3a.dat

System State: Fig_3a.png

Parameters:

             mu         : 4.0
             Nl         : 100
             Kl         : 64       
             G          : Adjacency Matrix of Regular Network with Nl,Kl
             eps        : overall coupling strength      
-----------------------------------------------------------------------------------------------------------------
