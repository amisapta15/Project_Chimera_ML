# Coupled Logistic Map

    Li (t+1) = logs(Li(t)) +  eps/sum (aij) * Sum [ Gij * {logs(Lj(t)) - logs(Li(t))} ]
    logs(x)=mu * x * (1-x)
  
Code: Logistic_Map_Dynamics_on_Single_Layer_Regular_Network.cpp
------------------------------------------------------------------------------------------------
Description of variables

Input:

            logs        : Logistic map function
             mu         : Bifurcation Parameter
             Nl           : Size of the regular Network
             Kl           : Node degree or average degree of the Network       
             G           : NxN Adjacency Matrix
             eps        : overall coupling strength      
             L          :  state variable 
Output: 

                ---  time x N matrix of the state variable data 

