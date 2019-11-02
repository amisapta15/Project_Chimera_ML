Code: FitzHugh–Nagumo_Oscilator_on_single_layer_ring_network.jl

Model:FitzHugh-Nagumo (FHN) Oscillators with Rotational Coupling 

epss * (dui/dt) = ui − (ui^3)/3 − vi + coup/sum (aij) * Sum [ Gij * {buu (uj-ui) + buv (vj-vi)} ]

dvi/dt = ui + a + coup/sum (aij) * Sum [ Gij * {bvu (uj-ui) + bvv (vj-vi)} ]

----------------------------------------------------------------------------------------------------------------------
Description of variables

Coupling matrix : 

           [buu,buv ; bvu,bvv] = [cosφ , sinφ ; -sinφ , cosφ]
           
Input:

             N           : Size of the regular Network
             
             k           : Node degree or average degree of the Network
             
             r           : Coupling radius = k/2N 
             
             G           : NxN Adjacency Matrix
             
             coup        : overall coupling strength    
             
             'φ' (phii)  : rotation angle
             
              a          :  excitability threshold
              
              u          :  excitatory varriable  represented as U[1:N]
              
              v          :  inhibitory varriable  represented as U[N:2N]
              
              Initial u,v chosen such that they are randomly distributed on a circle with u^2 + v^2 = 4 

Output: 

           sol     --- 2N x time matrix of the excitatory [1:N] and inhibitory [N:2N] Data
           
                       remember to exclude initial transient time when extracting
                       
-----------------------------------------------------------------------------------------------------------------
Fig 2a : Incoherent state for FHN Oscilators  
Code : FitzHugh–Nagumo_Oscilator_on_single_layer_ring_network.jl
Parameters:
N=300;          
r=0.35;         
coup=0.1       
phii=(pi/2.0)-0.1 #Parameter for Chimera
aa=0.5            
epss=0.05

             N           : 300
             
             k           : 210
             
             r           : 0.35
             
             G           : Regular Network with N,r
             
             coup        : 0.1
             
             'φ' (phii)  : rotation angle
             
              a          :  0.5
              
              u          :  excitatory varriable  represented as U[1:N]
              
              v          :  inhibitory varriable  represented as U[N:2N]
