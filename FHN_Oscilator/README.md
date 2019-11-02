# Model:FitzHugh-Nagumo (FHN) Oscillators with Rotational Coupling 

           epss * (dui/dt) = ui − (ui^3)/3 − vi + coup/sum (aij) * Sum [ Gij * {buu (uj-ui) + buv (vj-vi)} ]

           dvi/dt = ui + a + coup/sum (aij) * Sum [ Gij * {bvu (uj-ui) + bvv (vj-vi)} ]

Code: FitzHugh–Nagumo_Oscilator_on_single_layer_ring_network.jl
-----------------------------------------------------------------------------------------
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
# Fig 2a : Incoherent state for FHN Oscilators  

Code : FitzHugh–Nagumo_Oscilator_on_single_layer_ring_network.jl

Data File: Fig_2a.dat

System State: Fig_2a.png

Parameters:


             N           : 300
             
             k           : 210
             
             r           : 0.35
             
             G           : Regular Network with N,r
             
             coup        : 0.1
             
             'φ' (phii)  : pi/2 + 0.35
             
              a          :  0.5
              
              epss       : 0.05

-----------------------------------------------------------------------------------------------------------------
# Fig 2b : Chimera state for FHN Oscilators  

Code : FitzHugh–Nagumo_Oscilator_on_single_layer_ring_network.jl

Data File: Fig_2b.dat

System State: Fig_2b.png

Parameters:


             N           : 300
             
             k           : 210
             
             r           : 0.35
             
             G           : Regular Network with N,r
             
             coup        : 0.1
             
             'φ' (phii)  : pi/2 - 0.1
             
              a          :  0.5
              
              epss       : 0.05
              
-----------------------------------------------------------------------------------------------------------------
# Fig 2c : Chimera state for FHN Oscilators  

Code : FitzHugh–Nagumo_Oscilator_on_single_layer_ring_network.jl

Data File: Fig_2c.dat

System State: Fig_2c.png

Parameters:


             N           : 300
             
             k           : 210
             
             r           : 0.35
             
             G           : Regular Network with N,r
             
             coup        : 0.1
             
             'φ' (phii)  : pi/2 - 0.15
             
              a          :  0.5
              
              epss       : 0.05
              
-----------------------------------------------------------------------------------------------------------------
# Fig 2d : Coherent state for FHN Oscilators  

Code : FitzHugh–Nagumo_Oscilator_on_single_layer_ring_network.jl

Data File: Fig_2d.dat

System State: Fig_2d.png

Parameters:


             N           : 300
             
             k           : 210
             
             r           : 0.35
             
             G           : Regular Network with N,r
             
             coup        : 0.1
             
             'φ' (phii)  : pi/2 - 0.35
             
              a          :  0.5
              
              epss       : 0.05
