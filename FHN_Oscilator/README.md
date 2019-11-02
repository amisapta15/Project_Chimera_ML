Code: FitzHugh–Nagumo_Oscilator_on_single_layer_ring_network.jl

Model:FitzHugh–Nagumo_Oscilator

Coupled FitzHugh-Nagumo (FHN) Oscillators with Rotational Coupling matrix##############
epss * (dui/dt) = ui − (ui^3)/3 − vi + coup/sum (aij) * Sum [ Gij * {buu (uj-ui) + buv (vj-vi)} ]
dvi/dt = ui + a + coup/sum (aij) * Sum [ Gij * {bvu (uj-ui) + bvv (vj-vi)} ]

#Description of variables
coupling matrix
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

Output: 
           sol     --- 2N x time matrix of the excitatory [1:N] and inhibitory [N:2N] Data
                       remember to exclude initial transient time when extracting


