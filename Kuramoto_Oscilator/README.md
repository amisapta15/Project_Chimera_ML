# Model:Kuramoto Oscillator (KOSC)

    dui/dt = omegai − epss/sum (aij) * Sum [ Gij * sin{ uj - ui + alpha} ]

Code: Kuramoto_Oscilator_on_single_layer_regular_network.jl

-----------------------------------------------------------------------------------------------------------------

Description of variables

Input:

                                                        Nl           : Size of the regular Network
                                                        kk           : Node degree or average degree of the Network
                                                        r           : Coupling radius = k/2N   
                                                        G           : NxN Adjacency Matrix
                                                        epss        : overall coupling strength 
                                                        alpha       : phase lag parameter
                                                        omega       : Natural Frequency
                                                        u          : phase (state) variable represented as U[1:N]

Output:

                                                    sol     --- N x time matrix of the phase [1:N] time data
                                                    # remember to exclude initial transient time when extracting the full time-phase matrix

----------------------------------------------------------------------------------------------------------------
# Fig 1a : Incoherent state for Kuramoto Oscilators  

Code : Kuramoto_Oscilator_on_single_layer_regular_network.jl

Data File: Fig_1a.dat

System State: Fig_1a.png

Parameters:

             Nl           : 250
             
             kk           : 160
             
             r            : 0.32
             
             G           : Regular Network with N,r
             
             epss            : 0.1  
             
             Omega       : 0.01 for all N
             
             alpha           :1.65
             
----------------------------------------------------------------------------------------------------------------
# Fig 1b : Chimera state for Kuramoto Oscilators  

Code : Kuramoto_Oscilator_on_single_layer_regular_network.jl

Data File: Fig_1b.dat

System State: Fig_1b.png

Parameters:


             Nl           : 250
             
             kk           : 160
             
             r            : 0.32
             
             G           : Regular Network with N,r
             
             epss            : 0.1   
             
             Omega       : 0.01 for all N
             
             alpha           :1.4

----------------------------------------------------------------------------------------------------------------
# Fig 1c : Chimera state for Kuramoto Oscilators  

Code : Kuramoto_Oscilator_on_single_layer_regular_network.jl

Data File: Fig_1c.dat

System State: Fig_1c.png

Parameters:


             Nl           : 250
             
             kk           : 160
             
             r            : 0.32
             
             G           : Regular Network with N,r
             
             epss            : 0.1   
             
                omega       : 0.01 for all N
             
                alpha           :1.5
        

----------------------------------------------------------------------------------------------------------------

# Fig 1d : Coherent state for Kuramoto Oscilators  

Code : Kuramoto_Oscilator_on_single_layer_regular_network.jl

Data File: Fig_1d.dat

System State: Fig_1d.png

Parameters:


                Nl           : 250
                kk           : 160
                r            : 0.32
                G           : Regular Network with N,r
                epss            : 0.1      
                omega       : 0.01 for all N
                alpha           :0.2
