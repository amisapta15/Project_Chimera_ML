# Characterization of Chimera using Machine Learning

Ref: 

## Model:Kuramoto Oscillator (Eq.1)

Code: Kuramoto_Oscilator_on_single_layer_regular_network.jl

Input:   'α' (alpha) : phase lag parameter (Eq.1) taken from the range [0 < α < π]

Output: Time Series of phases of all N nodes in N x time (single-time or snapshot) matrix. All phase values are wrapped (mod 2π) in the [-π,π] envelop. 

In the following, we write specific values of α used to generate figures displayed in manuscript:

* **Fig 1a** depicting *Incoherent state* for kuramoto Oscilators for *α = 0.2* (Data File: Fig_1a.dat). Dataset for the traininng (table 1) and testing (table 2) ML algorithims are generated by using different initial conditions, considering states at different time-steps and changing α values in the neibourhood of the value mentioned above.

* **Fig 1b** and **Fig 1c** depicting *Chimera state* for kuramoto Oscilators for simulated for *α = 1.4* (Data File: Fig_1b.dat), *α = 1.5* (Data File: Fig_1c.dat), respectively. Dataset for the traininng (table 1) and testing (table 2) ML algorithims are generated by using different initial conditions, considering states at different time-steps and changing α values in the neibourhood of the value mentioned above.

* **Fig 1d** depicts Coherent state for kuramoto Oscilators for α = 1.65 (Data File: Fig_1d.dat). Dataset for the traininng (table 1) and testing (table 2) ML algorithims are generated by considering states at different time-steps, and changing α values in the neibourhood of the value mentioned above.

## Model:FitzHugh-Nagumo (FHN) Oscillators with Rotational Coupling (Eq.4)

Code: FitzHugh–Nagumo_Oscilator_on_single_layer_ring_network.jl

Input:   'φ' (phii)  : rotation angle of the coupling matrix (Eq.5) taken from the range [π/2-0.35 < φ < π/2 +0.35]

Output: Time Series of all N nodes in N (Only excitatory varriable u) x time (single or time window) matrix. Each node has two variables u,v, only u variable has been selected to display the figures, resulting in N time series.

In the following, we write specific values of φ used to generate figures displayed in manuscript:

* **Fig 2a** depicting *Incoherent state* for FHN Oscilators for *φ = π/2 + 0.35* (Data File: Fig_2a.dat). Dataset for the traininng (table 1) and testing (table 2) ML algorithims are generated by using different initial conditions, considering states at different time-steps and changing φ values in the neibourhood of the value mentioned above.

* **Fig 2b** and **Fig 2c** depicts *Chimera state* for FHN Oscilators for *φ = π/2 - 0.1* (Data File: Fig_2b.dat), *φ = π/2 - 0.15* (Data File: Fig_2c.dat), respectively. Dataset for the traininng (table 1) and testing (table 2) ML algorithims are generated by using different initial conditions, considering states at different time-steps and changing φ values in the neibourhood of the value mentioned above.

* **Fig 2d** depicts *Coherent state* for FHN Oscilators for *φ = π/2 - 0.35* (Data File: Fig_2d.dat). Dataset for the traininng (table 1) and testing (table 2) ML algorithims are generated by considering states at different time-steps, and changing φ values in the neibourhood of the value mentioned above.

## Model: Coupled Logistic Map (Eq.6)

Code: Logistic_Map_Dynamics_on_Single_Layer_Regular_Network.cpp

Input:   'ϵ' (eps)  : overall coupling strength taken from the range [0 < ϵ < 1]

Time Series of all N nodes in N x time (single-time or snapshot) matrix. 

In the following, we write specific values of ϵ used to generate figures displayed in manuscript:

* **Fig 3a** depicting *Incoherent state* for coupled Logistic map for *ϵ = 0.1* (Data File: Fig_3a.dat). Dataset for the traininng (table 1) and testing (table 2) ML algorithims are generated by considering states at different time-steps and changing ϵ values in the neibourhood of the value mentioned above.

* **Fig 3b** and **Fig 3c** depicts *Chimera state* for coupled Logistic map for *ϵ = 0.34* (Data File: Fig_3b.dat), *ϵ = 0.37* (Data File: Fig_3c.dat), respectively. Dataset for the traininng (table 1) and testing (table 2) ML algorithims are generated by considering states at different time-steps and changing ϵ values in the neibourhood of the value mentioned above.

* **Fig 3d** depicts *Coherent state* for coupled Logistic map for *ϵ = 0.7* (Data File: Fig_3d.dat). Dataset for the traininng (table 1) and testing (table 2) ML algorithims are generated by considering states at different time-steps, and changing ϵ values in the neibourhood of the value mentioned above.

# Coupled Henon Map

Code: Henon_Map_Dynamics_on_Single_Layer_Regular_Network.cpp

Input:   'ϵ' (eps)  : overall coupling strength taken from the range [0 < ϵ < 1]

Time Series of all N nodes in N x time (single-time or snapshot) matrix. 

In the following, we write specific values of ϵ used to generate figures displayed in manuscript:

* **Fig 4a** depicting *Incoherent state* for coupled Henon map for *ϵ = 0.1* (Data File: Fig_4a.dat). Dataset for the traininng (table 1) and testing (table 2) ML algorithims are generated by considering states at different time-steps and changing ϵ values in the neibourhood of the value mentioned above.

* **Fig 4b** and **Fig 4c** depicts *Chimera state* for coupled Henon map for *ϵ = 0.26* (Data File: Fig_4b.dat), *ϵ = 0.3* (Data File: Fig_4c.dat), respectively. Dataset for the traininng (table 1) and testing (table 2) ML algorithims are generated by considering states at different time-steps and changing ϵ values in the neibourhood of the value mentioned above.

* **Fig 4d** depicts *Coherent state* for coupled Henon map for *ϵ = 0.7* (Data File: Fig_4d.dat). Dataset for the traininng (table 1) and testing (table 2) ML algorithims are generated by considering states at different time-steps, and changing ϵ values in the neibourhood of the value mentioned above.

