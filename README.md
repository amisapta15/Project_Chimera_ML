# Identification of Chimera using Machine Learning

Ref: arXiv:2001.08985 [nlin.AO]

## Model:Kuramoto Oscillator (Eq.4)

Code File: Kuramoto_Oscilator_on_single_layer_regular_network.jl

To execute: "julia filename.jl"; 
Note the dependencies.
  
Input: alpha, the phase lag parameter ('α' in Eq.3), taken from the range [0 < α < π]

Output: Time Series of phases of all *N* nodes in [*N* x 1] vector for a particular time-step in the steady state. Note that, in the code, all phase values are wrapped (mod 2π) in the [-π,π] envelop. 

In the following, we write specific values of α for generateing figures of the manuscript:

* **Fig 1a** depicting *Incoherent state* for kuramoto Oscilators for *α = 0.2* (Data File: Fig_1a.dat). Different [*N*x1] vectors, used for the traininng (table 1) and testing (table 2) ML algorithims, are generated by using different initial conditions, by considering states at different time-steps and by changing α values in the neibourhood of the value mentioned above.

* **Fig 1b** and **Fig 1c** depicting *Chimera state* for kuramoto Oscilators for simulated for *α = 1.4* (Data File: Fig_1b.dat), *α = 1.5* (Data File: Fig_1c.dat), respectively. Different [*N*x1] vectors, used for the traininng (table 1) and testing (table 2) ML algorithims are generated by using different initial conditions, by considering states at different time-steps and by changing α values in the neibourhood of the value mentioned above.

* **Fig 1d** depicts Coherent state for kuramoto Oscilators for α = 1.65 (Data File: Fig_1d.dat). Different [*N*x1] vectors, used for the traininng (table 1) and testing (table 2) ML algorithims are generated by considering states at different time-steps, and changing α values in the neibourhood of the value mentioned above.

## Model:FitzHugh-Nagumo (FHN) Oscillators with Rotational Coupling (Eq.5)

Code File: FitzHugh–Nagumo_Oscilator_on_single_layer_ring_network.jl

To execute: "julia filename.jl";  
Note the dependencies.

Input: phii, the rotation angle of the coupling matrix ('φ' in Eq.5), taken from the range [π/2-0.35 < φ < π/2 +0.35]

Output: Time Series of all N nodes in [*N* x 1] vector for excitatory varriable u, at a particular time-step in the steady state. Note that each node has two variables u,v, only u variable has been selected to display the figures, resulting in N time series.

In the following, we write specific values of φ for generating figures in the manuscript:

* **Fig 2a** depicting *Incoherent state* for FHN Oscilators for *φ = π/2 + 0.35* (Data File: Fig_2a.dat). 

* **Fig 2b** and **Fig 2c** depicts *Chimera state* for FHN Oscilators for *φ = π/2 - 0.1* (Data File: Fig_2b.dat), *φ = π/2 - 0.15* (Data File: Fig_2c.dat), respectively. 

* **Fig 2d** depicts *Coherent state* for FHN Oscilators for *φ = π/2 - 0.35* (Data File: Fig_2d.dat). 

Datasets are generated using the above mentioned code file in the same fashion as mentioned in Fig 1.


## Model: Coupled Logistic Map (Eq.7)

Code File: Logistic_Map_Dynamics_on_Single_Layer_Regular_Network.cpp

To compile: "g++ filename.cpp -o output"; 
To execute: "./output"

Input: eps, the overall coupling strength ('ϵ' in Eq.6) taken from the range [0 < ϵ < 1]

Output: Time Series of the state variable for all N nodes in [*N* x 1] vector at a particular time-step in the steady state.

In the following, we write specific values of ϵ for generating figures in the manuscript:

* **Fig 3a** depicting *Incoherent state* for coupled Logistic map for *ϵ = 0.1* (Data File: Fig_3a.dat). 

* **Fig 3b** and **Fig 3c** depicts *Chimera state* for coupled Logistic map for *ϵ = 0.34* (Data File: Fig_3b.dat), *ϵ = 0.37* (Data File: Fig_3c.dat), respectively. 

* **Fig 3d** depicts *Coherent state* for coupled Logistic map for *ϵ = 0.7* (Data File: Fig_3d.dat). 

Datasets are generated using the above mentioned code file in the same fashion as mentioned in Fig 1.

## Model: Coupled Henon Map (Eq.8)

Code File: Henon_Map_Dynamics_on_Single_Layer_Regular_Network.cpp

To compile: "g++ filename.cpp -o output"; 
To execute: "./output"

Input: eps, the overall coupling strength ( 'ϵ' in Eq.7) taken from the range [0 < ϵ < 1]

Output: Time Series of the state variable for all N nodes in [*N* x 1] vector at a particular time-step in the steady state.

In the following, we write specific values of ϵ for generating figures in the manuscript:

* **Fig 4a** depicting *Incoherent state* for coupled Henon map for *ϵ = 0.1* (Data File: Fig_4a.dat). 

* **Fig 4b** and **Fig 4c** depicts *Chimera state* for coupled Henon map for *ϵ = 0.26* (Data File: Fig_4b.dat), *ϵ = 0.3* (Data File: Fig_4c.dat), respectively. 

* **Fig 4d** depicts *Coherent state* for coupled Henon map for *ϵ = 0.7* (Data File: Fig_4d.dat). 


Datasets are generated using the above mentioned code file in the same fashion as mentioned in Fig 1.

## Machine Learning Algorithms Matlab Code:

This code is based on Zhang Le's implementation available from 
(https://github.com/P-N-Suganthan/CODES/blob/master/2015-TCyb-oblique-RF.rar

Paper:  L. Zhang, P. N. Suganthan, “Oblique Decision Tree Ensemble via Multisurface Proximal Support Vector Machine,” IEEE Trans on Cybernetics, DoI: 10.1109/TCYB.2014.2366468 , Vol. 45, No. 10, pp. 2165-2176, Oct 2015. 

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.)

This code is not optimized for efficiency. It has been tested on Windows 10 with MATLAB R2017a.

Please cite the following papers if you are using this code.

Reference: 
1. Ganaie, M. A., Saptarshi Ghosh, Naveen Mendola, M. Tanveer, and Sarika Jalan. "Identification of Chimera using Machine Learning." arXiv preprint arXiv:2001.08985 (2020).
2. Breiman, Leo. "Random forests." Machine learning 45.1 (2001): 5-32.
3. Zhang, Le, and Ponnuthurai N. Suganthan. "Oblique Decision Tree Ensemble via Multisurface Proximal Support Vector Machine."  IEEE Transactions on Cybernetics, Year: 2015, Volume: PP, Issue: 99(2014).
4. Zhang, Yongshan, Jia Wu, Zhihua Cai, Bo Du, and S. Yu Philip. "An unsupervised parameter learning model for RVFL neural network." Neural Networks 112 (2019): 85-97.
5. Ganaie, M. A., M. Tanveer, and P. N. Suganthan. "Oblique Decision Tree Ensemble via Twin Bounded SVM." Expert Systems with Applications 143 (2020): 113072.




%%%%%%%% How to run the code %%%%%%%
1) Generate the data corresponding to different models using the above given procedure.
2) Data_path: is the path of the dataset used in the main.m file
3) Run main.m file
4) Labels of data corresponding to four models is stored in the the variable ALL_Labels.mat
   - Note: The meaning of labels is given as follows:
     * Chimera: 1
     * Coherent: 2
     * InCoherent: 3

For more details please refer to the paper.

If you have problems about this software, please contact: phd1901141006@iiti.ac.in

04/05/2020
