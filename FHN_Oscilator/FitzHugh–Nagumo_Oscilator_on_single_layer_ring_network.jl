#########################################################################
# Time Dynamics Simulation of Coupled FitzHugh-Nagumo (FHN) Oscillators on a Single layer Regular Network
# Devloped at Complex Systems Lab, IIT Indore India
# October, 31 2019. 
# Ref:
#Contact: sarikajalan9@gmail.com
#########################################################################
####Coupled FitzHugh-Nagumo (FHN) Oscillators with Rotational Coupling matrix##############
# epss * (dui/dt) = ui − (ui^3)/3 − vi + coup/sum (aij) * Sum [ Gij * {buu (uj-ui) + buv (vj-vi)} ]
# dvi/dt = ui + a + coup/sum (aij) * Sum [ Gij * {bvu (uj-ui) + bvv (vj-vi)} ]

## Description of variables
#coupling matrix
#           [buu,buv ; bvu,bvv] = [cosφ , sinφ ; -sinφ , cosφ]
#Input:
#             N           : Size of the regular Network
#             k           : Node degree or average degree of the Network
#             r           : Coupling radius = k/2N           
#             G           : NxN Adjacency Matrix
#             coup        : overall coupling strength      
#             'φ' (phii)  : rotation angle
#              a          :  excitability threshold
#              epss       :  Systems Parameter
#              u          :  excitatory varriable  represented as U[1:N]
#              v          :  inhibitory varriable  represented as U[N:2N]

#Output: 
#           sol     --- 2N x time matrix of the excitatory [1:N] and inhibitory [N:2N] Data
#                       # remember to exclude initial transient time when extracting

#######################Required Libraries########################################
using DifferentialEquations
using Random
using DataFrames
using Statistics
using IterableTables
using CSV

#Global varriable declaration
global coup,phii,aa,A,epss

################Regular Network G generating function###############################
function regular(n::Integer,k::Integer)
    !iseven(n*k) && throw(ArgumentError("n*k must be even!"))
    !(0<=k<n) &&  throw(ArgumentError("The 0<= k<n must be satisfied!!"))

    G=zeros((n,n))
    for i=1:n-1
        for j=i+1:n
            if(abs(i-j)<=k/2)
                G[i,j]=G[j,i]=1
            else
                if((n-abs(i-j))<=abs(k/2))
                    G[i,j]=G[j,i]=1
                end
            end
        end
    end

    return G

end
########################################################################
#Input Parameters
N=300          
r=0.35         
coup=0.1       
phii=(pi/2.0)-0.1 #Parameter for Chimera
aa=0.5            
epss=0.05

########Initial Condition
#chosen such that they are randomly distributed on a circle with u^2 + v^2 = 4
rng = MersenneTwister(599);
theta0=2pi*rand(rng,N)
uv0=zeros((2*N))
for s= 1:N
            uv0[s]= 2 * cos(theta0[s])
            uv0[s+N]= 2 * sin(theta0[s])
     end
#For value references: L Schülen, S Ghosh, et al.(2019) Chaos, Solitons & Fractals 128, 290-296
#########################################################################################
# Differential Equation solver Parameters (Stiff Differential Equation)
niter=21000     # Total Iteration
novp=1000       # Ovserved Time period after transient
dt = 0.01       # initial time step size
ti = 0.0; tf = niter*dt
tr = (niter-novp)*dt
tspan = (ti, tf)    # Time interval


#Network Building
kk = floor(Int64,r*N*2)
A=regular(N,kk)

#To record a picture of the initial condition
# uu=uv0[1:N]
# vv=uv0[N+1:2N]
# scatter(uu,vv)
# savefig("initial.png")

#FitzHugh-Nagumo (FHN) Osc Eqn.
function fhn(du,u,p,t)
    global coup,phii,aa,epss,A,N

    for i in 1:N
        du[i]=  (   u[i]  - (((u[i])^3)/ 3.0)  - u[i+N]  +  (coup  /(float(r*N*2)))   * sum(j-> (   A[i,j] * (    cos(phii) *   (u[j]-u[i])  +  sin(phii) * (u[j+N] - u[i+N])  )) ,        1:N ) ) /  float(epss)
        du[i+N]=  u[i]  + aa  + (coup /(float(r*N*2)))   * sum(k-> (     A[i,k] * (    - sin(phii) * (u[k] - u[i])  +  cos(phii) * (u[k+N] - u[i+N]) )) , 1:N )
    end
end

########Solving the stiff Differential Equation with Autovern method
prob = ODEProblem(fhn, uv0, tspan)
sol = solve(prob,AutoVern9(KenCarp4()),reltol=1e-6)

df = DataFrame(sol)
########################Extracting the whole 2N x time data
#CSV.write("time_data1.csv",df)
##########################Extracting for N (only excitatory varriable; u) x observed time data matrix
#CSV.write("time_data2.csv",df[1:N,size(df)[2]-novp:size(df)[2]],writeheader=false)
###################Extracting for a single time(final_state) N(only excitatory varriable; u) x 1 time data matrix
CSV.write("time_data.csv",df[1:N,[size(df)[2]]],writeheader=false)
## Extracing for both excitatory and inhibitory varriable for single time  2N x 1 time data matrix
#writedlm( "FileName.csv",  df[:,names(df)[end]], '\n')



