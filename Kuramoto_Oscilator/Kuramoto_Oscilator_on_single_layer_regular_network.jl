#########################################################################
# Time Dynamics Simulation of Coupled Kuramoto Oscillators on a Single layer Regular Network
# Devloped at Complex Systems Lab, IIT Indore India
# October, 31 2019. 
# Ref:
#Contact: sarikajalan9@gmail.com
#########################################################################
####Coupled Kuramoto Oscillators##############
# dui/dt = omegai − epss/sum (aij) * Sum [ Gij * sin{ uj - ui + alpha} ]
## Description of variables
#Input:
#             Nl           : Size of the regular Network
#             kk           : Node degree or average degree of the Network
#             r           : Coupling radius = k/2N           
#             G           : NxN Adjacency Matrix
#             epss        : overall coupling strength      
#             alpha       : phase lag parameter
#             omega       : Natural Frequency
#              u          : phase (state) variable represented as U[1:N]

#Output: 
#           sol     --- N x time matrix of the phase [1:N] time data
#                       # remember to exclude initial transient time when extracting thw full time-phase matrix

#######################Required Libraries########################################
using DifferentialEquations
using Random
using PyPlot
using Printf
using Statistics
using DataFrames
using IterableTables
using LinearAlgebra
using DelimitedFiles
using CSV

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
global Nl=250                                    #Network Size
global r=0.32                                    #Coupling Radius
global epss=0.1                                  #Overall Coupling Strength
global omega=0.01*ones(Nl)                       #Constant natural frequency for all oscilators
global alpha=1.47                                # Phase lag parameter for chimera state

######Initial Phases#####################################
#chosen such that they are randomly distributed with an gausian envelop
rng = MersenneTwister(58959);
theta0=[6*(rand(rng)-0.5)*exp(-30*(x/float(Nl)-0.5)*(x/float(Nl)-0.5))  for x in range(1,stop=Nl)]
writedlm("generated_initial.txt", theta0)
#For value references: DANIEL M. ABRAMS and STEVEN H. STROGATZ (2006) Int. Jour. of Bifurc. and Chaos, 16(1) 21–37.
#-----------------------------------------------------------------------------------------------
#Network Building#####################################################
kk = floor(Int64,r*Nl*2)
A=regular(Nl,kk)
######################################################################
# Differential Equation solver Parameters (Non-Stiff Differential Equation)
niter=31000                 # Total Iteration
novp=2000                   # Ovserved Time period after transient
dt = 0.01                   # initial time step size
ti = 0.0; tf = niter*dt
tr = (niter-novp)*dt
tspan = (ti, tf)             # Time interval
##############################################Kuramoto (KOSC) Osc Eqn.######################################
function kuramoto(du,u,p,t)
    global omega, A, epss, r,Nl
    for i = 1:Nl
        du[i] = omega[i]+ (epss/float(r*Nl*2)) * (sum(j->(A[i,j]* sin(u[j]-u[i]+alpha)),1:Nl))
    end
end
###########################Order parameter Calculation####################################################
function order_para(y)
    n = floor(Int64,length(y))
    real_sum = 0.0; imag_sum = 0.0
    for l = 1:n
        real_sum += cos(y[l])
        imag_sum += sin(y[l])
    end
    real_sum = real_sum/n
    imag_sum = imag_sum/n
    r = sqrt((real_sum)^2 + (imag_sum)^2)
    return r
end

###############################Solving the stiff Differential Equation with Autovern method########################################
#system integration
prob = ODEProblem(kuramoto,theta0,tspan)
sol = solve(prob,Tsit5(), reltol=1e-6, saveat=dt) 

##########################Extracting for Nl x observed time phase-data matrix (without fmod operation by 2pi)
df = DataFrame(sol[:,niter-novp:niter]);
#CSV.write("time_data.csv",df,writeheader=false)   ## Use CAUTION....fmod operation has not been done 
###################Extracting for a single time(final_state) Nl x 1 time phase-data matrix
writedlm("time_data.txt",[mod2pi(i)-pi for i in sol[:,niter]])







##########################################Experimental Use caution##############################
# #####################For Multiple alpha Values#####################################
# for Kinc = 0:0.05:2.0
#     global alpha = Kinc;

# #system integration
# prob = ODEProblem(kuramoto,theta0,tspan)
# sol = solve(prob,Tsit5(), reltol=1e-6, saveat=dt) 

# ##########################Writing the time phase-data matrix############################################
# df = DataFrame(sol[:,niter-novp:niter]);
# fname = @sprintf("out_alpha_%1.3f.csv",alpha)
# CSV.write(fname,df)

# #################################Plotting the final State################################################
# clf()
# x= [mod2pi(i)-pi for i in sol[:,niter]]
# ylim([-pi,pi])
# plot(x,"go",markersize=2)
# fname2 = @sprintf("final_state_%1.3f.png",alpha)
# savefig(fname2)

# ################################Calculating The Frequency#############################
# T=size(sol)[2]
# N=size(sol)[1]

# fovp=30000
# x = zeros(0)
# for n in range(1,stop=N)
#     sigp=[mod2pi(sol[n,t])-pi for t in range(T-fovp-1,stop=T-2)]
#     upcross=findall((sigp[1:end-1] .<=0) .& (sigp[2:end] .>0))
#     ff=2*pi*size(upcross)[1]/((upcross[end] - upcross[1])*0.01)
#     append!( x, ff )
# end

# #################################Ploting the Frequency################################
# clf()
# plot(x,"mo",markersize=2)
# fname3 = @sprintf("freq_%1.3f.png",alpha)
# savefig(fname3)
    
# end
# #############################################################################################
# y=[mod2pi(x) for x in sol[niter-novp]]
# println(order_para(y[1:100]))
# println(order_para(y[100:200]))



