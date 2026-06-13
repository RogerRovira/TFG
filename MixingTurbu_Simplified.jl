# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.19.3
#   kernelspec:
#     display_name: Julia 1.9.0
#     language: julia
#     name: julia-1.9
# ---

using DelimitedFiles
using Plots
using FFTW

# +
function update_turbu_and_concentration(theta, conc, A, R, S, Diff, ops)
    
    #START FILTERING THETA
    theta_tr = ops.dealias .* fft(theta);
    theta_DA = real.(ifft(theta_tr))
    
    sin2t_tr =  fft( sin.(2 * theta_DA) );
    cos2t_tr =  fft( cos.(2 * theta_DA) );
    
    #COMPUTE THE DERIVATIVES IN FOURIER SPACE AND THEN BACKTRANSFORM
    dxthet = real.(ifft( im * ops.dealias .* ops.kx .* theta_tr ))
    dythet = real.(ifft( im * ops.dealias .* ops.ky .* theta_tr ))
    
    
    D_tr =   ops.dealias .* (ops.kx2_ky2 .* sin2t_tr ./2 .- ops.kxy .* cos2t_tr ) ;
    
    #THIS TERM IS CALLED ERIKSEN STRESS.
    eric = real.(ifft( im * ops.dealias .* ops.kx .* ops.k2 .* theta_tr ) ) .*  dythet .- real.( ifft( im * ops.dealias .* ops.ky .* ops.k2 .* theta_tr  ) ) .*  dxthet;
    eric_tr = ops.dealias .* fft( eric );
    
    #SOLUTION OF THE EQUATION FOR PSI, INVOLVING NABLA^4
    psi_tr =  ( S  * D_tr .+ R/A * eric_tr .- R/(2 * A ) * ops.dealias .*  ops.k4_im .* theta_tr )./ops.k4_im;
    psi_tr = ops.dealias .* psi_tr
    
    #THEN PREPARING TERMS FOR THE EVOLUTION OF THETA
    dxpsi = real.(ifft( im * ops.dealias .* ops.kx .* psi_tr) )
    dypsi = real.(ifft( im * ops.dealias .* ops.ky .* psi_tr) )
    
    cov_der_the = dxpsi .* dythet .- dypsi  .* dxthet;
    cov_der_the_tr = ops.dealias .* fft( cov_der_the )
    
    #IN THIS STEP WE FINALLY UPDATE THETA!
    #theta_new_tr = theta_tr.+ dt.*( cov_der_the_tr  .+ dealias .* k2 .* (psi_tr ./2 .- theta_tr./A) );

    #ALTERNATIVE!!! SEMI-IMPLICIT EULER (IMEX) TO INCREASE STABILITY
    theta_new_tr = (theta_tr .+ ops.dt .* (cov_der_the_tr .+ ops.dealias .* ops.k2 .* psi_tr ./ 2)) ./ (1 .+ ops.dt .* ops.k2 ./ A);

    theta_out = real.( ifft( theta_new_tr ) ) ;
   
    #UPDATE THE CONCENTRATION FIELD
    conc_tr = ops.dealias .* fft(conc);
    #conc = real.(ifft( conc_tr ));
   
    dxconc  = real.( ifft(  im .* ops.kx .* conc_tr ));
    dyconc  = real.( ifft(  im .* ops.ky .* conc_tr ));
    #d2conc = real.( ifft( - k2 .* conc_tr));
   
    covc_temp = dypsi .* dxconc .- dxpsi .* dyconc;
    #ALTERNATIVE!!! SEMI-IMPLICIT EULER (IMEX) TO INCREASE STABILITY
    covc_tr   = ops.dealias .* fft( covc_temp )
    conc_tr   = ( conc_tr .- ops.dt .* covc_tr ) ./ ( 1 .+ ops.dt .* (Diff/A) .* ops.k2 )
    conc_out  = real.( ifft( conc_tr ) )
   
    #covc = real.(ifft( dealias .* fft( covc_temp) ));
   
    #conc_out = conc .+ dt .* ( -covc .+ Diff/A .* d2conc );
    #conc_out = conc .+ dt .* (Diff/A .* d2conc)
    
    
    return theta_out, psi_tr, conc_out;
    
end

# +
# Numerical params
global L=Int(256);    # number of elements in row/column
 
global A = 10000;     # activity
global R = 1.;          # viscosity ratio
global S = 1.;          # 1. = conractile, -1 = extensile
global Diff = 1;

# +
## SAVING PARAMS.

global T=Int(200);       # integration time
global dt=0.002;         # timestep

global n_saving = 200;    # number of frames to be saved
global n_time_steps = Int(T/dt);
global time_gap = Int(n_time_steps/n_saving);


#CHANGE THIS PATH TO SAVE IN ANOTHER FOLDER
global fold_out = "COLLAPSE_L_$L"*"_A_$A/";

#CREATING THE FOLDER, IF IT DOESN'T EXIST
try
    mkpath(fold_out);
catch   
end


#SAVE THE PARAMETERS YOU HAVE CHOSEN.
par_names = "A,  R, S, L, T, dt, time_gap \n";
num_params = 7

par_vals = Array{Float64, 1}(undef, num_params)
par_vals = [ A,  R, S, L, T, dt, time_gap ];

open(fold_out*"params.txt", "w") do file
    write(file, par_names)
    
    for k = 1:num_params
        write(file, string(par_vals[k]) )
        write(file, ", ")
    end
    
end

# +
#DEFINING THE OPERATORS YOU WILL USE IN THE FOURIER SPACE
global kx = reshape(fftfreq(L, 2π*L), (L, 1));
global ky = reshape(fftfreq(L, 2π*L), (1, L));

global kxy = @. kx .* ky ;

global k2 = @. kx^2 + ky^2;
global kx2_ky2 = @. kx^2 - ky^2;
global ky2_kx2 = @. ky^2 - kx^2;
global k4 = @. kx^4 + ky^4 .+ 2 .* kx^2 .* ky^2 ;

#REGULARIZING K4, THAT ENTERS IN AN INVERSION!
k_small = copy(k4);
k_small = fill!(k_small, 0)
k_small[1] = 0.0000000000001

global k4_im = k4 + im * k_small;


#PREPARING THE DE-ALIASING OPERATOR
kmax = maximum(kx)*2/ 3;
ind_x_alias = reshape( abs.(kx) .> kmax, L );
ind_y_alias = reshape( abs.(ky) .> kmax, L) ;

dealias = copy(kxy)
dealias = fill!(dealias, 1.)
dealias[ind_x_alias,:] .= 0.;
dealias[:, ind_y_alias] .= 0.;
# -

# BUNDLE THE FOURIER OPERATORS (AND dt) INTO A NAMEDTUPLE.
# Passing this `ops` to update_turbu_and_concentration avoids reading untyped globals
# inside the hot loop: a NamedTuple has concrete, inferred field types, so the function
# specializes and stays type-stable. Rerun this cell whenever the operators or dt change.
ops = (
    kx      = kx,
    ky      = ky,
    k2      = k2,
    kx2_ky2 = kx2_ky2,
    kxy     = kxy,
    k4_im   = k4_im,
    dealias = dealias,
    dt      = dt,
)

# +
# PREPARING SOME ARRAYS, SOME OF THEM WILL BE USED FOR SAVING THE ITERATIONS OVER TIME, 
# SOME JUST FOR STORAGE IN THE FUNCTIONS. I DON'T DEFINE THEM INSIDE FUNCTIONS NOT TO ALLOCATE
# TOO MUCH MEMORY, BUT MAYBE IN JULIA THIS IS NOT A PROBLEM. THERE ARE SOME PRO TIPS ON THE INTERNET, IF YOU 

# Snapshots are real-valued (theta, velocity, concentration), so the storage containers
# are concretely typed as Vector{Matrix{Float64}}. This removes the boxing of the old
# abstract Array{Array} and stores real (not complex) data, speeding up access and saving.
global theta_time = Vector{Matrix{Float64}}(undef, n_saving+1)
global Ux_time    = Vector{Matrix{Float64}}(undef, n_saving+1)
global Uy_time    = Vector{Matrix{Float64}}(undef, n_saving+1)
global conc_time  = Vector{Matrix{Float64}}(undef, n_saving+1)

global theta = Array{Complex{Float64}}(undef,L,L);


#VARIABLES FOR THE TIME LOOP
global theta_new = copy(theta);
global theta_tr_new = copy(theta);
global psi_tr_new = copy(theta);


global theta_out = copy(theta);

global conc = copy(theta);
global conc_new = copy(theta);

# +
#RANDOM INITIAL PERTURBATION OF THE THETA = 0 STATE

ranAmpli = 2 .* ( rand(25) .-0.5 ) .* pi./5000;

theta = fill!(theta, 0)

for i =1:L
    for j =1:L
        theta[i,j] =  ranAmpli[1] .* sin( 2π /L *i  ) .+ ranAmpli[2] .* cos( 2π /L *i  ) +
                      + ranAmpli[3] .* sin( 2π /L *j  ) .+ ranAmpli[4] .* cos( 2π /L *j )+
                      + ranAmpli[5] .* sin( 2π /L *(i+j)  ) .+ ranAmpli[6] .* cos( 2π /L *(i+j) )+ 
                      + ranAmpli[7] .* sin( 2π /L *(i-j)  ) .+ ranAmpli[8] .* cos( 2π /L *(i-j) )+
                      + ranAmpli[9] .* sin( 2π /L *2*i  ) .+ ranAmpli[10] .* cos( 2π /L *2*i  ) +
                      + ranAmpli[11] .* sin( 2π /L *2*j  ) .+ ranAmpli[12] .* cos( 2π /L *2*j )+
                      + ranAmpli[13] .* sin( 2π /L *(2*i+j)  ) .+ ranAmpli[14] .* cos( 2π /L *(2*i+j) )+ 
                      + ranAmpli[15] .* sin( 2π /L *(i+2*j)  ) .+ ranAmpli[16] .* cos( 2π /L *(i+2*j) );
    end
end

theta = real.(theta)

#imshow(theta, cmap ="twilight_shifted")
#colorbar()

# +
conc .= 0


###################
#TANH PROFILE:
###################

# interface steepness
global del = 50/L;


global r0 = 0.05*L;
global cx = 0.5*L;
global cy = 0.5*L;


for i =1:L
    for j =1:L
        conc[i,j] =  (1 - tanh( del*( sqrt( (i-cx)^2 + (j-cy)^2 ) -r0 ) ) )/2
    end
end


#imshow(real.(conc), cmap ="inferno")
#colorbar()
# -

# LOOP IN TIME, EACH TIME INVOKING THE FUNCTION TURBU STEP TO PERFORM AN EULER STEP
print( "Activity: $A, Diffusion: $Diff, time: $T, ")
for counter = 0:n_time_steps
    
    theta_new, psi_tr_new, conc_new  = update_turbu_and_concentration( theta, conc, A, R, S, Diff, ops )
    theta .= theta_new
    conc .= conc_new
    
    if mod(counter, time_gap) == 0 
        print( round(counter/n_time_steps, digits = 3) , ", ")
        
        # store real-valued snapshots into the concretely-typed Vector{Matrix{Float64}}
        theta_time[trunc(Int,counter/time_gap)+1] = real.(theta)
        conc_time[trunc(Int,counter/time_gap)+1] = real.(conc)
        Ux_time[trunc(Int,counter/time_gap)+1] = real.( ifft( im * ky .* psi_tr_new) )
        Uy_time[trunc(Int,counter/time_gap)+1] = real.( ifft( -im * kx .* psi_tr_new) )
        
        if any(isnan.(theta))
            print("NaNaNaN BATMAN");
            break
        end
    end
end
