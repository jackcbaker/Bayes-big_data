#### Implementation of the SGHMC algorithm for a multivariate t distribution 
#### with unknown location parameter. 
#### For more information see http://lancs.ac.uk/~bakerj1/pdfs/big_data/big_data.pdf

### Libraries
using StatsBase
using Distributions


### Main function
function hmc( t_init::Array{Float64,1}, max_iter::Int64, dataset::Array{Float64,2}, sample_size::Int64, n_obs::Int64, dimension::Int64, alpha::Float64, eta::Float64, L::Int64 )
    # t_init - initial value for \theta
    # max_iter - number of iterations
    # dataset - dataset
    # sample_size - the subsample size n
    # n_obs - number of observations
    # dimension - dimension of the problem
    # alpha - Momentum rate
    # eta - Learning rate
    # L - Trajectory length
    
    # Allocate memory
    stepsize( iter ) = step_const
    sample_out = zeros( Float64, (max_iter, dimension) )
    subsample_arr = zeros( Float64, (sample_size, dimension) )

    theta = t_init
    for ( i in 1:max_iter )
        # Print progress
        if ( i % 1000 == 0 )
            print( "$i " )
        end
        # Generate new momentum vars
        r = rand( Normal(0,1), dimension )
        # Reparameterise
        v = sqrt(eta).*r
        # Generate subsample
        subsample_arr = dataset[ sample( 1:n_obs, sample_size, replace=false ), : ]
        (v,theta) = sim_sghmc( v, theta, alpha, eta, L, subsample_arr, dimension )
        sample_out[i,:] = theta
    end
    print("\n")

    return( sample_out )
end

## Generate SGHMC dynamics
function sim_sghmc(  v::Array{Float64,1},        # Momentum
                    theta::Array{Float64,1},    # Position
                    alpha::Float64,             # Momentum param
                    eta::Float64,               # Learning param
                    L::Int64,                   # Trajectory
                    x::Array{Float64,2},        # Data
                    d::Int64                    # Dimension
                 )
    
    ## Helper functions
    function ll_diff( x, t, d )
        out = ( (3+d)/15 )*(x - t)
        out /= ( 1 + (1/15)*( dot(x-t,x-t) ) )
        return( out )
    end

    ## Constants
    batch_size = size(x)[1]
    noise_param = alpha*eta/4

    ## Generate SGHMC dynamics
    for ( i = 1:L )
        # Update theta
        theta += v
        # Calculate log lik
        llik = 0
        for ( j in 1:batch_size )
            llik += ll_diff ( vec( x[j,:] ), vec(theta), d )
        end
        # Momentum step 
        v += eta*llik - alpha*v  + rand( Normal(0,noise_param), d )
    end

    return ( v, theta )
end
