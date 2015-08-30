#### Implementation of the SGLD algorithm for a multivariate t distribution 
#### with unknown location parameter. 
#### For more information see http://lancs.ac.uk/~bakerj1/pdfs/big_data/big_data.pdf

### Libraries
using StatsBase
using Distributions


### Main function
function stoch_grad_desc( t_init::Array{Float64,1}, max_iter::Int64, dataset::Array{Float64,2}, sample_size::Int64, n_obs::Int64, dimension::Int64, step_const::Float64 )
    # t_init - initial value for \theta
    # max_iter - number of iterations
    # dataset - dataset
    # sample_size - the subsample size n
    # n_obs - number of observations
    # dimension - dimension of the problem
    # step_const - constant a for tuning the stepsize
    

    ## Gradient of the pdf of a multivariate t distribution with location parameter t.
    function ll_diff( x, t, d )
        out = ( (3+d)/15 )*(x - t)
        out /= ( 1 + (1/15)*( dot(x-t,x-t) ) )
        return( out )
    end

    ## Stepsize function
    stepsize( iter ) = step_const*1/n_obs*(iter+1)^(-.33)


    ## Allocate memory
    sample_out = zeros( Float64, (max_iter, dimension+1) )


    ## Main Body
    theta_current = t_init
    sample_set = zeros( Float64, (sample_size,2) )
    sample_indices = zeros( Int64, sample_size )
    for ( iter = 1:max_iter )
        # Print progress
        if ( iter % 100 == 0 )
            print("$iter ")
        end
        # Generate subsample
        sample_indices = sample( 1:n_obs, sample_size, replace=false )
        sample_set = dataset[sample_indices,:]
        # Simulate SGLD
        theta_diff = zeros( Float64, dimension )
        for ( elem = 1:sample_size )
            theta_diff += ll_diff( vec(sample_set[elem,:]), vec(theta_current), dimension )
        end
        step_current = stepsize( iter )
        theta_diff *= step_current*n_obs/(2*sample_size)
        theta_current += theta_diff + rand( Normal(0,step_current), dimension )
        # Store theta and stepsize
        sample_out[iter,:] = [theta_current, step_current]
    end
    print("\n")

    return( sample_out )
end
