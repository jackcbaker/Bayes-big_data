#### Methods to recombine MCMC which has been parallelised. 
#### For more details see http://lancs.ac.uk/~bakerj1/pdfs/big_data/big_data.pdf.

library(coda)
library(reshape)
library(MASS)
library(parallelMCMCcombine)
library(mnormt)
library(ks)

recomb_scott = function( s, batch_container ) {
    # Parametric recombination algorithm developed by Scott et al. 2013
    # s - integer, number of batches.
    # batch_container - array of approximate sample from each batch of dimension (MCMC sample size, dimension, s)


    ## Constants
    p0 = ncol( batch_container[,,1] )      # Dimension of theta
    sample_size = nrow( batch_container[,,1] )

    temp_container = array( rep( NA, s*p0*resamp_size ), c(p0,resamp_size,s) )
    for ( i in 1:s ) {
        temp_container[,,i] = t( batch_container[,,i] )
    }

    sample_consens = consensusMCindep( temp_container )

    return( t( sample_consens ) )
}


recomb_para = function( s, batch_container ) {
    # Parametric recombination algorithm developed by Neiswanger et al. 2014
    # s - integer, number of batches.
    # batch_container - array of approximate sample from each batch of dimension (MCMC sample size, dimension, s)

    ## Constants
    p0 = ncol( batch_container[,,1] )      # Dimension of theta
    sample_size = nrow( batch_container[,,1] )

    ## Variables
    mu_M = rep( 0, p0 )        # Normal approximation mean
    sigma_M = matrix( rep(0, p0*p0), ncol=p0 )     # Normal approximation variance

    # Calculate sigma_M and sigma_i (correcting for correlation in MCMC chain)
    inv_sigma_batch = array( rep( NA, p0*p0*s ), c( p0,p0,s ) )
    sigma_i = matrix( rep( NA, p0*p0 ), ncol=p0 )
    for ( i in 1:s ) {
        sigma_i = var( batch_container[,,i] )
        inv_sigma_batch[,,i] = solve( sigma_i )
        sigma_M = sigma_M + inv_sigma_batch[,,i]
    }
    sigma_M = solve( sigma_M )

    # Calculate mu_M
    for ( i in 1:s ) {
        mu_i = colMeans( batch_container[,,i] )
        mu_M = mu_M + inv_sigma_batch[,,i] %*% mu_i
    }
    mu_M = sigma_M %*% mu_M

    return( list( "mean"= as.vector(mu_M), "var"=sigma_M ) )

}


recomb_npara = function( s, batch_container, bw_init, anneal=TRUE ) {
    # Nonparametric recombination algorithm developed by Neiswanger et al. 2014
    # s - integer, number of batches.
    # batch_container - array of approximate sample from each batch of dimension (MCMC sample size, dimension, s)
    # bw_init - initial bandwidth to use in kde. Must be a d x d matrix, where d is the dimension of the problem. Non-diagonal entries are ignored.
    # anneal - logical, should the algorithm be annealed.

    ## Constants
    p0 = ncol( batch_container[,,1] )      # Dimension of theta
    sample_size = nrow( batch_container[,,1] )
    # Ensure bw_init is diagonal
    bw_init = diag( diag( bw_init ) )

    ## Helper Functions
    calc_weights = function( t_lab, batch_container, s, bw ) {
        theta_vals = matrix( rep( NA, s*p0 ), ncol=p0 )
        for ( i in 1:s ) {
            theta_vals[i,] = batch_container[t_lab[i],,i]
        }
        theta_mean = colMeans( theta_vals )
        # Store log weights for numerical stability
        weight_val = sum( log( dmnorm( theta_vals, mean = theta_mean, var = bw^2 ) ) ) 

        return( weight_val )
    }

    calc_mean = function( t_lab, batch_container, s, bw ) {
        theta_mean = 0
        for ( i in 1:s ) {
            theta_mean = theta_mean + batch_container[t_lab[i],,i]/s
        }

        return( theta_mean )
    }
    
    ## Initialise
    sample_container = matrix( rep( NA, p0*sample_size ), ncol=p0 )
    T_current = sample( sample_size, s, replace=TRUE )
    bw = bw_init
    w_t = calc_weights( T_current, batch_container, s, bw )
    for ( j in 1:sample_size ) {
        if ( anneal == TRUE ) {
            h = j^(-1/(4+p0))
            bw = h*bw_init
        }
        else {
            bw = bw_init
        }
        for ( i in 1:s ) {
            C_current = T_current
            C_current[i] = sample( sample_size, 1 )
            u = runif(1)
            w_c = calc_weights( C_current, batch_container, s, bw  )
            if ( log(u) < w_c - w_t ) {
                T_current = C_current
                w_t = w_c
            }
        }
        mean_current = calc_mean( T_current, batch_container, s, bw )
        sample_container[j,] = rmnorm( 1, mean=mean_current, var=bw^2/s )
    }

    return( sample_container )
}


recomb_semipara = function( s, batch_container, bw_init, anneal = TRUE ) {
    # Semiparametric recombination algorithm developed by Neiswanger et al. 2014
    # s - integer, number of batches.
    # batch_container - array of approximate sample from each batch of dimension (MCMC sample size, dimension, s)
    # bw_init - initial bandwidth to use in kde. Must be a d x d matrix, where d is the dimension of the problem. Non-diagonal entries are ignored.
    # anneal - logical, should the algorithm be annealed.

    ## Constants
    p0 = ncol( batch_container[,,1] )      # Dimension of theta
    sample_size = nrow( batch_container[,,1] )
    # Ensure bw_init is diagonal
    bw_init = diag( diag( bw_init ) )
    samp_var = array( rep( NA, p0*p0*s ), c(p0,p0,s) )
    samp_mean = matrix( rep( NA, p0*s ), ncol=p0 )
    for ( i in 1:s ) {
        samp_var[,,i] = var( batch_container[,,i] )
        samp_mean[i,] = colMeans( batch_container[,,i] )
    }
    var_tot_inv = 0
    mean_tot = 0
    for ( i in 1:s ) {
        current_var_inv = solve( samp_var[,,i] )
        var_tot_inv = var_tot_inv + current_var_inv
        mean_tot = mean_tot + current_var_inv %*% samp_mean[i,]
    }
    mean_tot = solve( var_tot_inv ) %*% mean_tot

    ## Helper Functions
    calc_weights = function( t_lab, batch_container, s, bw ) {
        theta_vals = matrix( rep( NA, s*p0 ), ncol=p0 )
        for ( i in 1:s ) {
            current_index = t_lab[i]
            theta_vals[i,] = batch_container[current_index,,i]
        }
        theta_mean = colMeans( theta_vals )
        # Store logs for numerical stability
        weight_val = sum( log( dmnorm( theta_vals, mean = theta_mean, var = bw^2 ) ) )

        return( weight_val )
    }

    calc_mean = function( t_lab, batch_container, s, bw ) {
        theta_mean = 0
        for ( i in 1:s ) {
            current_index = t_lab[i]
            theta_mean = theta_mean + batch_container[current_index,,i]/s
        }

        return( theta_mean )
    }
    
    ## Initialise
    sample_container = matrix( rep( NA, p0*sample_size ), ncol=p0 )
    T_current = sample( sample_size, s, replace=TRUE )
    C_current = T_current
    bw = diag( rep( 1, 2 ) )
    w_t = calc_weights( C_current, batch_container, s, bw )
    # Begin MCMC
    for ( j in 1:sample_size ) {
        if ( anneal == TRUE ) {
            h = j^(-1/(4+p0))
            bw = h*bw_init
        }
        else {
            bw = bw_init
        }
        for ( i in 1:s ) {
            C_current = T_current
            C_current[i] = sample( sample_size, 1 )
            u = runif(1)
            w_c = calc_weights( C_current, batch_container, s, bw  )
            if ( log(u) < w_c - w_t ) {
                T_current = C_current
                w_t = w_c
            }
        }
        kern_var = s * solve( bw )
        var_current = solve( kern_var + var_tot_inv )
        mean_current = calc_mean( T_current, batch_container, s, bw )
        mean_current = var_current %*% ( kern_var %*% mean_current + var_tot_inv %*% mean_tot )
        sample_container[j,] = rmnorm( 1, mean=mean_current, var=var_current )
    }

    return( sample_container )
}


recomb_semipara_full = function( s, batch_container, bw, anneal = TRUE ) {
    # Semiparametric recombination algorithm developed by Neiswanger et al. 2014 
    # (algorithm NeisSemiparaFull in the report http://lancs.ac.uk/~bakerj1/pdfs/big_data/big_data.pdf)
    # s - integer, number of batches.
    # batch_container - array of approximate sample from each batch of dimension (MCMC sample size, dimension, s)
    # bw_init - initial bandwidth to use in kde. Must be a d x d matrix, where d is the dimension of the problem. Non-diagonal entries are ignored.
    # anneal - logical, should the algorithm be annealed.

    bw = diag( bw ) 
    temp_container = array( rep( NA, s*p0*resamp_size ), c(p0,resamp_size,s) )
    for ( i in 1:s ) {
        temp_container[,,i] = t( batch_container[,,i] )
    }

    sample_out = semiparamDPE( temp_container, bandw=bw, anneal = anneal )

    return( t( sample_out ) )
}


recomb_npara_grid = function( s, batch_container, bw, gridsize=30 ) {
    # Nonparametric recombination algorithm developed by Neiswanger et al. 2014
    # Method is implemented exactly on a grid, rather than using MCMC.
    # Referred to as the grid methd in http://lancs.ac.uk/~bakerj1/pdfs/big_data/big_data.pdf
    # s - integer, number of batches.
    # batch_container - array of approximate sample from each batch of dimension (MCMC sample size, dimension, s)
    # bw_init - initial bandwidth to use in kde. Must be a d x d matrix, where d is the dimension of the problem. Non-diagonal entries are ignored.
    # anneal - logical, should the algorithm be annealed.

    ## Constants
    bw = diag( diag( bw ) )
    p0 = ncol( batch_container[,,1] )      # Dimension of theta
    sample_size = nrow( batch_container[,,1] )


    ## Generate limits for kde grid
    kde_limits = mat.or.vec(2,p0)
    extremes = mat.or.vec(2,p0)
    for ( i in 1:s ) {
        extremes[1,] = apply( batch_container[,,i], 2, min )
        extremes[2,] = apply( batch_container[,,i], 2, max )
        
        for ( dim_index in 1:p0 ) {
            if ( extremes[1,dim_index] < kde_limits[1,dim_index] ) {
                kde_limits[1,dim_index] = extremes[1,dim_index]
            }
            if ( extremes[2,dim_index] > kde_limits[2,dim_index] ) {
                kde_limits[2,dim_index] = extremes[2,dim_index]
            }
        }
    }


    ## Approximate using kde for each batch and combine
    dens_index = p0 + 1
    kde_container = matrix( rep( 1, dens_index*gridsize^p0 ), ncol=dens_index )
    temp_container = list()
    for ( dim_index in 1:p0 ) {
        temp_container[[paste("theta",dim_index,sep="")]] = seq( from=kde_limits[1,dim_index], to=kde_limits[2,dim_index], length.out=gridsize )
        kde_container[,1:p0] = as.matrix( expand.grid( temp_container ) )
    }

    
    ## Run kde estimation using each batch
    for ( i in 1:s ) {
        # Estimate bandwidth for batch
        cat(paste(i," ",sep=""))
        kde_container[,dens_index] = kde_container[,dens_index] * apply( kde_container[,1:p0], 1, 
              function(x) 1/sample_size * sum( dmnorm( matrix( rep( x, each=sample_size ), ncol=p0 ), 
                                                  mean=batch_container[,,i], var=bw ) ) )
        # Take out common factor and normalise
        kde_container[,dens_index] = exp( log( kde_container[,dens_index] ) 
                                - max( log( kde_container[,dens_index] ) ) )
        kde_container[,dens_index] = kde_container[,dens_index] / sum(kde_container[,dens_index])
    }
    cat("\n")

    return( kde_container )

}
