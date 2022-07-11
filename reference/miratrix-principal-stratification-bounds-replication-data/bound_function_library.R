
##
## Collection of functions to calculate bounds on complier groups given data
##
## Also functions to calculate observe statistics and parameter values given complete data
##
## (C) Luke Miratrix, 2018   lmiratrix@gse.harvard.edu


library( xtable )


STRATA_NAMES = c(  "eat", "hqat", "hqc", "lqat", "lqc" )

# this messyness is to deal with empty categories.  Has to be an easier
# way, but what?
SHELL_FRAME = data.frame( Z=rep( c(0,1), each=3 ), 
                    S = rep( c("e", "hq", "lq" ), 2 ),
                    W = 0,
                    n = 0 )
SHELL_FRAME$cat = with( SHELL_FRAME, paste( Z, S, sep="." ) )



# Calculate the weighted means and relative sizes of each group
calc.statistics = function( Yobs, S, Z, w ) {
    require( plyr )

    if ( is.data.frame( Yobs ) ) {
        df = Yobs
    } else {
        df = data.frame( Yobs = Yobs, S = S, Z = Z, w = w )
    }
    
    sumstat = ddply( df, .( Z, S ), summarize,
                     n = length( w ),
                     W = sum( w ),
                     Y.bar = sum( w * Yobs ) / sum( w ) )
    sumstat$cat = paste( sumstat$Z, sumstat$S, sep="." )
    missing = setdiff( c( "0.e", "0.hq", "0.lq", "1.e", "1.hq", "1.lq" ), sumstat$cat )
    
    sumstat = rbind.fill( sumstat, subset( SHELL_FRAME, cat %in% missing ) )
    rownames(sumstat) = sumstat$cat

    sumstat = ddply( sumstat, .(Z), function( dd ) {
        dd$p = dd$W / sum(dd$W)
        dd
    })
    rownames( sumstat ) = sumstat$cat
    sumstat$cat = NULL
    sumstat
}


# Calculate the weighted means and relative sizes of each group
# given  _Complete_ info, i.e.
# Y0, Y1, S0, S1 
calc.statistics.oracle = function( Y0, Y1, S0, S1 ) {
    df = data.frame( Z = rep( c(0,1), each=length(Y0) ),
                     Yobs = c( Y0, Y1 ),
                     S = c( as.character(S0), as.character( S1 ) ),
                     w = 0.5 )
    calc.statistics( df )
}


# Given full knowledge of strata membership and all P.Os,
# calculate parameters of interest.
calc.params = function( dat ) {
    
    rs = ddply( dat, .(R), summarize, mu.0 = mean( Y0 ),
                mu.1 = mean( Y1 ),
                N = length( Y1 ) )
    rs$pi = rs$N / sum( rs$N )
    rs$ITT = rs$mu.1 - rs$mu.0 
    rs
}


# Make sure non-NA elements of vector v lies in [0,1], truncating if not
# @return vector with each element in [0,1] or NA
constrain.interval = function( v, low = 0, high=1 ) {
    pmax( low, pmin( v, high ) )
}


# Take the estimated means of all the groups and return a list of possible
# ITT.lqc, ITT.hqc pairs that are possible
# @param stats  Dataframe with column 'strata' corresponding to the strata and a ITT.hat.high
#    and ITT.hat.low for the high and low bounds for the strata.  Looks for lqc and hqc
get.tradeoff.vector = function( stats ) {
    rownames( stats ) = stats$strata
    highs = stats[ c("hqc","lqc" ), "ITT.hat.high" ]
    lows = stats[ c("hqc", "lqc"), "ITT.hat.low" ]
    
    rs = rbind( c( highs[1], lows[2] ), c( lows[1], highs[2] ) )
    colnames( rs ) = c( "ITT.hqc", "ITT.lqc" )
    rs
}


# Given stats describing the min and max, plot a line indicating all
# possible bounds under those constraints.
# @param stats  Dataframe with a 'strata' column corresponding to the strata and a ITT.hat.high
#    and ITT.hat.low for the high and low bounds for the strata.  Looks for lqc and hqc
plot.bounds = function( stats, xlim=c(-1,1), ylim=c(-1,1), asp=1,
                        xlab=expression("ITT"["lqc"]), ylab=expression("ITT"["hqc"]), bty="n",
                        add=FALSE, col="black", ... ) {
    rownames( stats ) = stats$strata
    tr = get.tradeoff.vector( stats )
    
    if ( !add ) {
        plot( NA, xlim=xlim, ylim=ylim, asp=asp, xlab=xlab, ylab=ylab, ... )
    }
    
    #if ( stats["hqc","N.s"] == 0 || stats["lqc","N.s"]==0 ) {
    #    text( 0, 0, "no compliers of one or both groups\nto plot" )
    #} else {
    mu.hat.0 = stats$mu.hat.0
    names( mu.hat.0 ) = rownames( stats )
    if ( stats["hqc","N.s"] == 0 ) {
        ITT.hqc.high = 1
        ITT.hqc.low = -1
        tr[,"ITT.hqc"] = c(-1,1)
    } else {
        ITT.hqc.high = stats[ "hqc", "ITT.hat.high"]
        ITT.hqc.low = stats[ "hqc", "ITT.hat.low"]
    }
    if ( stats["lqc","N.s"]==0 ) {
        ITT.lqc.high = 1
        ITT.lqc.low = -1
        tr[,"ITT.lqc"] = c(-1,1)
    } else {
        ITT.lqc.high = stats[ "lqc", "ITT.hat.high"]
        ITT.lqc.low = stats[ "lqc", "ITT.hat.low"]
    }
        abline( v = 0, col="grey" )
        abline( h = 0, col="grey" )
        abline( h=c( ITT.hqc.high, ITT.hqc.low ), col="red")
        abline( v=c( ITT.lqc.high, ITT.lqc.low ), col="red" )
        
        lines( tr[,1] ~ tr[,2], type="l", lwd=6, col=col )    
   
}


# As plot bounds, but for the treatment means not the ITT
plot.bounds.mean = function( stats, xlim=c(0,1), ylim=c(0,1), asp=1,
                             xlab="LQC mean (tx)", ylab="HQC mean (tx)", add=FALSE, col="black", ... ) {
    rownames( stats ) = stats$strata
    
    if ( !add ) {
        plot( NA, xlim=xlim, ylim=ylim, asp=asp, xlab=xlab, ylab=ylab, ... )
    }
    
    if ( stats["hqc","N.s"] == 0 || stats["lqc","N.s"]==0 ) {
        text( mean(xlim), mean(ylim), "no compliers of one or both groups\nto plot" )
    } else {
        
        hqc.high = stats[ "hqc", "mu.hat.1.high"]
        hqc.low = stats[ "hqc", "mu.hat.1.low"]
        lqc.high = stats[ "lqc", "mu.hat.1.high"]
        lqc.low = stats[ "lqc", "mu.hat.1.low"]
        
        abline( v = 0, col="grey" )
        abline( h = 0, col="grey" )
        abline( h=c( hqc.high, hqc.low ), col="red")
        abline( v=c( lqc.high, lqc.low ), col="red" )
        
        segments( x0 = lqc.low, y0 = hqc.high, x1 = lqc.high, y1 = hqc.low, lwd=6, col=col )    
    }
}


calc.bounds.from.table = function( stats ) {
    p = stats$p
    names( p ) = rownames(stats)
    Y.bar = stats$Y.bar
    names( Y.bar ) = rownames(stats)
    
    ## Calculate proportions
    pis = c( p[ "0.e" ],
             p[ "1.hq" ],
             p[ "0.hq" ] - p[ "1.hq" ],
             p[ "1.lq" ],
             p[ "0.lq" ] - p[ "1.lq" ] )
    names( pis ) = STRATA_NAMES
    pis[ is.na( pis ) | (pis < 0 ) ] = 0
    
    #stopifnot( all( pis > 0 ) )
    #if ( any( pis <= 0 ) ) {
    #    browser()
    #}
    
    # don't let nonexistence of ATs kill any bounds
    # (But make value big to catch bugs if we multiply by nonzero mass)
    Y.bar[ p == 0 ] = 100000
    
    # Calculate identifiable group means
    mu.hat.0 = rep( NA, 5 )
    names( mu.hat.0 ) = STRATA_NAMES
    
    mu.hat.0[ "lqat" ] = Y.bar[ "1.lq" ]
    mu.hat.0[ "hqat" ] = Y.bar[ "1.hq" ]
    mu.hat.0[ "eat" ]  = Y.bar[ "0.e" ]
    mu.hat.1 = mu.hat.0
    
    mu.hat.0[ "lqc" ] = ( p[ "0.lq" ] / (p["0.lq"] - p["1.lq"] )) * Y.bar["0.lq"] - (p["1.lq"]/(p["0.lq"] - p["1.lq"] ))*Y.bar["1.lq"]
    mu.hat.0[ "hqc" ] = ( p[ "0.hq" ] / (p["0.hq"] - p["1.hq"] )) * Y.bar["0.hq"] - (p["1.hq"]/(p["0.hq"] - p["1.hq"] ))*Y.bar["1.hq"]
    
    mu.hat.0 = constrain.interval( mu.hat.0 )
    
    # Calculate the bounds using formulae
    B.lqc = (p["1.e"] / (p["0.lq"]-p["1.lq"])) * Y.bar["1.e"] - (p["0.e"]/(p["0.lq"]-p["1.lq"])) * Y.bar["0.e"]
    mu.1.lqc.low = B.lqc - (p["0.hq"] - p["1.hq"])/(p["0.lq"]-p["1.lq"])
    mu.1.lqc.high = B.lqc
    
    B.hqc = (p["1.e"] / (p["0.hq"]-p["1.hq"])) * Y.bar["1.e"] - (p["0.e"]/(p["0.hq"]-p["1.hq"])) * Y.bar["0.e"]
    mu.1.hqc.low = B.hqc - (p["0.lq"] - p["1.lq"])/(p["0.hq"]-p["1.hq"])
    mu.1.hqc.high = B.hqc
    
    mu.hat.1.high = mu.hat.1.low = mu.hat.1
    mu.hat.1.low[ c("lqc","hqc") ] = c( mu.1.lqc.low, mu.1.hqc.low )
    mu.hat.1.high[ c("lqc","hqc") ] = c( mu.1.lqc.high, mu.1.hqc.high )
    mu.hat.1.low = constrain.interval( mu.hat.1.low )
    mu.hat.1.high = constrain.interval( mu.hat.1.high )
    
    ITT.hat.low = mu.hat.1.low - mu.hat.0
    ITT.hat.high = mu.hat.1.high - mu.hat.0
    
    # Force undefined bounds to be maximal
    ITT.hat.low[ is.na( ITT.hat.low ) ] = -1
    ITT.hat.high[ is.na( ITT.hat.high ) ] = 1
    
    stat2 = data.frame( strata = STRATA_NAMES,
                        pi = pis,
                        mu.hat.0 = mu.hat.0,
                        mu.hat.1 = mu.hat.1,
                        mu.hat.1.low = mu.hat.1.low,
                        mu.hat.1.high = mu.hat.1.high,
                        ITT.hat.low = ITT.hat.low,
                        ITT.hat.high = ITT.hat.high,
                        N.s = pis * sum( stats$W ),
                        N = sum( stats$W ) )

    # Table of bounds given the observed statistics
    stat2
}


# Main function: Calculate the bounds for a collection of observations
# @param Yobs 0/1 outcome.  S compliance behavior.  Z treatment assignment,
#          w is sampling weight of unit.
# @return Dataframe with the various quantities calculated for the strata
#    of interest.
calc.bounds = function( Yobs, S, Z, w, 
                        return.type=c("main","intervals") ) {
    return.type = match.arg(return.type)
        
    # hack so we can pass a dataframe with the Yobs, S, Z, w columns defined
    if ( is.data.frame( Yobs ) ) {
         S = Yobs$S
         Z = Yobs$Z
         w = Yobs$w
         Yobs = Yobs$Yobs
    }
    
    stats = calc.statistics( Yobs, S, Z, w )
    
    stat2 = calc.bounds.from.table( stats )
    
    if ( return.type=="intervals" ) {
        get.tradeoff.vector( stat2 )   
    } else {
        stat2
    }      
}






# Calculate the bounds for a population with complete potential outcome table
# @param Y0, Y1 0/1 outcomes.  S0, S1 compliance behavior under Z=0, 1. 
# @return Dataframe with the various quantities calculated for the strata
#    of interest.
calc.bounds.oracle = function( Y0, Y1, S0, S1, 
                        return.type=c("main","intervals") ) {
    return.type = match.arg(return.type)
    
    # hack so we can pass a dataframe with the Yobs, S, Z, w columns defined
    if ( is.data.frame( Y0 ) ) {
        Y1 = Y0$Y1
        S0 = Y0$S0
        S1 = Y0$S1
        Y0 = Y0$Y0
    }
    
    stats = calc.statistics.oracle( Y0, Y1, S0, S1 )

    stat2 = calc.bounds.from.table( stats )
    
    if ( return.type=="intervals" ) {
        get.tradeoff.vector( stat2 )   
    } else {
        stat2
    }          

}



# Plot some results!
plot.bounds.strat = function( stats, ..., nrow=2 ) {
    stats = subset( stats, B != "TOTAL" )
        K = length( unique( stats$B ) )
        if ( K > 4*nrow ) {
            plot = FALSE
            warning( "Too many subpanels to plot" )
        } else {
            par( mfrow=c(nrow, trunc((K+nrow-1)/nrow) ) )
        }
    
    d_ply( stats, .(B), function( df ) {
        plot.bounds( df, main=df$B[[1]], ... )
    } )
}


# Calculate bounds when we stratify and then re-weight the strata.
calc.bounds.strat = function( Yobs, S, Z, w, B, 
                              return.type=c("main","full","components","intervals") ) {
    return.type = match.arg(return.type)
    
    require( plyr)

    if ( is.data.frame( Yobs ) ) {
        df = Yobs
        stopifnot( !is.null( df$S ) && !is.null( df$Z ) && !is.null( df$w ) )
        stopifnot( !is.null( df$B) )
    } else {
        df = data.frame( Yobs=Yobs, S=S, Z=Z, w=w, B=B )
    }
    
    # Calculate bounds for each slice
    bnds = ddply( df, .(B), function( dd ) {
        #cat( "Number of rows: ", nrow(dd), "\n" )
        rs = with( dd, calc.bounds( Yobs, S, Z, w ) )
        rs$N = nrow( dd )
        rs
    } )
    
    # add in strata weights for posterity
    bnds = ddply( bnds, .(strata), transform, weight = N.s / sum( N.s ) )

    
    #browser()
    #mg[ c( 5, 3, 1, 4, 2 ), ]
    if ( return.type=="components" ) {
        bnds
    } else {
        mg = ddply( bnds, .(strata), function( dd ) {
            
            a = sapply( dd[c("mu.hat.0","mu.hat.1","mu.hat.1.low","mu.hat.1.high","ITT.hat.low","ITT.hat.high")],
                        function( d ) { 
                            sum( (dd$weight * d)[dd$weight > 0] ) 
                        } )
            a["N"] = sum( dd$N )
            a["N.s"] = sum( dd$N.s )
            a["pi"] = sum( dd$N.s ) / sum(dd$N)
            a["weight"] = sum( dd$weight )
            a
        } )
        
        if ( return.type=="full" ) {
            mg$B = "TOTAL"
            rs = rbind.fill( bnds, mg )
            rs$B = factor( rs$B, levels = c( levels( bnds$B ), "TOTAL" ) )
            rs
        } else if ( return.type=="intervals") {
            get.tradeoff.vector( mg )
        } else {
            mg
        }
    }
}

# For simulation, if we know all potential outcomes we can calculate the true bounds.
calc.bounds.strat.oracle = function( Y0, Y1, S0, S1, B, 
                              return.type=c("main","full","components","intervals") ) {
    return.type = match.arg(return.type)
    
    require( plyr)
    
    if ( is.data.frame( Y0 ) ) {
        df = Y0
    } else {
        df = data.frame( Y0 = Y0, Y1 = Y1, S0=S0, S1=S1, B=B )
    }
    
    # Calculate bounds for each slice
    bnds = ddply( df, .(B), function( dd ) {
        #cat( "Number of rows: ", nrow(dd), "\n" )
        rs = with( dd, calc.bounds.oracle( dd ) )
        rs$N = nrow( dd )
        rs
    } )
    
    # add in strata weights for posterity
    bnds = ddply( bnds, .(strata), transform, weight = N.s / sum( N.s ) )
    
    if ( return.type=="components" ) {
        bnds
    } else {
        mg = ddply( bnds, .(strata), function( dd ) {
            
            a = sapply( dd[c("mu.hat.0","mu.hat.1","mu.hat.1.low","mu.hat.1.high","ITT.hat.low","ITT.hat.high")],
                        function( d ) { 
                            sum( (dd$weight * d)[dd$weight > 0] ) 
                        } )
            a["N"] = sum( dd$N )
            a["N.s"] = sum( dd$N.s )
            a["pi"] = sum( dd$N.s ) / sum(dd$N)
            a["weight"] = sum( dd$weight )
            a
        } )
        
        if ( return.type=="full" ) {
            mg$B = "TOTAL"
            rs = rbind.fill( bnds, mg )
            rs$B = factor( rs$B, levels = c( levels( bnds$B ), "TOTAL" ) )
            rs
        } else if ( return.type=="intervals") {
            get.tradeoff.vector( mg )
        } else {
            mg
        }
    }
}


# @param dataframe 'data' with columns given by the variables Yobs, S, Z, w, B
#     (It will look up variables automatically.)
# @param R number of bootstrap iterations
boot.bounds = function( Yobs, S, Z, w, B = NULL, data, R=100, save.samples=TRUE ) {
    
    Yobs = eval( substitute( Yobs ), data )
    S = eval( substitute( S ), data )
    Z = eval( substitute( Z ), data )
    w = eval( substitute( w ), data )
    B = eval( substitute( B ), data )
    stratify = FALSE
    if ( !is.null(B) ) {
        stratify=TRUE
    }
    if ( stratify ) {
        cb_func = calc.bounds.strat
        dat = data.frame( Yobs=Yobs, S=S, Z=Z, w=w, B=B )
    } else {
        cb_func = calc.bounds        
        dat = data.frame( Yobs=Yobs, S=S, Z=Z, w=w )
    }
    
    bnd = cb_func( dat, return.type="main" )
    
    runs = 1:R
    names(runs) = runs
    res = laply( runs, function( index ) {
        dat.star = dat[ sample( nrow(dat), replace=TRUE ), ]
        rsf = cb_func( dat.star, return.type="main" )
        rs =         get.tradeoff.vector( rsf )
        pies =rsf$pi
        names( pies ) = rownames( rsf )
        c( ITT.hqc.low = rs[[2,1]], ITT.hqc.high=rs[[1,1]], 
           ITT.lqc.low=rs[[1,2]], ITT.lqc.high=rs[[2,2]], pi=pies )
        #rs
        #cbind( Bnd=c("A","B"), as.data.frame( rs ) )
    } )
    
    head( res )
    res = cbind( res, hqc.width=res[,2]-res[,1], 
                      lqc.width=res[,4]-res[,3],
                      disjoint = ( res[,2] < res[,3] | res[,4] < res[,1] ) )
    
    CIs = adply( res, 2, function( v ) {
        c( mn = mean(v), q=quantile(v, c(0.05,0.95) ), SE=sd( v ) )
    }, .id = "bound" )
    CIs


    
    CIs$pe = c( bnd$ITT.hat.low[3], bnd$ITT.hat.high[3], bnd$ITT.hat.low[5], bnd$ITT.hat.high[5], bnd$pi, NA, NA, NA)
    CIs$pe[10] = CIs$pe[2] - CIs$pe[1]
    CIs$pe[11] = CIs$pe[4] - CIs$pe[3]
    CIs$strata = c( rep( c("hqc","lqc"), each=2 ),  as.character( bnd$strata ), "hqc", "lqc", "joint")
    CIs$type = c( rep( c("low","high"), 2 ), rep( "pi", 5 ), "width", "width", "overlap" )

    attr( CIs, "bounds" ) <- bnd
    if ( save.samples ) {
        attr( CIs, "samples" ) <- res
    }
    CIs
}

plot.boot.bounds = function( CIs, xlim=c(-1,1), ylim=c(-1,1), asp=1,
                             xlab=expression("ITT"["lqc"]), ylab=expression("ITT"["hqc"]), bty="n",
                             add=FALSE, col="black", ... ) {
    bnd = attr( CIs, "bounds" )
    plot( NA, xlim=xlim, ylim=ylim, asp=asp, xlab=xlab, ylab=ylab, ... )

    abline( v=c( CIs$`q.5%`[3], CIs$`q.95%`[4] ), col="green" )
    abline( h=c( CIs$`q.5%`[1], CIs$`q.95%`[2] ), col="green" )
    res = attr( CIs, "samples" )
    if ( !is.null( res ) ) { 
        a_ply( res, 1, function( v ) {
            lines( v[1:2] ~ v[c(4,3)], col="grey" )
        } )
    }
    plot.bounds( bnd, add=TRUE )
}


# Take list of CIs on the endpoints and make a nice display of the confidence
# intervals.
convert.boot.bounds = function( CIs ) {
    ddply( subset( CIs, strata %in% c("lqc","hqc") ), .(strata), function( df ) {
        pes = df$pe
        SEs = df$SE
        names( SEs ) = names( pes ) = df$type
        c( low=df[1,3], high=df[2,4], pe=pes, SE=SEs )
    })
}


# Clean up table and print it with only essential information.
nice.print = function( res ) {
    rs = subset( res, strata %in% c( "hqc", "lqc" ) )
    rs$N = rs$weight = NULL
    rs$mu.hat.1 = NULL
    print( rs, digits=3, row.names=FALSE )
}


# Takes the result of a sliced bound calculation (with the individual 
# strata information kept via return.type="full") and generate a table of 
# the estimated proportions of each strata in each slice.
get.pi.by.strata = function( estimate.table ) {
    
    # library( plyr )
    # 
    # ddply( stat2b, .(B), function( df ) {
    #     pis = as.numeric( df$pi )
    #     names( pis ) = df$strata
    #     pis
    # })
    # 
    
    require( reshape2 )
    dcast( estimate.table[ c(1:3) ],  B ~ strata, value.var="pi" )
    
}