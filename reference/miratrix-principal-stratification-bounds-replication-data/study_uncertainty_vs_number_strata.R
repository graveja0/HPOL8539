

## Simulation to look at how uncertainty in bound estimation grows with
## increased levels of stratification.
##
## We then compare the estimated uncertainty to the _actual_ uncertainty
## generated from oracle estimators.
##
## WARNING: This is a long script to run as it does several replications of a
## bootstrap procedure which is time consuming
##
## (C) Luke Miratrix, 2018   lmiratrix@gse.harvard.edu

set.seed( 10191 )

library( xtable )
library( dplyr )

source( "bound_function_library.R" )
source( "data_generators.R" )
source( "set_paper_parameters.R" )

RUN_SIM_ON_SOURCE = TRUE

# if we run a sim on source, which sim?
RUN_BOOT_SIM = TRUE  

Ks = c( 1, 3, 5, 10, 15, 30 )  #c( 1:15, 20, 25, 30, 35)
#names( Ks ) = Ks

# number of boot replicates within each simulation
BR = 10  # WARNING: for real life, set this to 100 or more

# number of simulation runs
N_SIM = 10  # WARNING: for real life, set this to 100 or more


noise = 0.5
#dat = gen.preset.dat( pop.K=K, pop.N=N * 10, 
#                      noise.sd.1 = noise, noise.sd.2 = noise, noise.sd.3 = noise,
#                      message=FALSE)

#get.R2s(dat)


# Take data and generate bounds for a sequence of Ks
# if 'include.oracle' then also do it with the potential outcomes, ignoring tx
# assignment.  Otherwise just estimates.
one.run = function( dat, Ks, include.oracle = TRUE ) {
    
    # get bound estimates for each K
    rs = ldply( Ks, function( k ) {
        dat = strat.dat( dat, k )
        stat2b = with( dat, calc.bounds.strat( Yobs=Yobs, S=S, Z=Z, w=w, B=B1 ) )
        rs = (stat2b$ITT.hat.high - stat2b$ITT.hat.low)[c(3,5)]
        #names(rs) = c( "hqc.width", "lqc.width" )
        #rs
        a = stat2b[ c(3,5), c( "strata", "ITT.hat.low", "ITT.hat.high")]
        data.frame( K=k,
                    strata=c( "hqc", "lqc", "hqc", "lqc" ), 
                    side=c("low","low","high","high"), 
                    value=c(a[[2]], a[[3]]),
                    stringsAsFactors = FALSE )        
    }, .id=NULL )
    #rs$kind = "est"
    

    # calculate parameters    
    pr = calc.params( dat )
    #pr = pr[ rep( c(3,5), length(Ks)), c("R", "ITT", "ITT" ) ]
    #pr$kind = "param"
    #pr = cbind( K=rep(Ks,each=2), pr )
    #names(pr) = names(rs)
    
    if ( include.oracle ) {
        # get oracle bounds for each K
        rs.o = ldply( Ks, function( k ) {
            dat = strat.dat( dat, k )
            stat2b = with( dat, calc.bounds.strat.oracle( Y0, Y1, S.0, S.1, B=B1 ) )
            a = stat2b[ c(3,5), c( "strata", "ITT.hat.low", "ITT.hat.high")]
            data.frame( K=k, 
                        strata=c( "hqc", "lqc", "hqc", "lqc" ), 
                        side=c("low","low","high","high"), 
                        value=c(a[[2]], a[[3]]),
                        stringsAsFactors = FALSE )
        }, .id=NULL )
        #rs.o$kind = "oracle"
        rs$oracle = rs.o$value #rs = rbind( rs, rs.o )
    }
    rs$param = pr$ITT[c(3,5)]
    
    rs
    
    #rs = transform( rs, err.low = ITT.hat.low - pr$ITT.hat.low,
    #                    err.high =  pr$ITT.hat.high - ITT.hat.high )
    
    #rs = merge( rs, rs.o, by="K", all=TRUE )
    #rs
    #rbind.fill( rs, pr )
}

if ( FALSE ) {
    dat = gen.preset.dat( pop.K=K, pop.N=N, 
                          noise.sd.1 = noise, noise.sd.2 = noise, noise.sd.3 = noise )
    
    one = one.run( dat, Ks )
    one
}

######### Bootstrap check ##################

one.run.boot = function( Ks ) {
    
    dat = gen.preset.dat( pop.K=K, pop.N=N, 
                          noise.sd.1 = noise, noise.sd.2 = noise, noise.sd.3 = noise,
                          message=FALSE)

    main = one.run( dat, Ks )
    
    # get bootstrap adjusted bounds
    runs = 1:BR
    names(runs) = runs
    res = ldply( runs, function( index ) {
        dat.star = dat[ sample( nrow(dat), replace=TRUE ), ]
        rs = one.run( dat.star, Ks, include.oracle=FALSE )
    } )
    
    r2 = ddply( res, .(K,strata,side), summarize, boot.SE = sd( value ),
                                             boot.mean = mean( value ),
                                             boot.CI.l = quantile( value, c(0.05) ),
                                             boot.CI.h = quantile( value, c(0.95 ) ) )
    
    r3 = merge( main, r2 )
    r3 = r3[ order( r3$K ), ]
    
    r3
}

if ( FALSE ) {
   
    system.time( one.b <- one.run.boot( Ks ) )
    one.b
    
    plot( value ~ oracle, data=one.b, col=ifelse( side=="high", "red","green" ), pch=ifelse(strata=="hqc",21,22) )    
}


################# Main Simulation Code ##################################

##
## Do the sim study (with bootstrapping)
##
if ( RUN_SIM_ON_SOURCE && RUN_BOOT_SIM ) {
    cat( "\nRunning simulation (with bootstrap)\n" )
    
    rns = ldply( 1:N_SIM, function( i ) { 
        one.run.boot( Ks )
    }, .progress = "text" )
    
    head( rns )
    
    # are the point estimates biased?
    sts = ddply( rns, .(K,strata,side), summarize, bias = mean( value - oracle ),
                 SE = sd( value ),
                 SE.hat = mean( boot.SE ),
                 miss.below = mean( boot.CI.l >= oracle ),
                 miss.above = mean( boot.CI.h <= oracle ),
                 miss = mean(  boot.CI.l >= oracle |  boot.CI.h <= oracle ),
                 p.below.tau = mean( value <= param ),
                 p.above.tau = mean( value >= param ))
    subset( sts, K==3 )
    
    
    # reshape so we can examine intervals
    head( rns )
    highs = 2 * (1:(nrow(rns)/2)) - 1
    lows = 2 * (1:(nrow(rns)/2))
    
    rns.lows = rns[ lows, -c(1,2,3,6,7,8,10) ]
    head( rns.lows )
    rns.highs = rns[ highs, -c(1,2,3,6,7,8,9) ]
    head( rns.highs )
    head( rns[ highs, c(1,2,6) ] )
    rnsInt = cbind( rns[ highs, c(1,2,6) ],  
                    rns.lows, 
                    rns.highs )
    head( rnsInt[ c(4,5,7,8) ] )
    names( rnsInt )[ c(4,5,7,8) ] = c( "est.l", "oracle.l", "est.h", "oracle.h" )
    head( rnsInt )
    
    # full coverage of entire interval
    #rnsInt$inCI = with( rnsInt, boot.CI.l <= oracle.l & oracle.h <= boot.CI.h )
    
    # coverage of parameter
    rnsInt$inCI.boot = with( rnsInt, boot.CI.l <= param & param <= boot.CI.h )
    rnsInt$inCI = with( rnsInt, est.l <= param & param <= est.h )
    
    head( rnsInt )
    rnsInt = mutate( rnsInt, width = est.h - est.l,
                     oracle.width = oracle.h - oracle.l,
                     boot.width = boot.CI.h - boot.CI.l )
    head( rnsInt )
    
    
    
    
    cat( "Saving results of all the hard work\n" )
    save.image( file=gen_file_name("study_uncertainty_sim_results_boot" ) )
}


##
## Further analysis of sim results
##
if ( FALSE ) {
    load( file=gen_file_name( "study_uncertainty_sim_results_boot" ) )

    stsInt = ddply( rnsInt, .(K, strata), summarize, 
                    coverage=mean( inCI ),
                    mean.width = mean( width ),
                    mean.oracle = mean( oracle.width ),
                    mean.boot = mean( boot.width ),
                    sd.width = sd( width ),
                    sd.oracle = sd( oracle.width ) )
    
    stsInt = mutate( stsInt, per.off = mean.width / mean.oracle )
    arrange( stsInt, strata, K )
    
    library( tidyverse )
    s2 = gather( stsInt, mean.width, mean.oracle, mean.boot, key="type", value="width" )
    tb = arrange( stsInt, strata, K )
    
    head( s2 )
    plt <- ggplot( s2, aes(K, width, col=type) ) +
        facet_wrap(~ strata, scales="free" ) +
        geom_point() + geom_line( size = 1) +
        labs( x="K (number of slices)", y="mean interval width")
    print( plt )
    filename = gen_file_name( "strata_study")
    cat( "Saving to", filename,"\n" )
    ggsave( plt, file=filename )

    
        
    plot( bias ~ K, data=sts ) 
    abline( h=0 )

    plot( SE.hat ~ K, data=sts ) 
    plot( SE  ~ K, data=sts ) 
    
    head( sts )
    ggplot( filter( s2, K>=2 ), aes( width, sd.width ) ) + 
        facet_wrap( ~ strata, scales="free" ) + geom_point()
    ggplot( filter( s2, K>=2 ), aes( y=sd.width, x=K ) ) + 
        facet_wrap( ~ strata, scales="free" ) + geom_point()
    
    
    
    hist( subset( rns, side=="high" & strata=="hqc" & K==3 )$value, breaks=30 )
    
    library( ggplot2 )
    library( ggthemes )
    my_t = theme_calc() + theme( legend.position="bottom", legend.direction="horizontal", legend.key.width=unit(1,"cm")  )
    theme_set( my_t )
    ggplot( data=subset( rns, strata=="hqc"), aes(value) ) +
            facet_grid( side ~ K, scales = "free" ) +
                geom_histogram( bins=10) +
        labs( title="Distribution of the bound endpoint estimates" )
    

    head( sts )
    plt = ggplot( sts, aes(K, SE, group=side, color=side ) ) + 
        facet_grid( . ~ strata, scales="free_y" ) + 
        geom_line( ) + geom_point() +
        xlab("Number of slices") + ylab("SE of bound estimate" )
    print( plt )
    
    head( stsInt )
    nrow( stsInt )
    plt = ggplot( stsInt, aes(K, mean.width ) ) + 
        facet_grid( . ~ strata, scales="free_y" ) + 
        geom_point( ) +
        xlab("Number of slices") + ylab("Width of Interval" )
    print( plt )
    
}





##
## Just sim study without bootstrapping
##
if ( RUN_SIM_ON_SOURCE && !RUN_BOOT_SIM ) {
    cat( "\nRunning simulation (with NO bootstrap)\n" )
    
    # Without boot
    rns = rdply( 2, { 
        dat = gen.preset.dat( pop.K=K, pop.N=N, 
                              noise.sd.1 = noise, noise.sd.2 = noise, noise.sd.3 = noise )
        
        one.run( dat, Ks )
    }, .id="runID" )
    
    head( rns )
    
    library( tidyverse )
    
    rns = gather( rns, value, oracle, param, key="kind",value="value" )
    head( rns )
    rns = spread( rns, side, value )
    
    # Before I did the gather above, I used to do this.
    # this fails?
    #reshape2::dcast( r2, runID+K+strata~side, value.var=c("value", "oracle","param"), sep="")
    #rns = reshape( rns, idvar=c("runID","K","strata"),
    #        v.names=c("value","oracle","param"),
    #        timevar="side",
    #        direction="wide" )

    head( rns )
    rns$width = with( rns, high - low )
    
    library( ggplot2 )
    library( ggthemes )
    # theme_bw()
    my_t = theme_calc() + theme( legend.position="bottom",
                                 legend.direction="horizontal", 
                                 legend.key.width=unit(1,"cm")  )
    theme_set( my_t )
    plt = ggplot( rns, aes(K, width, group=kind, color=kind ) ) + 
        facet_grid( . ~ strata, scales="free_y" ) + 
        geom_point( ) +
        xlab("Number of slices") + ylab("Width of Interval" )
    
    print( plt )
    
}



