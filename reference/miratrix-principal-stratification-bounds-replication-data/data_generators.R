##
## Functions to generate simulation data for doing bounds calculations on
##
## Primary function is gen.dat.param() and gen.strat.dat.outcome()
##
## (C) Luke Miratrix, 2018   lmiratrix@gse.harvard.edu

source( "bound_function_library.R" )
library( plyr )



# Generate data with predictive U1, U2, U3.
# This is the DGP in the paper simulation study.
#
# These U's predict different things.
# U1 - High or low quality 
# U2 - being an always-taker or not (high U2 is echs always-taker) and 
#     low U2 is always-taker of other schools
# U3 - Predicts outcome (logistic model)
# @param beta  The slope for the response of Y to U1, U2, U3
# Note: Nonzero beta_1, beta_2 will make the U1, U2 variables also predictive of outcome.
# mean.outcome.co/tx - mean outcomes for different strata under tx and co.  Could violate exclusion
#    restriction with these parameter settings.
# Removed-> School shift (impact of going to a school type) in order of lq, hq, echs
gen.dat.param = function( N, 
                          pies = rep( 0.2, 5 ),
                          mean.outcome.co = c(1, 0, 0, -1, -1),
                          mean.outcome.tx = c(1, 0, 1, -1, 1),
                          beta = c( 0, 0, 2 ),
                          prop.tx = 0.5 ) {


    #if ( is.null( names( mean.outcome.co ) ) ) {
    #    names( mean.outcome.co ) = STRATA_NAMES    
    #}
    #if ( is.null( names( mean.outcome.tx ) ) ) {
    #    names( mean.outcome.tx ) = STRATA_NAMES    
    #}
    names( mean.outcome.co ) = names( mean.outcome.tx ) = STRATA_NAMES    
    
    delm = mean.outcome.tx - mean.outcome.co
    if( delm["eat"] != 0 || delm["hqat"] != 0 || delm["lqat"] != 0 ) {
        warning( "Exclusion restriction is violated!" )
    }
    #cat( "Effects:\n")
    #print( delm )
    
    if ( is.null( names( pies ) ) ) {
        names( pies ) = STRATA_NAMES
        print( pies )
    }
    #names( school.shift ) = c("lq","hq","e" )
    
    
    U1 = rnorm( N )
    U2 = rnorm( N )
    U3 = rnorm( N )
    
    
    # Given not EAT, go to hq school?
    phq = (pies["hqc"]+pies["hqat"])/(1-pies["eat"])
    h_or_l = as.numeric( U1 > qnorm( 1 - (pies["hqc"]+pies["hqat"])/(1-pies["eat"]) ) )
    
    # are they EATs?
    go_echs.0 = as.numeric( U2 > qnorm( 1 - pies["eat"] ) )  
    
    pi.hq_cut = pies["hqc"]/(pies["hqc"]+pies["hqat"]) * (1 - pies["eat"] ) + pies["eat"]
    hq_cut = qnorm( 1 - pi.hq_cut )
    
    pi.lq_cut = pies["lqc"]/(pies["lqc"]+pies["lqat"]) * (1 - pies["eat"] ) + pies["eat"]
    lq_cut = qnorm( 1 - pi.lq_cut )

    # complier or always-taker?
    go_echs.1 = as.numeric( ifelse( h_or_l, U2 > hq_cut, U2 > lq_cut ) )
                                    #U2 > qnorm( 1 - pies["lqc"] - pies["eat"]*(1-phq) ) ) )
                                    #U2 > qnorm( 1 - pies["lqc"]/(pies["lqc"]+pies["lqat"]) - pies["eat"] ) ) )
    
    S.0 = ifelse( go_echs.0, "e", ifelse( h_or_l, "hq", "lq" ) )
    S.1 = ifelse( go_echs.1, "e", ifelse( h_or_l, "hq", "lq" ) )
    R = ifelse( S.0=="e", "eat", ifelse( S.0=="hq", ifelse( S.1=="e", "hqc", "hqat" ),
                                         ifelse( S.1=="e", "lqc", "lqat" ) ) )
    
    Y0.p = beta[[1]] * U1 + beta[[2]] * U2 + beta[[3]] * U3 + mean.outcome.co[R]
    Y1.p = beta[[1]] * U1 + beta[[2]] * U2 + beta[[3]] * U3 + mean.outcome.tx[R]
    
    Y0.p = exp(Y0.p)/(1+exp(Y0.p))
    Y1.p = exp(Y1.p)/(1+exp(Y1.p))
    
    n.tx = round( N * prop.tx )
    Z = rep( c(0,1), c(N-n.tx, n.tx ) )
    
    S = ifelse( Z, S.1, S.0 )

    # get actual POs
    Us = runif( length( Y0.p ) )
    Y0 = as.numeric( Us <= Y0.p )
    Y1 = as.numeric( Us <= Y1.p )
    #  Yobs = 0 + as.numeric( runif( length( Yobs ) ) <= Yobs )
    
    Yobs = ifelse( Z, Y1, Y0 )
    
    
    data.frame( U1=U1, U2=U2, U3=U3, Z = Z, S=S, Yobs=Yobs, w=1, R=R, 
                Y0 = Y0, Y1=Y1, Y0.p = Y0.p, Y1.p = Y1.p, S.0 = S.0, S.1 = S.1 )
}






strat.dat = function( dat, K ) {
    dat$B1 = cut( dat$X1, breaks=quantile( dat$X1, 0:K / K ), include.lowest=TRUE )
    dat$B2 = cut( dat$X2, breaks=quantile( dat$X2, 0:K / K ), include.lowest=TRUE )
    dat$B3 = cut( dat$X3, breaks=quantile( dat$X3, 0:K / K ), include.lowest=TRUE )
    
    dat
    
}

# Generate parameterized data and then generate observed
# strata based on noisy versions of the latent parameters
gen.dat.strat = function( K, N, pies, 
                                  mean.outcome.co = c(1, 0, 0, -1, -1),
                                  mean.outcome.tx = c(1, 0, 1, -1, 1),
                                  beta,
                                  prop.tx = 0.5,
                                  noise.sd.1=0, noise.sd.2=0, noise.sd.3=0 ) {
    
    dat = gen.dat.param( N, pies, mean.outcome.co=mean.outcome.co, mean.outcome.tx=mean.outcome.tx,
                         beta=beta, prop.tx=prop.tx )
    
    dat$X1 = round( dat$U1 + rnorm( nrow(dat), sd=noise.sd.1), digits = 2 )
    dat$X2 = round( dat$U2 + rnorm( nrow(dat), sd=noise.sd.2 ), digits = 2 )
    dat$X3 = round( dat$U3 + rnorm( nrow(dat), sd=noise.sd.3 ), digits = 2 )

    # and slice and return resulting data.
    strat.dat( dat, K )
}






# Reassign treatment and regenerate the observed behavior and outcome
# Utility for simulation studies
rerandomize.dat = function( dat ) {
    dat$Z = sample( dat$Z )
    dat$Yobs = ifelse( dat$Z, dat$Y1, dat$Y0 )
    dat$S = ifelse( dat$Z, as.character(dat$S.1), as.character( dat$S.0 ) )
    dat
}



# testing
if ( FALSE ) {
    td = gen.strat.dat( 3, mu.0 = rep( 0.1,5), mu.1 = rep(0.4,5) )
    head( td )
    
    calc.params( td )
}



if  (FALSE ) {
    
    # testing
    ig  = gen.dat( pi=c(0.1,0.1,0.1,0.1,0.6), mu.0 = rep(0.1,5), mu.1=rep(0.9,5) )
    table( ig$R )
    calc.params( ig )
    head( ig )
    
    dt = gen.dat.param( 20 )
    head( dt )
    calc.params( dt ) 
    
    dt = gen.dat.param( 100000, pies = c( 0.1, 0.40, 0.3, 0.15, 0.05 ) )
    rs = calc.params( dt ) 
    rs
    rs$pi[2]+rs$pi[3]
    rs$pi[4]+rs$pi[5]
    head( dt )
    
    ddply( dt, .(R), summarize, Y0 = mean(Y0), 
                    Y1 = mean( Y1 ),
                    tau = Y1 - Y0 )
    
    boxplot( U2 ~ R, data=dt )
    boxplot( U1 ~ R, data=dt )
    boxplot( U3 ~ R, data=dt )
    
    
}

