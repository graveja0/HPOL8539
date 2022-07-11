
## Script to examine at how uncertainty in bound estimation grows with increased
## levels of stratification
##
## (C) Luke Miratrix, 2018   lmiratrix@gse.harvard.edu


set.seed( 10191 )

library( xtable )
library( dplyr )

source( "bound_function_library.R" )
source( "data_generators.R" )
source( "set_paper_parameters.R" )

##
## gen a finite sample population of data
##
N = 1000
dat = gen.dat.strat( K=1, N=N, pies = pi.vec, beta=beta.vec, 
                 noise.sd.1 = noise.1, noise.sd.2 = noise.2, noise.sd.3 = noise.3 )



params = calc.params( dat )
params

Ks = c( 1:15, 20, 25, 30, 35)
rs = ldply( Ks, function( k ) {
    dat = strat.dat( dat, k )
    stat2b = with( dat, calc.bounds.strat( Yobs, S, Z, w, B=B1 ) )
    rs = (stat2b$ITT.hat.high - stat2b$ITT.hat.low)[c(3,5)]
    names(rs) = c( "hqc.width", "lqc.width" )
    rs
} )
rs$K = Ks

plot( hqc.width ~ K, data=rs )


rs2 = ldply( Ks, function( k ) {
    dat = strat.dat( dat, k )
    dat$B = dat$B1
    CIs = boot.bounds( Yobs, S, Z, w, B, data=dat, R=20  )
    rs = c( hqc.width = CIs[ 2, "q.95%" ] - CIs[ 1, "q.5%" ],
            lqc.width = CIs[ 4, "q.95%" ] - CIs[ 3, "q.5%" ] )
    # rs = stat2b$`q.95%`[c(5,6)]
    #names( rs ) = c( "hqc.width", "lqc.width" )
    rs
} )
rs2
rs2$K = Ks

par( mfrow=c(1,2), mar=c(3,3,2,0.5) )
plot( hqc.width ~ K, ylim=range( rs$hqc.width, rs2$hqc.width), data=rs, main="HQC", bty="n" )
points( hqc.width ~ K, data=rs2, pch=19 )

plot( lqc.width ~ K, ylim=range( rs$lqc.width, rs2$lqc.width), data=rs, main="LQC", bty="n" )
points( lqc.width ~ K, data=rs2, pch=19 )




##
## Investigations
##

dat = strat.dat( dat, 10 )
stat2b = with( dat, calc.bounds.strat( Yobs, S, Z, w, B=B1, return.type = "full" ) )
ss = subset( stat2b, strata=="hqc" )
head( ss )
nrow( ss )

par( mfrow=c(1,1) )
plot( ITT.hat.low ~ as.numeric(B), data=ss,ylim=c(-1,1), type="b" )
lines( ITT.hat.high ~ as.numeric(B), data=ss, type="b",ylim=c(-1,1), col="red" )
