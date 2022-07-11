

## Paper Illustration of the different bounds for each of the 6 strata
## to show how the trade-off of tightness of bounds varies by strata.
##
## For a specific dataset, slice based on covariates and make a plot of how
## wide the bounds are within each slice.
## 
## Do this for principle and prognostic scores.
##
## Note: we changed this script to one that uses the ECHS data so the illustration
## is more interesting.
##
## (C) Luke Miratrix, 2018   lmiratrix@gse.harvard.edu


##
## Exploring bounds
##

set.seed( 1019 )

library( xtable )
library( dplyr )

source( "bound_function_library.R" )

source( "data_generators.R" )

source( "set_paper_parameters.R" )

##
## gen a finite sample population of data
##

dat = gen.preset.dat()



params = calc.params( dat )
params

table( dat$B1, dat$R, dat$Z )
table( dat$B1 )


stats = with( dat, calc.statistics( Yobs, S, Z, w ) )

print( xtable( stats, caption="Observed quantities", 
               label="tab:observed_data", digits=c(0,0,0,0,0,2,2 ) ) )
stats



##
## Calculate bounds (without stratification)!
##

stat2 = with( dat, calc.bounds( Yobs, S, Z, w ) )
stat2




##
## principal score stratification!
##


cat( "Using the strata\n" )
stat2b = with( dat, calc.bounds.strat( Yobs, S, Z, w, B=B1, 
                                       return.type = "full" ) )

pdf( file=gen_file_name( "demo_plot"), width=8, height=5 )
plot.bounds.strat( stat2b, bty="n" )
dev.off()



########################################################################
##
## Make a table showing relative sizes of principal strata
## within each slice
##
########################################################################

stat2b = with( dat, calc.bounds.strat( Yobs, S, Z, w, B=B1, 
                                       return.type = "full" ) )

head( stat2b )

tb = ddply( stat2b, .(B), function( dd ) { 
    n = dd$weight
    # browser()
    names(n) = as.character( dd$strata )
    n
} )

tb

tb2 = ddply( stat2b, .(B), function( dd ) {
    pis = dd$pi
    names(pis) = dd$strata
    c( pi.hqc = pis[["hqc"]], pi.lqc = pis[["lqc"]], ratio = pis[[ "hqc" ]] / pis[[ "lqc" ]] )
})
tb2
names( tb2 ) = c( "B", "pi.hqc", "pi.lqc", "ratio" )

lvl = levels( tb$B )
tb = merge( tb, tb2, by="B" )
tb$B = factor( tb$B, levels=lvl )
tb = tb[ order(tb$B), ]
tb

tb$pi.hqc = paste( round( 100*tb$pi.hqc), "%", sep="" )
tb$pi.lqc = paste( round( 100*tb$pi.lqc), "%", sep="" )
tb$ratio = round( tb$ratio, digits= 1)

tb

library( xtable )
print( xtable( tb, caption="Weights of each strata for each slice", 
               label="tab:weight-chart" ),
       include.rownames=FALSE )




##
## Plot the bounds given stratification
##

stat2b
par( mfrow=c( 1, 2), mar=c(3,3,0.5,0.5), mgp=c(1.7,0.8,0))
plot.bounds( stat2 )
plot.bounds( subset( stat2b, B=="TOTAL" ) )

# Plot on top of each other to show shrinkage
par( mfrow=c( 1, 1) )
plot.bounds( stat2 )
plot.bounds( subset( stat2b, B=="TOTAL" ), col="blue", add=TRUE )






##
## Prognostic Example: when we slice based on prognostic score how do the
## different bounds within a slice look?
##

# re-gen data with no noise for prognostic
# noise.2 = 0
#dat = gen.strat.dat.outcome( K=6, N=N, pies = pi.vec, beta=beta, 
#                             noise.sd.1 = noise.1, noise.sd.2 = noise.2, noise.sd.3 = noise.3 )


cat( "Using the prognostic strata\n" )

stat2b = with( dat, calc.bounds.strat( Yobs, S, Z, w, B=B3, 
                                       return.type = "full" ) )

pdf( file=gen_file_name( "strat_demo_prognostic_plot" ), width=8, height=5 )
plot.bounds.strat( stat2b, bty="n" )
dev.off()

subset( stat2b, strata=="lqc" )


# Look at mean outcome by strata and treatment assignment
library( reshape2 )
dd = ddply( dat, .(B3, Z), summarize, N = length( Yobs ), Ybar = mean( Yobs ) )
head( dd )
dcast( dd, B3 ~ Z, value.var="Ybar" )


