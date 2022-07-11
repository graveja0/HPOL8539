##
##
## This script allows one to explore and investigate the properties of a simulated dataset
## with ease.
##
## This demo code also shows how the bundle of code works that accompanies the paper
## Bounding, An Accessible Method for Estimating Principal Causal Effects, Examined and Explained
## by Luke Miratrix, Jane Furey, Avi Feller, Todd Grindal & Lindsay C. Page
##
## It explores three different stratification methods (principle, complier, and prognostic)
##
## It generates a data set and then calculates different bounds using different 
## covariates.
##
## This is a good script to play around with the bounds for a single data set.
##
## (C) Luke Miratrix, 2018   lmiratrix@gse.harvard.edu


library( xtable )
library( dplyr )

source( "bound_function_library.R" )

source( "data_generators.R" )

source( "set_paper_parameters.R" )


##
## gen a finite sample population of data
##

noise = 0
dat = gen.preset.dat( pop.K=K, pop.N=N, 
                      noise.sd.1 = noise, noise.sd.2 = noise, noise.sd.3 = noise )



##
## Look at population parameters
##

params = calc.params( dat )
params

prop.hqc = params$pi[3]
prop.lqc = params$pi[5]
prop.hqc / prop.lqc
prop.lqc / prop.hqc
1 - prop.hqc / (prop.hqc + prop.lqc )



##
## Look at observed statistics
##

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
## Now try bounds with stratification!
##


cat( "Using the strata\n" )
table( dat$B1, dat$R, dat$Z )
table( dat$B1 )
stat2b = with( dat, calc.bounds.strat( Yobs, S, Z, w, B=B1, return.type = "full" ) )
stat2b

get.pi.by.strata( stat2b )


cat( "Using the compliance predictor\n" )
table( dat$B2, dat$R )
table( dat$B2 )
stat2c = with( dat, calc.bounds.strat( Yobs, S, Z, w, B=B2, return.type = "full" ) )
stat2c



cat( "Using the Y (prognostic predictor)\n" )
table( dat$B3, dat$R )

# different levels of outcome for different strata
group_by( dat, B3 ) %>% summarize( mean.Y = mean( Yobs ) )

stat2d = with( dat, calc.bounds.strat( Yobs, S, Z, w, B=B3, return.type = "full" ) )
stat2d
subset( stat2d, strata=="hqc" )

library( reshape2 )
dd = ddply( dat, .(B3, Z), summarize, N = length( Yobs ), Ybar = mean( Yobs ) )
dd
dcast( dd, B3 ~ Z, value.var="Ybar" )



cat( "Using both\n" )
dat$Bx = with( dat, interaction( B1, B3 ) )
table( dat$Bx )
table( dat$Bx, dat$R )
stat2e = with( dat, calc.bounds.strat( Yobs, S, Z, w, B=Bx, return.type = "full" ) )
stat2e


##
## Compare bounds for compliers across the different methods
##


stat2$method = "raw"
stat2$B = "TOTAL"
stat2b$method = "principal"
stat2c$method = "compliance"
stat2d$method = "prognostic"
stat2e$method = "princ+prog"

stat.all = rbind.fill( stat2b, stat2c, stat2d, stat2, stat2e  )
stat.all

stat.all$weight = NULL
stat.all$width = with( stat.all, ITT.hat.high - ITT.hat.low )

print( params )

print( subset( stat.all, strata =="hqc" & B=="TOTAL" ) )

print( subset( stat.all, strata == "lqc" & B=="TOTAL") )



