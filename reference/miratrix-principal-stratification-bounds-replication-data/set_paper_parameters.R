
##
## This script sets the parameter settings that all scripts use.
##
## This makes it easy to run different scripts on the same DGP.
##
## (C) Luke Miratrix, 2018   lmiratrix@gse.harvard.edu

# I'm putting these variables in one file so I don't accidentally make mistakes
# with my examples not aligning with the simulation study

source( "data_generators.R" )

VERSION = c( "default", "ECHS", "simple", "manystrata" )
VERSION = VERSION[[1]]

cat( "\n\nSimulation variant: ", VERSION, "\n\n" )

#           "eat"  "hqat"   "hqc"  "lqat"   "lqc" 
pi.vec = c( 0.1,    0.05,   0.45,   0.15,    NA )
pi.vec[5] = 1 - sum( pi.vec[1:4] )
names( pi.vec ) = STRATA_NAMES

#  coefficient for latent variables to outcome model (for U1, U2, and U3)
beta.vec = c( 0, 0, 4 )

#  impact of different school environments
#school.shift = c( -1, 0, 1 )
mu.Y.co = c( 1, 0, 0, -1, -1 )
mu.Y.tx = c( 1, 0, 1, -1,  1 )

#  sample size
N = 4000

# proportion to treat
prop.tx.default = 0.5

##
## Baseline noise (sd of measurement error on the latent variables) when noise not varied.
##

# Principal
noise.1 = 0

# compliance
noise.2 = 0

# prognostic
noise.3 = 0



##
## For data analysis
##

# Number of slices
K = 6

## 
## For saving results
##

FILE_ROOT = "results/"
FILE_SUFFIX = ""

gen_file_name = function( name ) {
    paste( FILE_ROOT, name, FILE_SUFFIX, ".pdf", sep="" )
}




##
## ECHS-like version
## These parameter settings mimic the moment estimates from the ECHS study
##
if ( VERSION == "ECHS" ) {
    cat( "\nDoing ECHS version\n" )
    K = 4
    
    #           eat     hqat      hqc     lqat      lqc 
    mu.Y.co = c( 17,     8,   13,      6,     6.5 ) 
    mu.Y.tx = mu.Y.co + c( 0, 0, -2, 0, +3.1 )
    
    #school.shift = c( 2.5, 4, 5 )
    
    pi.vec = c( 0.03, 0.03, 0.11, 0.11, 0.72 )
    names( pi.vec ) = STRATA_NAMES

    # make beta very predictive
    beta.vec = c( 0, 0, 5.65 )
    
    prop.tx.default = 0.58
    N = 3820
    
    FILE_SUFFIX = "_echs"
}

if ( VERSION == "simple" ) {
    cat( "\nDoing simple version (no always-takers)\n" )
    
    pi.vec = c( 0, 0, 0.7, 0, 0.3 )
    
    mu.Y.co = c( 5, 3, 4, 1.5, 2.5 )
    mu.Y.tx = c( 5, 3, 6, 1.5, 4 )
    
    names( pi.vec ) = STRATA_NAMES

    FILE_SUFFIX = "_simple"
    
}

if ( VERSION == "manystrata" ) {
    K = 15
    FILE_SUFFIX = "_many"
}



## Utility to generate a population and randomization with the set parameters above
gen.preset.dat = function( pop.K = K, pop.N = N, pies = pi.vec, beta = beta.vec,
                           prop.tx = prop.tx.default,
                           mean.outcome.co = mu.Y.co, mean.outcome.tx = mu.Y.tx,
                           noise.sd.1 = noise.1, noise.sd.2 = noise.2, noise.sd.3 = noise.3,
                           message=TRUE) {
    if ( message ) {
        cat( "Generating the default dataset (size = ", pop.N, ")...\n", sep="" )
    }
    
    gen.dat.strat( K=pop.K, N=pop.N, pies = pies, beta = beta, 
                   mean.outcome.co = mean.outcome.co, mean.outcome.tx = mean.outcome.tx,
                   prop.tx = prop.tx,
                   noise.sd.1 = noise.sd.1, noise.sd.2 = noise.sd.2, noise.sd.3 = noise.sd.3 )
}


##
## Get some parameters on an extra-large sample and print
## them out for reference.
##

dat = gen.preset.dat( pop.N=N*10 )


cat( "Parameters of the population being generated\n" )

cat( "pi:\n" )
print( pi.vec )

cat( "all parameters:\n" )
params = calc.params( dat )
params$pi = round( params$pi, digits=3 )
print( params )

cat( "ITT = ", mean( dat$Y1 - dat$Y0 ), "\n" )
cat( "M0 = ", with( subset( dat, R %in% c("lqc","hqc") ), mean( Y0 ) ), "\n" )
cat( "M1 = ", with( subset( dat, R %in% c("lqc","hqc") ), mean( Y1 ) ), "\n" )
cat( "CACE = ", with( subset( dat, R %in% c("lqc","hqc") ), mean( Y1 - Y0 ) ), "\n" )

cat( "\nSample filename: '", gen_file_name( "sample"), "'\n", sep="" )

head( dat )

get.R2s = function( dat ) {
    dd = dplyr::filter( dat, Z==0, S !="e" )
    table( dd$R )
    dd$goLQ = dd$S == "lq"
    mod <- glm(goLQ~X1, data=dd, family="binomial")
    nullmod <- glm(goLQ~1, data=dd, family="binomial")
    R2.1 = 1-logLik(mod)/logLik(nullmod)
    
    dd = dplyr::filter( dat, Z==1 )
    dd$isComp = dd$S == "e"
    mod <- glm(isComp~X2, data=dd, family="binomial")
    nullmod <- glm(isComp~1, data=dd, family="binomial")
    R2.2 = 1-logLik(mod)/logLik(nullmod)
    
    dd = dplyr::filter( dat, Z==0 )
    #dd = dplyr::filter( dat, Z==0, S=="e" )
    mod <- glm(Yobs~X3+S, data=dd, family="binomial")
    nullmod <- glm(Yobs~S, data=dd, family="binomial")
    R2.3 = 1-logLik(mod)/logLik(nullmod)
    
    c( R2.princ = R2.1, R2.comp = R2.2, R2.prog = R2.3 )
}

if (FALSE ) {
    # examine no-noise strength of covariates
    summary( lm( Yobs ~ X1, data=dat ) )
    summary( lm( Yobs ~ X2, data=dat ) )
    summary( lm( Yobs ~ X3, data=dat ) )

    noise = 0
    dat = gen.preset.dat( pop.N=N*10, noise.sd.1=noise, noise.sd.2=noise, noise.sd.3=noise )
    
    get.R2s( dat )
    
    library( tidyverse )
    head( dat )
    dd = filter( dat, S != "e", Z==0 )
    dd = dplyr::filter( dat, Z==0, S=="e" )
    table( dd$R )
    ggplot( dd, aes( X3, Y0.p, col=S ) ) + geom_point( alpha=0.3 ) + geom_smooth()
    
    
    group_by( dat, Z, B3 ) %>% summarize( meanY = mean( Yobs ) )
}

cat( "Prop Tx = ", mean( dat$Z ), "\n" )
cat( "R2 with noise = ", noise.1, noise.2, noise.3, "\n" )
print( get.R2s( dat ) )
