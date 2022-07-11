
##
## Exploring bounds
## Looking at stratification on the three different Xs
## 
## See how bound width changes as function of the noise (i.e., the strength
## of connection between observed covariates and outcome, strata membership, etc.)
##
## This script generates the primary figures in the ECHS paper.
##
## (C) Luke Miratrix, 2018   lmiratrix@gse.harvard.edu

library( xtable )
#library( dplyr )

source( "bound_function_library.R" )

source( "data_generators.R" )

source( "set_paper_parameters.R" )

RUN_ON_SOURCE = TRUE


# The various degrees of noise to explore
#noise = round( seq( 0, 3, length.out=1000 ), digits=2 )

# Prior simulations found that the relationship between R2 and noise was
#  (1/sqrt(R2)) = a + b * noise (i.e., was linear in shape, roughly)
# To get equal intervals we then did the following (for ECHS data)
a = 1 #0.77
b = 1.19
desired.R2 = #round( 4 * seq( 0.01, 1.00, length.out=40 ), digits = 1 ) / 4
desired.R2 = seq( 0.01, 1.00, by=0.01 )
    
desired.R2[ desired.R2 == 0 ] = 0.01
length( unique( desired.R2 ))
noise = round( (1/b) * (1/sqrt(desired.R2) - a ), digits=2 )
tail( noise )
#noise = c( noise, 0 )
length( noise )
noise = rep( noise, 25 )
#noise = rep( noise, 2 )

one.run = function( noise, K ) {

    cat( "\nNoise level = ", noise, "\n" )
    dat = gen.preset.dat( pop.K=K, pop.N=N, 
                          noise.sd.1 = noise, noise.sd.2 = noise, noise.sd.3 = noise )

    #print( table( dat$Yobs, dat$Z ) )

    stat = with( dat, calc.bounds( Yobs, S, Z, w ) )
    stat

    stat.1 = with( dat, calc.bounds.strat( Yobs, S, Z, w, B=B1  ) )
    stat.2 = with( dat, calc.bounds.strat( Yobs, S, Z, w, B=B2  ) )
    stat.3 = with( dat, calc.bounds.strat( Yobs, S, Z, w, B=B3 ) )
    
    # both
    dat$Bx = with( dat, interaction( B1, B3 ) )
    stat.4 = with( dat, calc.bounds.strat( Yobs, S, Z, w, B=Bx ) )

    stat$method = "raw"
    stat.1$method = "X1"
    stat.2$method = "X2"
    stat.3$method = "X3"
    stat.4$method = "X1+X3"
    
    stat.all = rbind.fill( stat, stat.1, stat.2, stat.3, stat.4 )
    stat.all
    stat.all$method = factor( stat.all$method, levels = c( "raw","X1","X2","X3","X1+X3" ))
    
    stat.all$weight = NULL
    stat.all$width = with( stat.all, ITT.hat.high - ITT.hat.low )

    stat.all
}






## Do single testing bit
if ( FALSE ) {
    dd = gen.dat.strat( K=6, N=10000, pies = pi.vec, 
                                mean.outcome.co = mean.outcome.co, mean.outcome.tx=mean.outcome.tx, beta=beta,
                                noise.sd.1 = 1, noise.sd.2 = 1, noise.sd.3 = 1 )
    head( dd )
    
    dat = gen.dat.strat( K=6, N=3000, pies = pi.vec, beta = beta, 
                                 mean.outcome.co=mean.outcome.co,
                                 mean.outcome.tx=mean.outcome.tx,
                                 noise.sd.1 = noise, noise.sd.2 = noise, noise.sd.3 = noise )
    calc.params( dat )
    
    stat.all = one.run( 0, K )

    subset( stat.all, strata=="hqc" )
}





########   Main Simulation Running Code (with plot-making code ) #######

if ( RUN_ON_SOURCE ) {
    library( tidyverse )
    
    names( noise ) = noise 
    rs = ldply( noise, one.run, K=K, .id="noise" )
    rs$noise = as.numeric( as.character( rs$noise ) )
    levels( rs$method ) = c( "No strat", "X1: principal", "X2: compliance", "X3: prognostic", "X1+X3: princ + prog" )
    head( rs )
    
    # calculate the R2s for the different levels of noise to give a better
    # x-axis to our plot.
    unoise = unique( noise )
    R2s = ldply( unoise, function( noise ) {
        cat( "Gen R2 for noise = ", noise, "\n" )
        dat = gen.preset.dat( pop.K=K, pop.N=10*N, 
                              noise.sd.1 = noise, noise.sd.2 = noise, noise.sd.3 = noise )
        get.R2s(dat)
    }, .id="noise" )
    R2s$noise  = unoise
    #R2s$noise = as.numeric( as.character( R2s$noise ) )
    
    subset( R2s, noise %in% c( 0, 1, 2, 3 ) )
    
    # examine R2 for noise
    R22 = tidyr::gather( R2s, R2.princ, R2.comp, R2.prog, key="X", value="R2" )
    nrow(R22)
    head(R22)
    ggplot( R22, aes(noise, R2, col=X ) ) + geom_point() + geom_smooth( se=FALSE) 
        #scale_y_log10()
    ggplot( R22, aes(y=noise, x=R2, col=X ) ) + geom_point() #+ geom_smooth( se=FALSE,  alpha=0.2, )
    #scale_y_log10()
    ggplot( R22, aes(noise, I(1/sqrt(R2)), col=X ) ) + geom_point() + geom_smooth( se=FALSE) # +
    summary( lm( I(1/sqrt(R2)) ~ noise + X, data = R22 ) )                              
    
    
    
    print( subset( rs, as.numeric( noise ) == 3 ) )
    
    rs = merge( rs, R2s, by="noise", all.x = TRUE, all.y=FALSE )
    rs$R2 = NULL
    levels( rs$method )
    head( select( rs, method, starts_with("R2" ) ) )
    rs = rs %>% group_by( method ) %>% 
                mutate( R2 = switch( as.numeric( method[[1]] ), 
                              (R2.princ + R2.comp + R2.prog) / 3,
                              R2.princ,
                              R2.comp,
                              R2.prog,
                              (R2.princ + R2.prog)/2 ) )
    ggplot( rs, aes( x=noise, y=R2, col=method ) ) + geom_line()
    head( rs )                  
                    
    rss = subset( rs, strata %in% c("lqc","hqc" ) )
    head( rss )
    
    # theme_bw()
    library( ggthemes )
                    
    my_t = theme_calc() + theme( legend.position="bottom", legend.direction="horizontal", legend.key.width=unit(1,"cm"),
                                 panel.border = element_blank() )
    theme_set( my_t )
    plt = ggplot( rss, aes(R2, width, group=method, color=method, lty=method ) ) + 
        facet_grid( . ~ strata, scales="free_y" ) +
        geom_jitter( cex=0.8, alpha=0.05 ) +
        geom_smooth( size=1, se=FALSE) +
        xlab("Approximate R2 of Covariates") + ylab("Width of Interval" )
    
    #plt + guides(fill=guide_legend(
    #    keywidth=30,
    #    default.unit="cm")
    #)
    print( plt )
    

    filename = gen_file_name( "performance" )
    cat( "Saving to", filename,"\n" )
    ggsave( plt, file=filename )
    
    filename = gsub( ".pdf", ".RData", filename )
    filename
    save.image( file=filename )
}

# percent improvement plot
if ( RUN_ON_SOURCE ) {
    head( rss )    
    levels( rss$method )
    table( rss$method )
    widths = filter( rss, as.numeric(method) == 1 ) %>% 
        group_by( noise, strata ) %>% 
        dplyr::summarise( max.width = mean( width ) )
    widths    
    ggplot( widths, aes( noise, max.width, col=strata ) ) + geom_line()
    rss = merge( rss, widths, by=c("noise","strata"))
    head(rss)
    rss = mutate( rss, per.width = 100 * width / max.width )
    
    library( ggthemes )
    
    my_t = theme_calc() + theme( legend.position="bottom", legend.direction="horizontal", legend.key.width=unit(1,"cm"),
                                 panel.border = element_blank() )
    theme_set( my_t )
    plt = ggplot( filter( rss, method != "No strat"), aes(R2, per.width, group=method, color=method ) ) + 
        facet_grid( . ~ strata, scales="free_y" ) +
        geom_jitter( cex=0.8, alpha=0.10 ) +
        geom_hline( yintercept=100) +
        geom_smooth( size=1, se=FALSE) +
        xlab("Approximate R2 of Covariates") + ylab("Percent width of interval compared to no adjustment" ) +
        scale_y_continuous( breaks=seq(0,125,by=25))
    print( plt )
    filename = gen_file_name( "performance_percent")
    cat( "Saving to", filename,"\n" )
    ggsave( plt, file=filename )
    
}


# Presentation plots
if ( FALSE ) {
    
    # Plots for presentations/talks (wider lines, etc)
    # (Seperate LQC and HQC plots.)
    ggplot( subset( rss, strata=="lqc" ), aes(noise, width, group=method, color=method) ) + 
                      #facet_grid( . ~ strata, scales="free_y" ) +
                      geom_jitter( size=1.2 ) +
                      geom_smooth( size=2, se=FALSE) +
                        xlab("Amount of noise in covariates") + ylab("Width of Interval" )
    
    ggplot( subset( rss, strata=="hqc" ), aes(noise, width, group=method, color=method) ) + 
        #facet_grid( . ~ strata, scales="free_y" ) +
        geom_jitter( size=1.2 ) +
        geom_smooth( size=2, se=FALSE) +
        xlab("Amount of noise in covariates") + ylab("Width of Interval" )
    
}

