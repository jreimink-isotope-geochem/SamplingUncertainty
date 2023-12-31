
## Range of Sample Sizes

fte_theme_white <- function() {
  library( RColorBrewer )
  palette          <- brewer.pal( "Greys", n = 9 )
  color.grid.major = palette[ 4 ]
  color.axis.text  = palette[ 7 ]
  color.axis.title = palette[ 8 ]
  color.title      = palette[ 9 ]
  theme_bw( base_size = 9 ) +
    theme( panel.grid.major = element_line( color = color.grid.major, size = 0.25 ) ) +
    theme( axis.line.x      = element_line( color = color.grid.major, size = 0.35 ) ) +
    theme( axis.line.y      = element_line( color = color.grid.major, size = 0.35 ) ) +
    theme( panel.grid.minor = element_blank( ) ) +
    theme( axis.ticks       = element_blank( ) ) +
    theme( legend.text       = element_text( size = 7, color = color.axis.title ) ) +
    theme( plot.title   = element_text( color = color.title, size = 25, vjust = 1.25, hjust = 0.5 ) ) +
    theme( plot.subtitle   = element_text( color = color.title, size = 16, vjust = 1.25, hjust = 0.5 ) ) +
    theme( axis.text.x  = element_text( size = 28, color = color.axis.text ) ) +
    theme( axis.text.y  = element_text( size = 28, color = color.axis.text ) ) +
    theme( axis.title.x = element_text( size = 36, color = color.axis.title, vjust = 0 ) ) +
    theme( axis.title.y = element_text( size = 36, color = color.axis.title, vjust = 1.25 ) ) +
    theme( plot.margin = unit( c( 0.35, 0.2, 0.3, 0.35 ), "cm" ) )
}

shiny.min.parameters.raw <- read.csv( "Shiny_OutPut.csv", stringsAsFactors = FALSE)

shiny.min.parameters <- subset( shiny.min.parameters.raw, Modes > 0 )

#############################  DATA IMPORT #############################
# Import the modes, densities, and molar masses
modes <- t( shiny.min.parameters$Modes )
densities <- t( shiny.min.parameters$Densities )
# mol.masses <- t( shiny.min.parameters$Weights )
grain.dims <- t( cbind( shiny.min.parameters$Dimension1, shiny.min.parameters$Dimension2, shiny.min.parameters$Dimension3 )  )
min.dists <- t( shiny.min.parameters$Distribution )
min.comps <- t( cbind( shiny.min.parameters[, c( which( colnames(shiny.min.parameters)=="SiO2" ):ncol(shiny.min.parameters))] ) )

## Set column names to be minerals
colnames(modes) = shiny.min.parameters$X
colnames(densities) = shiny.min.parameters$X
# colnames(mol.masses) = shiny.min.parameters$X
colnames(grain.dims) = shiny.min.parameters$X
colnames(min.dists) = shiny.min.parameters$X
colnames(min.comps) = shiny.min.parameters$X



## set sample mass
sample.masses <- c( 625, 1250, 2500, 5000, 10000, 20000 ) #in grams
N_n <- 1000


for( j in 1:length(sample.masses)) {
  # set sample mass
  
  sample.mass = sample.masses[j]
  ## this all replicates the first columns in Stanley spreadsheet
  # inverse.molar.volume <- 1 / ( mol.masses / densities ) 
  modes.internal <- modes * densities 
  modes.wt.total <- rowSums( modes.internal )
  modes.wt <- modes.internal / modes.wt.total * 100 
  # modes.moles.internal <- modes * inverse.molar.volume 
  # tot.modes.moles.internal <- rowSums( modes.moles.internal )
  # modes.mole <- modes.moles.internal / tot.modes.moles.internal * 100 
  grain.vols <- apply( grain.dims, 2, prod )
  grain.mass <- grain.vols * densities
  tot.min.mass <- sample.mass * modes.wt / 100 
  num.min.grains <- tot.min.mass / grain.mass
  tot.num.min.grains <- num.min.grains/( modes/100 )
  ## calculate error, using the hypergeometric calc in Stanley spreadsheet
  
  err.num.min.grains.h <- sqrt(( 100-modes)/100*modes/100*tot.num.min.grains*
                                 (N_n*tot.num.min.grains-tot.num.min.grains)/(N_n*tot.num.min.grains-1))
  err.num.min.grains.b <- sqrt(( 100-modes)/100*modes/100*tot.num.min.grains)
  err.num.min.grains.p <- sqrt( num.min.grains )
  
  rsd.err.num.min.grains.h <- err.num.min.grains.h / num.min.grains * 100
  rsd.err.num.min.grains.b <- err.num.min.grains.b / num.min.grains * 100
  rsd.err.num.min.grains.p <- err.num.min.grains.p / num.min.grains * 100
  
  ## create a random set of mineral grains in a given rock, then recalculate the mode in wt% of each grain
  nreps = 10000
  grain.num.model <- data.frame( matrix( nrow = nreps, ncol = ncol( num.min.grains ) ) )
  colnames( grain.num.model ) <- colnames( num.min.grains )
  
  
  
  
  ## model the number of mineral grains in a synthetic rock
  for (col in 1:ncol( grain.num.model ) ) {
    # col=11
    # Extract mean and sd for current column
    mean_val <- num.min.grains[ 1, col ]
    # calculate SD based on the distribution of the grains
    sd_val <- ifelse( min.dists[col] == "h", err.num.min.grains.h[1, col],
                      ifelse( min.dists[col] == "b", err.num.min.grains.b[1, col],
                              err.num.min.grains.p[1, col] ) )
    
    # Generate resampled value based on mean and sd - using Poisson but others are normally
    # distributed with the SD values being defined by the Stanley spreadsheet. So, effectively normal
    if( min.dists[col] == "p" ){
      resampled_val <- rpois(n = nreps, lambda = mean_val)
    } else {
      resampled_val <- rnorm( nreps, mean = mean_val, sd = sd_val )
    } 
    
    # resampled_val <- rnorm( nreps, mean = mean_val, sd = sd_val )
    # Assign resampled value to corresponding column in the empty dataframe
    grain.num.model[, col] <- resampled_val
  }
  ## convert number of grains to mass of each grain
  grain.mass.model <- sweep( grain.num.model, 2, unlist( grain.mass[1, ] ), `*` )
  
  grain.num.summary <- data.frame( means = apply( grain.num.model, 2, mean ),
                                   stdev = apply( grain.num.model, 2, sd ) )
  
  ## remove all negative values
  grain.mass.model <- replace( grain.mass.model, grain.mass.model < 0, NA )
  
  
  ## calculate the wt% of each mineral
  # Calculate the column totals
  mineral.model.tot.mass <- rowSums( grain.mass.model, na.rm = T )
  
  # Normalize the columns and add them to the dataframe
  mineral.model.normalized <- grain.mass.model %>%
    mutate( across( everything(), ~ . / mineral.model.tot.mass ) )
  
  ## summarize the grain masses
  grain.mass.summary <- data.frame( means = apply( mineral.model.normalized, 2, mean ),
                                    stdev = apply( mineral.model.normalized, 2, sd ) )
  grain.mass.summary$rsd <- grain.mass.summary$stdev / grain.mass.summary$means * 100
  
  ## calculate mineral comps - this is in the same order as the elemental oxides listed in the min.comps file (row orders)
  minerals.comp.model <- list( SiO2 = sweep( mineral.model.normalized, 2, unlist( min.comps["SiO2", ] ), `*` ),
                               TiO2 = sweep( mineral.model.normalized, 2, unlist( min.comps["TiO2", ] ), `*` ),
                               Al2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps["Al2O3", ] ), `*` ),
                               Cr2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps["Cr2O3", ] ), `*` ),
                               Fe2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps["Fe2O3", ] ), `*` ),
                               FeO = sweep( mineral.model.normalized, 2, unlist( min.comps["FeO", ] ), `*` ),
                               MnO = sweep( mineral.model.normalized, 2, unlist( min.comps["MnO", ] ), `*` ),
                               MgO = sweep( mineral.model.normalized, 2, unlist( min.comps["MgO", ] ), `*` ),
                               CaO = sweep( mineral.model.normalized, 2, unlist( min.comps["CaO", ] ), `*` ),
                               Na2O = sweep( mineral.model.normalized, 2, unlist( min.comps["Na2O", ] ), `*` ),
                               K2O = sweep( mineral.model.normalized, 2, unlist( min.comps["K2O", ] ), `*` ),
                               O = sweep( mineral.model.normalized, 2, unlist( min.comps["O", ] ), `*` ),
                               H2O = sweep( mineral.model.normalized, 2, unlist( min.comps["H2O", ] ), `*` ),
                               P2O5 = sweep( mineral.model.normalized, 2, unlist( min.comps["P2O5", ] ), `*` ),
                               ZrO2 = sweep( mineral.model.normalized, 2, unlist( min.comps["ZrO2", ] ), `*` ),
                               CO2 = sweep( mineral.model.normalized, 2, unlist( min.comps["CO2", ] ), `*` ),
                               SO3 = sweep( mineral.model.normalized, 2, unlist( min.comps["SO3", ] ), `*` ),
                               Li2O = sweep( mineral.model.normalized, 2, unlist( min.comps["Li2O", ] ), `*` ),
                               ThO2 = sweep( mineral.model.normalized, 2, unlist( min.comps["ThO2", ] ), `*` ),
                               BaO = sweep( mineral.model.normalized, 2, unlist( min.comps["BaO", ] ), `*` ),
                               La2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps["La2O3", ] ), `*` ),
                               Ce2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps["Ce2O3", ] ), `*` ),
                               Pr2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps["Pr2O3", ] ), `*` ),
                               Nd2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps["Nd2O3", ] ), `*` ),
                               Sm2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps["Sm2O3", ] ), `*` ),
                               Eu2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps["Eu2O3", ] ), `*` ),
                               Gd2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps["Gd2O3", ] ), `*` ),
                               Tb2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps["Tb2O3", ] ), `*` ),
                               Dy2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps["Dy2O3", ] ), `*` ),
                               Ho2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps["Ho2O3", ] ), `*` ),
                               Er2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps["Er2O3", ] ), `*` ),
                               Tm2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps["Tm2O3", ] ), `*` ),
                               Yb2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps["Yb2O3", ] ), `*` ),
                               Lu2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps["Lu2O3", ] ), `*` ),
                               Y2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps["Y2O3", ] ), `*` ),
                               Sb2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps["Sb2O3", ] ), `*` ),
                               CuO = sweep( mineral.model.normalized, 2, unlist( min.comps["CuO", ] ), `*` ),
                               Au2O = sweep( mineral.model.normalized, 2, unlist( min.comps["Au2O", ] ), `*` ),
                               Ag2O = sweep( mineral.model.normalized, 2, unlist( min.comps["Ag2O", ] ), `*` ),
                               As2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps["As2O3", ] ), `*` ),
                               SnO2 = sweep( mineral.model.normalized, 2, unlist( min.comps["SnO2", ] ), `*` ),
                               Sc2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps["Sc2O3", ] ), `*` ),
                               UO2 = sweep( mineral.model.normalized, 2, unlist( min.comps["UO2", ] ), `*` ),
                               V2O5 = sweep( mineral.model.normalized, 2, unlist( min.comps["V2O5", ] ), `*` ),
                               WO3 = sweep( mineral.model.normalized, 2, unlist( min.comps["WO3", ] ), `*` ),
                               ZnO = sweep( mineral.model.normalized, 2, unlist( min.comps["ZnO", ] ), `*` ),
                               SrO = sweep( mineral.model.normalized, 2, unlist( min.comps["SrO", ] ), `*` ),
                               BeO = sweep( mineral.model.normalized, 2, unlist( min.comps["BeO", ] ), `*` ),
                               PbO = sweep( mineral.model.normalized, 2, unlist( min.comps["PbO", ] ), `*` ),
                               Rb2O = sweep( mineral.model.normalized, 2, unlist( min.comps["Rb2O", ] ), `*` ),
                               Cs2O = sweep( mineral.model.normalized, 2, unlist( min.comps["Cs2O", ] ), `*` ),
                               Nb2O5 = sweep( mineral.model.normalized, 2, unlist( min.comps["Nb2O5", ] ), `*` ),
                               B2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps["B2O3", ] ), `*` )
  )

  
  minerals.comp.model$Th <- minerals.comp.model$ThO2 * 0.878809 * 10000
  minerals.comp.model$Zr <- minerals.comp.model$ZrO2 * 0.740318 * 10000
  minerals.comp.model$Ba <- minerals.comp.model$BaO * 0.895651 * 10000
  minerals.comp.model$Li <- minerals.comp.model$Li2O * (1/2.1527) * 10000
  minerals.comp.model$La <- minerals.comp.model$La2O3 * (1/1.1728) * 10000
  minerals.comp.model$Ce <- minerals.comp.model$Ce2O3 * (1/1.1713) * 10000
  minerals.comp.model$Pr <- minerals.comp.model$Pr2O3 * (1/1.1703) * 10000
  minerals.comp.model$Nd <- minerals.comp.model$Nd2O3 * (1/1.1664) * 10000
  minerals.comp.model$Sm <- minerals.comp.model$Sm2O3 * (1/1.1596) * 10000
  minerals.comp.model$Eu <- minerals.comp.model$Eu2O3 * (1/1.1579) * 10000
  minerals.comp.model$Gd <- minerals.comp.model$Gd2O3 * (1/1.1526) * 10000
  minerals.comp.model$Tb <- minerals.comp.model$Tb2O3 * (1/1.1510) * 10000
  minerals.comp.model$Dy <- minerals.comp.model$Dy2O3 * (1/1.1477) * 10000
  minerals.comp.model$Ho <- minerals.comp.model$Ho2O3 * (1/1.1455) * 10000
  minerals.comp.model$Er <- minerals.comp.model$Er2O3 * (1/1.1435) * 10000
  minerals.comp.model$Tm <- minerals.comp.model$Tm2O3 * (1/1.1421) * 10000
  minerals.comp.model$Yb <- minerals.comp.model$Yb2O3 * (1/1.1387) * 10000
  minerals.comp.model$Lu <- minerals.comp.model$Lu2O3 * (1/1.1371) * 10000
  minerals.comp.model$Y <- minerals.comp.model$Y2O3 * (1/1.2699) * 10000
  minerals.comp.model$Sb <- minerals.comp.model$Sb2O3 * (1/1.3284) * 10000
  minerals.comp.model$Cu <- minerals.comp.model$CuO * (1/1.2518) * 10000
  minerals.comp.model$Au <- minerals.comp.model$Au2O * (1/1.0406) * 10000
  minerals.comp.model$Ag <- minerals.comp.model$Ag2O * (1/1.0741) * 10000
  minerals.comp.model$As <- minerals.comp.model$As2O3 * (1/1.3203) * 10000
  minerals.comp.model$Sn <- minerals.comp.model$SnO2 * (1/1.2696) * 10000
  minerals.comp.model$Sc <- minerals.comp.model$Sc2O3 * (1/1.5338) * 10000
  minerals.comp.model$U <- minerals.comp.model$UO2 * (1/1.1344) * 10000
  minerals.comp.model$V <- minerals.comp.model$V2O5 * (1/1.7852) * 10000
  minerals.comp.model$W <- minerals.comp.model$WO3 * (1/1.2610) * 10000
  minerals.comp.model$Zn <- minerals.comp.model$ZnO * (1/1.2448) * 10000
  minerals.comp.model$Sr <- minerals.comp.model$SrO * (1/1.1826) * 10000
  minerals.comp.model$Be <- minerals.comp.model$BeO * (1/2.7758) * 10000
  minerals.comp.model$Pb <- minerals.comp.model$PbO * (1/1.0772) * 10000
  minerals.comp.model$Rb <- minerals.comp.model$Rb2O * (1/1.0936) * 10000
  minerals.comp.model$Cs <- minerals.comp.model$Cs2O * (1/1.0602) * 10000
  minerals.comp.model$Nb <- minerals.comp.model$Nb2O5 * (1/1.4305) * 10000
  minerals.comp.model$B <- minerals.comp.model$B2O3 * (1/3.2202) * 10000
  minerals.comp.model$FeOt <- ( minerals.comp.model$Fe2O3 * 0.89981 ) + minerals.comp.model$FeO 
  
  
  
  
  # sum the rows in each list to get rock oxide compositions
  wr.comp.model <- lapply( minerals.comp.model, rowSums, na.rm = T )
  # unlist them into a dataframe
  wr.comp.model.comb <- bind_rows( wr.comp.model )
  
  ## summarize the grain numbers
  grain.number.summary <- data.frame( means = apply( grain.num.model, 2, mean ),
                                      stdev = apply( grain.num.model, 2, sd ) )
  grain.number.summary$rsd <- grain.number.summary$stdev / grain.number.summary$means * 100
  # grain.number.summary
  ## summarize the grain masses
  grain.mass.summary <- data.frame( means = apply( mineral.model.normalized, 2, mean ),
                                    stdev = apply( mineral.model.normalized, 2, sd ) )
  grain.mass.summary$rsd <- grain.mass.summary$stdev / grain.mass.summary$means * 100
  # grain.mass.summary
  ## summarize the wr.results
  wr.summary.table <- data.frame( means = apply( wr.comp.model.comb, 2, mean ),
                                  medians = apply( wr.comp.model.comb, 2, median ),
                                  stdev = apply( wr.comp.model.comb, 2, sd ) )
  wr.summary.table$rsd <- wr.summary.table$stdev / wr.summary.table$means * 100
  
  
  ## cbind into a large single data frame
  model.rock.data <- cbind( mineral.model.normalized, wr.comp.model.comb )
  
  # wr.summary.table
  assign( paste( "mineral.model.normalized.", sample.mass, ".g", sep="" ), mineral.model.normalized )
  assign( paste( "wr.comp.model.comb.", sample.mass, ".g", sep="" ), wr.comp.model.comb )
  assign( paste( "wr.summary.table.", sample.mass, ".g", sep="" ), wr.summary.table )
  assign( paste( "grain.mass.summary.", sample.mass, ".g", sep="" ), grain.mass.summary )
  assign( paste( "grain.number.summary.", sample.mass, ".g", sep="" ), grain.number.summary )
  assign( paste( "model.rock.data.", sample.mass, ".g", sep="" ), model.rock.data )
  
}


## compare the RSDs of the stdev for the mineral grain number used as an input to the 
#.  output normalized RSDs of the grain numbers
min.rsd.comp <- grain.mass.summary
min.rsd.comp$input.grain.rsd <- t( rsd.err.num.min.grains.h )
min.rsd.comp$diff <- min.rsd.comp$rsd - min.rsd.comp$input.grain.rsd 
## plot the mineral grain number RSDs change between input and modeled against the mineral mode
ggplot( min.rsd.comp, aes( x= means, y = diff, label = row.names(min.rsd.comp) ) ) +
  fte_theme_white() +
  # ylim( 3, -3) +
  # scale_x_log10() +
  labs( y = "RSD change",
        x = "Mode (fraction)") +
  geom_point() + geom_text( hjust=-0.2, vjust=-0.2 )



## combine RSDs from the datasets
summary.raw.sds <- data.frame( mass = sample.masses,
                               rsds = rbind( t( wr.summary.table.625.g$stdev ),
                                             t( wr.summary.table.1250.g$stdev ),
                                             t( wr.summary.table.2500.g$stdev ),
                                             t( wr.summary.table.5000.g$stdev ),
                                             t( wr.summary.table.10000.g$stdev ),
                                             t( wr.summary.table.20000.g$stdev ) ) )
colnames(summary.raw.sds) <- c( "mass", rownames(wr.summary.table.625.g))
summary.raw.sds <- data.frame(lapply(summary.raw.sds, function(x) round(x, digits = 3)))
write.csv( summary.raw.sds, "Run.sds.csv")

summary.rsds <- data.frame( mass = sample.masses,
                            rsds = rbind( t( wr.summary.table.625.g$rsd ),
                                          t( wr.summary.table.1250.g$rsd ),
                                          t( wr.summary.table.2500.g$rsd ),
                                          t( wr.summary.table.5000.g$rsd ),
                                          t( wr.summary.table.10000.g$rsd ),
                                          t( wr.summary.table.20000.g$rsd )) )
summary.rsds <- data.frame(lapply(summary.rsds, function(x) round(x, digits = 3)))

colnames(summary.rsds) <- c( "mass", rownames(wr.summary.table.625.g))
write.csv( summary.rsds, "Run.rsds.csv")


summary.means <- data.frame( mass = sample.masses,
                             means = rbind( t( wr.summary.table.625.g$means ),
                                            t( wr.summary.table.1250.g$means ),
                                            t( wr.summary.table.2500.g$means ),
                                            t( wr.summary.table.5000.g$means ),
                                            t( wr.summary.table.10000.g$means ),
                                            t( wr.summary.table.20000.g$means )) )
summary.means <- data.frame(lapply(summary.means, function(x) round(x, digits = 3)))

colnames(summary.means) <- c( "mass", rownames(wr.summary.table.625.g))
write.csv( summary.means, "Run.means.csv")



## summarize grain mass RSDS and grain number RSDs
summary.grain.mass.rsds <- data.frame( mass = sample.masses,
                                       rsds = rbind( t( grain.mass.summary.625.g$rsd ),
                                                     t( grain.mass.summary.1250.g$rsd ),
                                                     t( grain.mass.summary.2500.g$rsd ),
                                                     t( grain.mass.summary.5000.g$rsd ),
                                                     t( grain.mass.summary.10000.g$rsd ),
                                                     t( grain.mass.summary.20000.g$rsd ) ) )
colnames(summary.grain.mass.rsds) <- c( "mass", rownames(grain.mass.summary.625.g))
summary.grain.number.rsds <- data.frame( mass = sample.masses,
                                         rsds = rbind( t( grain.number.summary.625.g$rsd ),
                                                       t( grain.number.summary.1250.g$rsd ),
                                                       t( grain.number.summary.2500.g$rsd ),
                                                       t( grain.number.summary.5000.g$rsd ),
                                                       t( grain.number.summary.10000.g$rsd ),
                                                       t( grain.number.summary.20000.g$rsd ) ) )
colnames(summary.grain.number.rsds) <- c( "mass", rownames(grain.number.summary.625.g ) )

## calculate mean rock comps
summary.wr.comp.mean <- data.frame( mass = sample.masses,
                                    rsds = rbind( t( wr.summary.table.625.g$means ),
                                                  t( wr.summary.table.1250.g$means ),
                                                  t( wr.summary.table.2500.g$means ),
                                                  t( wr.summary.table.5000.g$means ),
                                                  t( wr.summary.table.10000.g$means ),
                                                  t( wr.summary.table.20000.g$means ) ) )
colnames(summary.wr.comp.mean) <- c( "mass", rownames( wr.summary.table.10000.g ) )
write.csv( summary.wr.comp.mean, "Run.mean.comps.csv" )


## read in the selected columns from the User input
oxide.order <- read.csv( "oxideschosen.csv" )
trace.order <- read.csv( "traceschosen.csv")
# transpose to get column names
oxide.order.plot <- t(oxide.order)
colnames(oxide.order.plot) <- oxide.order.plot[1,]

trace.order.plot <- t(trace.order)
colnames(trace.order.plot) <- trace.order.plot[1,]
full.plot.order <-  c(colnames(oxide.order.plot), colnames(trace.order.plot))


### define the order of the oxides for plotting and exporting 
oxide.order.rsdplot <-  colnames( oxide.order.plot )




summary.rsds.test <- summary.rsds[, c("mass", oxide.order.rsdplot ) ]


### plot  RSDs for the oxides
raw.rsd.plot <- summary.rsds.test %>% 
  gather( var, val , -mass, factor_key = TRUE ) %>% 
  # select( colnames( summary.raw.sds[ , 1:10] ) ) %>%
  ggplot( aes( x = var, y = val, color = as.factor( mass ), group = as.factor( mass ) ) ) +
  geom_line( linewidth = 1.5 ) +
  # ylim( 0, 10 ) +
  fte_theme_white()+
  theme( axis.text.x = element_text( size = 14 ),
         axis.title.x = element_blank() ) +
  labs( color = "Mass (g)",
        y = "RSD on Oxide (%)")
raw.rsd.plot

data_long <- reshape2::melt( summary.rsds.test, id.vars = "mass", variable.name = "oxide", value.name = "value")

# Plot the data using ggplot2
mass.oxide.plot <-  ggplot( data_long, aes( x = mass, y = value, color = oxide, group = oxide ) ) +
  fte_theme_white() +
  geom_line( linewidth = 1.5 ) +
  labs(x = "Mass", y = "Relative SD", color = "Variable") 
# scale_color_viridis( discrete = T ) 
mass.oxide.plot





data_long_means <- reshape2::melt( summary.means[ , c("mass", full.plot.order)], id.vars = "mass", variable.name = "oxide", value.name = "mean" )
data_long_stdev <- reshape2::melt( summary.raw.sds[ , c("mass", full.plot.order)], id.vars = "mass", variable.name = "oxide", value.name = "stdevs" )


merged_df <- merge(data_long_means, data_long_stdev, by = c("mass", "oxide"))
merged_df <- merged_df[order(merged_df$mass), ]
merged_df$oxide <- factor(merged_df$oxide, levels = full.plot.order)
merged_df <- merged_df[order(merged_df$oxide), ]
# merged_df$dodge <- rep(c(2,4,6,8,10,13), each = 4)


## table for exporting final data
wr.summary.table.export <- NULL
for( i in 1:length(sample.masses) ) {
  wr.summary.table.export <- rbind( wr.summary.table.export, 
                                    summary.means[ i, ], summary.raw.sds[i,], summary.rsds[i, ] )
  
}
wr.summary.table.export <- wr.summary.table.export[, c("mass", full.plot.order)]
wr.summary.table.export$value <- rep( c("Mean", "StDev", "RSD" ), times = length(sample.masses) )
wr.summary.table.export <- wr.summary.table.export[, c(ncol(wr.summary.table.export), 1:(ncol(wr.summary.table.export)-1))]



raw.sd.plot <- ggplot(merged_df, aes( x = mass, y = mean ) ) +
  fte_theme_white() +
  geom_point( color = "#645153", size = 4, position = position_dodge2(width = 1.5 ) ) +
  geom_linerange(aes(ymin = mean - stdevs, ymax = mean + stdevs),
                 position = position_dodge2(width = 1.5), color = "#645153", size = 0.6 ) +
  # viridis::scale_color_viridis() +
  # geom_errorbar(aes( color = mass ), linewidth = 1,width=0, position = position_jitterdodge(dodge.width = 2, jitter.width = 2, seed = 1)) +
  theme( axis.text.x = element_text( size = 10 ),
         axis.text.y = element_text( size = 10 ),
         axis.title.x = element_text( size = 20 ),
         axis.title.y = element_text( size = 20 ),
         strip.background = element_rect(
           color="black", fill="#213e47", size=0.2, linetype="solid"
         ),
         strip.text.x = element_text(
           size = 12, color = "white" ) ) +
  coord_cartesian() +
  labs( x = "Sample Mass (g)", y = "Wt% Oxide") +
  # expand_limits(y = 1) +
  facet_wrap(oxide ~ ., scales = "free" ) +
  theme(aspect.ratio = 1)

raw.sd.plot


# plot_data
## create fake data to plot
# fake.plot.data <- with(merged_df,
#                        data.frame( mean=c(mean+stdevs+0.05,mean-stdevs-0.05),
#                                    stdevs=c(stdevs,stdevs),
#                                    oxide=c(oxide,oxide),
#                                    mass = c(mass, mass) ) ) 
# 
# raw.sd.plot <- plot_data + geom_point(data=fake.plot.data,x=NA)
# raw.sd.plot



