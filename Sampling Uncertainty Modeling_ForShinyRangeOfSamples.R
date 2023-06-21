# Modeling U-Pb data for resolvability of U-Pb discordance modeling ages


############################# CLEAR THE ENVIRONMENT, RUN EVERY TIME! #############################
# rm( list = ls( ) )



fte_theme_white <- function() {
  
  library( RColorBrewer )
  # Generate the colors for the chart procedurally with RColorBrewer
  palette          <- brewer.pal( "Greys", n = 9 )
  # color.background = palette[ 2 ]
  color.grid.major = palette[ 4 ]
  color.axis.text  = palette[ 7 ]
  color.axis.title = palette[ 8 ]
  color.title      = palette[ 9 ]
  
  # Begin construction of chart
  theme_bw( base_size = 9 ) +
    
    # Set the entire chart region to a light gray color
    # theme( panel.background = element_rect( fill = color.background, color = color.background ) ) +
    # theme( plot.background  = element_rect( fill = color.background, color = color.background ) ) +
    # theme( panel.border     = element_rect( color = color.background ) ) +
    
    # Format the grid
    theme( panel.grid.major = element_line( color = color.grid.major, size = 0.25 ) ) +
    theme( axis.line.x      = element_line( color = color.grid.major, size = 0.35 ) ) +
    theme( axis.line.y      = element_line( color = color.grid.major, size = 0.35 ) ) +
    theme( panel.grid.minor = element_blank( ) ) +
    theme( axis.ticks       = element_blank( ) ) +
    # theme( legend.key = element_rect( fill = color.background ) ) +
    
    
    # Format the legend
    # theme( legend.position = "none" ) +
    # theme( legend.background = element_rect( fill = color.background ) ) +
    theme( legend.text       = element_text( size = 7, color = color.axis.title ) ) +
    
    # Set title and axis labels, and format these and tick marks
    theme( plot.title   = element_text( color = color.title, size = 25, vjust = 1.25, hjust = 0.5 ) ) +
    theme( plot.subtitle   = element_text( color = color.title, size = 16, vjust = 1.25, hjust = 0.5 ) ) +
    theme( axis.text.x  = element_text( size = 28, color = color.axis.text ) ) +
    theme( axis.text.y  = element_text( size = 28, color = color.axis.text ) ) +
    theme( axis.title.x = element_text( size = 36, color = color.axis.title, vjust = 0 ) ) +
    theme( axis.title.y = element_text( size = 36, color = color.axis.title, vjust = 1.25 ) ) +
    
    # Plot margins
    theme( plot.margin = unit( c( 0.35, 0.2, 0.3, 0.35 ), "cm" ) )
}

shiny.min.parameters <- read.csv("/Users/jessereimink/Google Drive (jxr1350@psu.edu)/Research/PSU/Projects/Treatise on Geochemistry Chapter/Modeling Uncertainties/Shiny_OutPut.csv", stringsAsFactors = FALSE)

  
#############################  DATA IMPORT #############################
# Import the modes, densities, and molar masses
modes <- t( shiny.min.parameters$Modes )
densities <- t( shiny.min.parameters$Densities )
mol.masses <- t( shiny.min.parameters$Weights )
grain.dims <- t( cbind( shiny.min.parameters$Dimension1, shiny.min.parameters$Dimension2, shiny.min.parameters$Dimension3 )  )
min.dists <- t( shiny.min.parameters$Distribution )
min.comps <- t( cbind( shiny.min.parameters[, c(9:ncol(shiny.min.parameters))] ) )

## Set column names to be minerals
colnames(modes) = shiny.min.parameters$X
colnames(densities) = shiny.min.parameters$X
colnames(mol.masses) = shiny.min.parameters$X
colnames(grain.dims) = shiny.min.parameters$X
colnames(min.dists) = shiny.min.parameters$X
colnames(min.comps) = shiny.min.parameters$X



## set sample mass
sample.masses <- c( 625, 1250, 2500, 5000, 10000, 20000, 40000 ) #in grams
N_n <- 1000


for( j in 1:length(sample.masses)) {
  # set sample mass
  
  sample.mass = sample.masses[j]
  ## this all replicates the first columns in Stanley spreadsheet
  inverse.molar.volume <- 1 / ( mol.masses / densities ) 
  modes.internal <- modes * densities 
  modes.wt.total <- rowSums( modes.internal )
  modes.wt <- modes.internal / modes.wt.total * 100 
  modes.moles.internal <- modes * inverse.molar.volume 
  tot.modes.moles.internal <- rowSums( modes.moles.internal )
  modes.mole <- modes.moles.internal / tot.modes.moles.internal * 100 
  grain.vols <- apply( grain.dims, 2, prod )
  grain.mass <- grain.vols * densities
  tot.min.mass <- sample.mass * modes.wt / 100 
  num.min.grains <- tot.min.mass / grain.mass
  tot.num.min.grains <- num.min.grains/( modes/100 )
  ## calculate error, using the hypergeometric calc in Stanley spreadsheet

  err.num.min.grains.h <- sqrt(( 100-modes)/100*modes/100*tot.num.min.grains*
                               (N_n*tot.num.min.grains-tot.num.min.grains)/(N_n*tot.num.min.grains-1))
  err.num.min.grains.b <- sqrt(( 100-modes)/100*modes/100*tot.num.min.grains)
  err.num.min.grains.p <- sqrt(tot.num.min.grains)
  
  rsd.err.num.min.grains.h <- err.num.min.grains.h / num.min.grains * 100
  rsd.err.num.min.grains.b <- err.num.min.grains.b / num.min.grains * 100
  rsd.err.num.min.grains.p <- err.num.min.grains.p / num.min.grains * 100
  
  ## create a random set of mineral grains in a given rock, then recalculate the mode in wt% of each grain
  nreps = 10000
  grain.num.model <- data.frame( matrix( nrow = nreps, ncol = ncol( num.min.grains ) ) )
  colnames( grain.num.model ) <- colnames( num.min.grains )
  
  
  
  
## model the number of mineral grains in a synthetic rock
  for (col in 1:ncol( grain.num.model ) ) {
    # Extract mean and sd for current column
    mean_val <- num.min.grains[ 1, col ]
    # calculate SD based on the distribution of the grains
    sd_val <- ifelse( min.dists[col] == "h", err.num.min.grains.h[1, col],
                      ifelse( min.dists[col] == "b", err.num.min.grains.b[1, col],
                              err.num.min.grains.p[1, col] ) )

    # Generate resampled value based on mean and sd
    resampled_val <- rnorm( nreps, mean = mean_val, sd = sd_val )
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
  minerals.comp.model <- list( SiO2 = sweep( mineral.model.normalized, 2, unlist( min.comps[1, ] ), `*` ),
                               TiO2 = sweep( mineral.model.normalized, 2, unlist( min.comps[2, ] ), `*` ),
                               Al2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps[3, ] ), `*` ),
                               Cr2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps[4, ] ), `*` ),
                               Fe2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps[5, ] ), `*` ),
                               FeO = sweep( mineral.model.normalized, 2, unlist( min.comps[6, ] ), `*` ),
                               MnO = sweep( mineral.model.normalized, 2, unlist( min.comps[7, ] ), `*` ),
                               MgO = sweep( mineral.model.normalized, 2, unlist( min.comps[8, ] ), `*` ),
                               CaO = sweep( mineral.model.normalized, 2, unlist( min.comps[9, ] ), `*` ),
                               Na2O = sweep( mineral.model.normalized, 2, unlist( min.comps[10, ] ), `*` ),
                               K2O = sweep( mineral.model.normalized, 2, unlist( min.comps[11, ] ), `*` ),
                               O = sweep( mineral.model.normalized, 2, unlist( min.comps[12, ] ), `*` ),
                               H2O = sweep( mineral.model.normalized, 2, unlist( min.comps[13, ] ), `*` ),
                               P2O5 = sweep( mineral.model.normalized, 2, unlist( min.comps[14, ] ), `*` ),
                               ZrO2 = sweep( mineral.model.normalized, 2, unlist( min.comps[15, ] ), `*` ),
                               Ce2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps[16, ] ), `*` ),
                               CO2 = sweep( mineral.model.normalized, 2, unlist( min.comps[17, ] ), `*` ),
                               SO2 = sweep( mineral.model.normalized, 2, unlist( min.comps[18, ] ), `*` ),
                               LiO2 = sweep( mineral.model.normalized, 2, unlist( min.comps[19, ] ), `*` ),
                               ThO2 = sweep( mineral.model.normalized, 2, unlist( min.comps[20, ] ), `*` )
  )
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
                                          t( wr.summary.table.20000.g$stdev ),
                                          t( wr.summary.table.40000.g$stdev ) ) )
colnames(summary.raw.sds) <- c( "mass", rownames(wr.summary.table.625.g))
write.csv( summary.raw.sds, "Run.sds.csv" )

summary.rsds <- data.frame( mass = sample.masses,
                            rsds = rbind( t( wr.summary.table.625.g$rsd ),
                                          t( wr.summary.table.1250.g$rsd ),
                                          t( wr.summary.table.2500.g$rsd ),
                                          t( wr.summary.table.5000.g$rsd ),
                                          t( wr.summary.table.10000.g$rsd ),
                                          t( wr.summary.table.20000.g$rsd ),
                                          t( wr.summary.table.40000.g$rsd )) )
colnames(summary.rsds) <- c( "mass", rownames(wr.summary.table.625.g))
write.csv( summary.rsds, "Run.rsds.csv" )


## summarize grain mass RSDS and grain number RSDs
summary.grain.mass.rsds <- data.frame( mass = sample.masses,
                            rsds = rbind( t( grain.mass.summary.625.g$rsd ),
                                          t( grain.mass.summary.1250.g$rsd ),
                                          t( grain.mass.summary.2500.g$rsd ),
                                          t( grain.mass.summary.5000.g$rsd ),
                                          t( grain.mass.summary.10000.g$rsd ),
                                          t( grain.mass.summary.20000.g$rsd ),
                                          t( grain.mass.summary.40000.g$rsd ) ) )
colnames(summary.grain.mass.rsds) <- c( "mass", rownames(grain.mass.summary.625.g))
summary.grain.number.rsds <- data.frame( mass = sample.masses,
                                       rsds = rbind( t( grain.number.summary.625.g$rsd ),
                                                     t( grain.number.summary.1250.g$rsd ),
                                                     t( grain.number.summary.2500.g$rsd ),
                                                     t( grain.number.summary.5000.g$rsd ),
                                                     t( grain.number.summary.10000.g$rsd ),
                                                     t( grain.number.summary.20000.g$rsd ),
                                                     t( grain.number.summary.40000.g$rsd ) ) )
colnames(summary.grain.number.rsds) <- c( "mass", rownames(grain.number.summary.625.g ) )

## calculate mean rock comps
summary.wr.comp.mean <- data.frame( mass = sample.masses,
                                         rsds = rbind( t( wr.summary.table.625.g$means ),
                                                       t( wr.summary.table.1250.g$means ),
                                                       t( wr.summary.table.2500.g$means ),
                                                       t( wr.summary.table.5000.g$means ),
                                                       t( wr.summary.table.10000.g$means ),
                                                       t( wr.summary.table.20000.g$means ),
                                                       t( wr.summary.table.40000.g$means ) ) )
colnames(summary.wr.comp.mean) <- c( "mass", rownames( wr.summary.table.10000.g ) )
write.csv( summary.wr.comp.mean, "Run.mean.comps.csv" )




### plot raw SDs for the oxides
raw.sd.plot <- summary.raw.sds[ , c(1:4,7:12,15)] %>% 
  gather( var, val , -mass, factor_key = TRUE ) %>% 
  # select( colnames( summary.raw.sds[ , 1:10] ) ) %>%
  ggplot( aes( x = var, y = val, color = as.factor( mass ), group = as.factor( mass ) ) ) +
  geom_line( linewidth = 1.5 ) +
  fte_theme_white()+
  theme( axis.text.x = element_text( size = 10 ),
         axis.title.x = element_blank() ) +
  labs( color = "Mass (g)",
        y = "SD on Oxide (wt%)") 


## plotting compositions vs the mineral modes
ggplot( model.rock.data, aes( x = Al2O3, y = FeO, color = bt ) ) +
  geom_point(alpha = 0.2 ) +
  fte_theme_white() +
  scale_color_viridis()


### plot raw SDs for the oxides
summary.raw.sds[ , c(1,20)] %>% 
  ggplot( aes( x = LiO2, y = mass / 1000 ) ) +
  geom_line( linewidth = 1.5 ) +
  fte_theme_white()+
  labs( y = "Mass (kg)",
        x = "SD on LiO2 (wt%)")


### plot  RSDs for the oxides
raw.rsd.plot <- summary.rsds[ , c(1:4,7:12)] %>% 
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


data_long <- reshape2::melt( summary.rsds[ , c(1:4,7:12)], id.vars = "mass", variable.name = "oxide", value.name = "value")

# Plot the data using ggplot2
mass.oxide.plot <-  ggplot( data_long, aes( x = mass, y = value, color = oxide, group = oxide ) ) +
  fte_theme_white() +
  geom_line( linewidth = 1.5 ) +
  labs(x = "Mass", y = "Relative SD", color = "Variable") 
  # scale_color_viridis( discrete = T ) 

summary.rsds[ , c(1:4,7:12)] %>% 
  # select( colnames( summary.raw.sds[ , 1:10] ) ) %>%
  ggplot( aes( x = mass, y = SiO2 ) ) +
  geom_line( linewidth = 1.5 ) +
  # ylim( 0, 10 ) +
  fte_theme_white()+
  theme( axis.text.x = element_text( size = 14 ),
         axis.title.x = element_blank() ) +
  labs( color = "Mass (g)",
        y = "RSD on Oxide (%)")




### plot  RSDs on the mineral masses
summary.grain.mass.rsds[ ,c( 1:4,10:11) ] %>% 
  gather( var, val , -mass, factor_key = TRUE ) %>% 
  # select( colnames( summary.raw.sds[ , 1:10] ) ) %>%
  ggplot( aes( x = var, y = val, color = as.factor( mass ), group = as.factor( mass ) ) ) +
  geom_line( linewidth = 1.5 ) +
  ylim(0, 15) +
  fte_theme_white()+
  theme( axis.text.x = element_text( size = 14 ),
         axis.title.x = element_blank() ) +
  labs( color = "Mass (g)",
        y = "RSD on Mineral Mass")

### plot  RSDs on the mineral grain number
summary.grain.number.rsds[ ,c( 1:4,10:11) ] %>% 
  gather( var, val , -mass, factor_key = TRUE ) %>% 
  # select( colnames( summary.raw.sds[ , 1:10] ) ) %>%
  ggplot( aes( x = var, y = val, color = as.factor( mass ), group = as.factor( mass ) ) ) +
  geom_line( linewidth = 1.5 ) +
  fte_theme_white()+
  ylim(0, 15) +
  theme( axis.text.x = element_text( size = 14 ),
         axis.title.x = element_blank() ) +
  labs( color = "Mass (g)",
        y = "RSD on Mineral Grain Number")


ggplot( wr.comp.model.comb.625.g, aes( x = SiO2 ) ) +
  fte_theme_white()+
  geom_density( fill = viridis(5)[1] ) +
  geom_density( data = wr.comp.model.comb.1250.g, fill = viridis(5)[2] ) +
  geom_density( data = wr.comp.model.comb.2500.g, fill = viridis(5)[3] ) +
  geom_density( data = wr.comp.model.comb.5000.g, fill = viridis(5)[4] ) +
  geom_density( data = wr.comp.model.comb.10000.g, fill = viridis(5)[5] ) 










