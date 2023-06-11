# Modeling U-Pb data for resolvability of U-Pb discordance modeling ages


############################# CLEAR THE ENVIRONMENT, RUN EVERY TIME! #############################
rm( list = ls( ) )


computer       <- "JesseReimink"
if( computer == "jxr1350" ) {
  google.drive.dir <- paste( "/Users/", computer, 
                             "/My Drive/", sep = "" )
} else {
  google.drive.dir <- paste( "/Users/", computer, 
                             "/Google Drive (jxr1350@psu.edu)/", sep = "" )
}
#############################  SET THE DIRECTORIES #############################
source.dir  <- paste( google.drive.dir, 
                      "Research/PostDoc/R scripts/UPb modeling Working", sep = "" )
lib.dir  <- paste( google.drive.dir, 
                   "Research/PostDoc/R scripts/", sep = "" )
data.dir		<- paste( google.drive.dir, 
                    "Research/PSU/Projects/Treatise on Geochemistry Chapter/Modeling Uncertainties", sep = "" )
export.dir		<- paste( google.drive.dir, 
                      "Research/PSU/Projects/Basin Fluid Flow/Renan Thesis/Modeling Sensitivity", sep = "" )
setwd( lib.dir )
source( "Libraries.R" )
source( "Constants.R" )
source( "Equations.R" )
source( "fte_theme_plotting.R" )
source( "GCDkit functions.R" )
setwd( data.dir )

#############################  DATA IMPORT #############################
# Import the modes, densities, and molar masses

modes <- read.csv( "min_modes.csv")
densities <- read.csv( "min_density.csv" )
mol.masses <- read.csv( "min_weights.csv" )
grain.dims <- read.csv( "min_dimensions.csv" )
min.dists <- read.csv( "min_distributions.csv" )
min.comps.raw <- read.csv( "min_comps.csv" ) # this has rows arranged in order of the Oxide Weight Percent Composition in the Stanley spreadsheet
min.comps <- min.comps.raw[ , -1 ]

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
  
  # iterate until the SDs match the Stanley one - accounting for internal 
  # normalization processes
# 
#   
#   tol   <- 0.1
#   loop.number <- 1
# 
#   stanley.rsds <- rsd.err.num.min.grains
#   original.norm.mass.rsds <- grain.mass.summary$rsd 
#   reset.norm.mass.rsds <- original.norm.mass.rsds 
#   reset.err.min.grains <- err.num.min.grains
#   # new.grain.num.err <- err.num.min.grains
#   
#   repeat {
#     reset.norm.mass.rsds <- grain.mass.summary$rsd
#     difs <- reset.norm.mass.rsds - stanley.rsds
#     print( loop.number )
#     print( sum( abs( difs ), na.rm = T ) )
#     if ( sum( abs( difs ), na.rm = T ) < tol ) {
#       break
#     }
#     reset.err.min.grains <- ( -1 * difs / 100 ) * reset.err.min.grains + reset.err.min.grains
# 
#     for (col in 1:ncol( grain.num.model ) ) {
#       # Extract mean and sd for current column
#       mean_val <- num.min.grains[1, col]
#       sd_val <- reset.err.min.grains[1, col]
#       # Generate resampled value based on mean and sd
#       resampled_val <- rnorm( nreps, mean = mean_val, sd = sd_val )
#       # Assign resampled value to corresponding column in the empty dataframe
#       grain.num.model[, col] <- resampled_val
#     }
#     ## convert number of grains to mass of each grain
#     grain.mass.model <- sweep( grain.num.model, 2, unlist( grain.mass[1, ] ), `*` )
# 
#     ## calculate the wt% of each mineral
#     # Calculate the column totals
#     mineral.model.tot.mass <- rowSums( grain.mass.model, na.rm = T )
# 
#     # Normalize the columns and add them to the dataframe
#     mineral.model.normalized <- grain.mass.model %>%
#       mutate( across( everything(), ~ . / mineral.model.tot.mass ) )
# 
#     ## summarize the grain masses
#     grain.mass.summary <- data.frame( means = apply( mineral.model.normalized, 2, mean ),
#                                       stdev = apply( mineral.model.normalized, 2, sd ) )
#     grain.mass.summary$rsd <- grain.mass.summary$stdev / grain.mass.summary$means * 100
# 
#     loop.number <- loop.number + 1
#     assign( paste( "grain.mass.summary.", loop.number, ".loop", sep="" ), grain.mass.summary )
#   }
#   
#   
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
                                  stdev = apply( wr.comp.model.comb, 2, sd ) )
  wr.summary.table$rsd <- wr.summary.table$stdev / wr.summary.table$means * 100
  # wr.summary.table
  assign( paste( "mineral.model.normalized.", sample.mass, ".g", sep="" ), mineral.model.normalized )
  assign( paste( "wr.comp.model.comb.", sample.mass, ".g", sep="" ), wr.comp.model.comb )
  assign( paste( "wr.summary.table.", sample.mass, ".g", sep="" ), wr.summary.table )
  assign( paste( "grain.mass.summary.", sample.mass, ".g", sep="" ), grain.mass.summary )
  assign( paste( "grain.number.summary.", sample.mass, ".g", sep="" ), grain.number.summary )
  
}


## compare the RSDs of the stdev for the mineral grain number used as an input to the 
#.  output normalized RSDs of the grain numbers
min.rsd.comp <- grain.mass.summary
min.rsd.comp$input.grain.rsd <- t( rsd.err.num.min.grains.h )
min.rsd.comp$diff <- min.rsd.comp$input.grain.rsd - min.rsd.comp$rsd 
## plot the mineral grain number RSDs change between input and modeled against the mineral mode
ggplot( min.rsd.comp, aes( x= means, y = diff, label = row.names(min.rsd.comp) ) ) +
  fte_theme_white() +
  ylim( 3, -3) +
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



## read in observed SB granite RSDs
sb.granite.observed.rsds <- read.csv( "SB Granite Observed RSDs.csv", header = T )
sb.granite.observed.sds <- read.csv( "SB Granite Observed SDs.csv", header = T )
sb.granite.predicted.rsds <- read.csv( "SB Granite Predicted RSDs.csv", header = T )
sb.granite.predicted.sds <- read.csv( "SB Granite Predicted SDs.csv", header = T )

## compare 625g sample to the observed
comparison.625g <- data.frame( rbind( sb.granite.observed.rsds[ 3, ],
                                      sb.granite.predicted.rsds[ 3, ],
                                      summary.rsds[ 3, c( "mass", "SiO2", "TiO2", "Al2O3", "FeO", "MnO", "MgO", "CaO", "Na2O", "K2O", "P2O5")] ) )
comparison.625g$category <- c( "Observed", "Stanley", "Predicted" )

comparison.625g[, -1] %>% 
  gather( var, val , -category, factor_key = TRUE ) %>% 
  # select( colnames( summary.raw.sds[ , 1:10] ) ) %>%
  ggplot( aes( x = var, y = val, color = as.factor( category ), group = as.factor( category ) ) ) +
  geom_line( linewidth = 1.5 ) +
  fte_theme_white()+
  theme( axis.text.x = element_text( size = 10 ),
         axis.title.x = element_blank() ) +
  labs( color = "Category",
        y = "RSD on Oxide (wt%)") 



### plot raw SDs for the oxides
summary.raw.sds[ , c(1:4,7:12,15)] %>% 
  gather( var, val , -mass, factor_key = TRUE ) %>% 
  # select( colnames( summary.raw.sds[ , 1:10] ) ) %>%
  ggplot( aes( x = var, y = val, color = as.factor( mass ), group = as.factor( mass ) ) ) +
  geom_line( linewidth = 1.5 ) +
  fte_theme_white()+
  theme( axis.text.x = element_text( size = 10 ),
         axis.title.x = element_blank() ) +
  labs( color = "Mass (g)",
        y = "SD on Oxide (wt%)") 



### plot raw SDs for the oxides
summary.raw.sds[ , c(1,20)] %>% 
  ggplot( aes( x = LiO2, y = mass / 1000 ) ) +
  geom_line( linewidth = 1.5 ) +
  fte_theme_white()+
  labs( y = "Mass (kg)",
        x = "SD on LiO2 (wt%)")


### plot  RSDs for the oxides
summary.rsds[ , c(1:4,7:12)] %>% 
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


ggplot( wr.summary.table.625.g, aes( x = sd ) ) 
  












