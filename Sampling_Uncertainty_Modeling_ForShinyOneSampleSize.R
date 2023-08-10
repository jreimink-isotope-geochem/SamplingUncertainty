
## Single Sample Size

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


shiny.min.parameters <- read.csv( "Shiny_OutPut.csv", stringsAsFactors = FALSE)
sample.mass.g <- read.csv("sample_size.csv", stringsAsFactors = FALSE)

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
colnames(grain.dims) = shiny.min.parameters$X
colnames(min.dists) = shiny.min.parameters$X
colnames(min.comps) = shiny.min.parameters$X


N_n <- 1000


# set sample mass
sample.mass <- sample.mass.g$sample.size[1]

## this all replicates the first columns in Stanley spreadsheet
# inverse.molar.volume <- 1 / ( mol.masses / densities ) 
modes.internal <- modes * densities 
modes.wt.total <- rowSums( modes.internal, na.rm = T )
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
err.num.min.grains.p <- sqrt(num.min.grains)

rsd.err.num.min.grains.h <- err.num.min.grains.h / num.min.grains * 100
rsd.err.num.min.grains.b <- err.num.min.grains.b / num.min.grains * 100
rsd.err.num.min.grains.p <- err.num.min.grains.p / num.min.grains * 100

## create a random set of mineral grains in a given rock, then recalculate the mode in wt% of each grain
nreps = 10000
grain.num.model <- data.frame( matrix( nrow = nreps, ncol = ncol( num.min.grains ) ) )
colnames( grain.num.model ) <- colnames( num.min.grains )



## model the number of mineral grains in a synthetic rock
for ( col in 1:ncol( grain.num.model ) ) {
  # col = 3
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
# Calculate the row totals
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
                             SO2 = sweep( mineral.model.normalized, 2, unlist( min.comps["SO2", ] ), `*` ),
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
                             Y2O3 = sweep( mineral.model.normalized, 2, unlist( min.comps["Y2O3", ] ), `*` )
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

wr.summary.table <- data.frame(lapply(wr.summary.table, function(x) round(x, digits = 4)))
wr.summary.table <- t( wr.summary.table )
colnames(wr.summary.table) <- colnames(wr.comp.model.comb)


## cbind into a large single data frame
model.rock.data <- cbind( mineral.model.normalized, wr.comp.model.comb )

# wr.summary.table
assign( paste( "mineral.model.normalized.", sample.mass, ".g", sep="" ), mineral.model.normalized )
assign( paste( "wr.comp.model.comb.", sample.mass, ".g", sep="" ), wr.comp.model.comb )
assign( paste( "wr.summary.table.", sample.mass, ".g", sep="" ), wr.summary.table )
assign( paste( "grain.mass.summary.", sample.mass, ".g", sep="" ), grain.mass.summary )
assign( paste( "grain.number.summary.", sample.mass, ".g", sep="" ), grain.number.summary )
assign( paste( "model.rock.data.", sample.mass, ".g", sep="" ), model.rock.data )








## summarize grain mass RSDS and grain number RSDs
summary.grain.mass.rsds <- data.frame( mass = sample.mass,
                                       rsds = rbind( t( grain.mass.summary$rsd ) ) )
colnames(summary.grain.mass.rsds) <- c( "mass", rownames(grain.mass.summary ) )
summary.grain.number.rsds <- data.frame( mass = sample.mass,
                                         rsds = rbind( t( grain.number.summary$rsd ) ) )
colnames(summary.grain.number.rsds) <- c( "mass", rownames(grain.number.summary ) )


wr.summary.table <- as.data.frame(wr.summary.table)



## plot the RSDs of the oxide
means_to_plot <- wr.summary.table[, ]
means_long <- data.frame(
  Category = factor(names(means_to_plot), levels = names(means_to_plot)),
  Mean = as.numeric(means_to_plot[ 'means',]),
  StDev = as.numeric(means_to_plot[ 'stdev',]),
  RSD = as.numeric(means_to_plot[ 'rsd',]),
  stringsAsFactors = FALSE
)
means_long$Group <- 1


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
oxide.order.rsdplot <-  colnames(oxide.order.plot)


means_long$Category <- factor( means_long$Category, levels = oxide.order.rsdplot)
means_long <- means_long[complete.cases(means_long), ]
means_long <- means_long[order(means_long$Category), ]


raw.rsd.plot <-ggplot( means_long, aes( x = Category, y = RSD, group = Group ) ) +
  fte_theme_white() +
  geom_line( color = "cadetblue3", linewidth = 2 ) +
  geom_point( color = "cadetblue4", size = 4 ) +
  theme( axis.text.x = element_text( size = 14 ),
         axis.text.y = element_text( size = 14 ),
         axis.title.x = element_blank() ) +
  labs(  y = "RSD on Oxide (%)")
raw.rsd.plot


full.plot.order <-  c(colnames(oxide.order.plot), colnames(trace.order.plot))

wr.summary.table.export <- wr.summary.table[,full.plot.order]
# Reshape the data frame to long format
df_long <- tidyr::gather(wr.comp.model.comb, key = "oxide", value = "percent")

df_long <- df_long[complete.cases(df_long), ]
df_long$oxide <- factor( df_long$oxide, levels = full.plot.order)
df_long <- df_long[complete.cases(df_long), ]
# df_long <- df_long[order(df_long$oxide), ]

wr.comp.model.comb.summarize <- cbind( col.var = rep( 1, times = nrow( wr.comp.model.comb )),
  wr.comp.model.comb )

summary_stats <- wr.comp.model.comb.summarize %>%
  summarise(across(everything(), list(mean = mean, sd = sd))) %>%
  pivot_longer(cols = -1, names_to = c("oxide", ".value"), names_sep = "_")
summary_stats <- summary_stats[-1,]

summary_stats_plotting <- summary_stats %>%
  filter(oxide %in% full.plot.order) %>%
  arrange(factor(oxide, levels = full.plot.order))
summary_stats_plotting <- summary_stats_plotting[order(summary_stats_plotting$oxide), ]


summary_stats_plotting$oxide <- factor(summary_stats_plotting$oxide, levels = full.plot.order)
summary_stats_plotting <- summary_stats_plotting[order(summary_stats_plotting$oxide), ]

# Plot facet_wrap with geom_density
raw.sd.plot <- ggplot(df_long, aes(x = percent)) +
  fte_theme_white() +
  theme( axis.text.x = element_text( size = 10 ),
         axis.text.y = element_text( size = 10 ),
         axis.title.x = element_text( size = 20 ),
         axis.title.y = element_text( size = 20 ),
         strip.background = element_rect(
           color="black", fill="#213e47", size=0.2, linetype="solid"
         ),
         strip.text.x = element_text(
           size = 12, color = "white" ) ) +
  geom_histogram( color = "black", fill = "#645153", size = 0.1, bins = 30 ) +
  geom_vline(data = summary_stats_plotting, aes( xintercept = mean), color = "#3A8997", linetype = "dashed", size = 1) +
  geom_vline(data = summary_stats_plotting, aes(xintercept = mean + 2 * sd), color = "#8f7767", linetype = "dotted", size = 1) +
  geom_vline(data = summary_stats_plotting, aes(xintercept = mean - 2 * sd), color = "#8f7767", linetype = "dotted", size = 1) +
  labs( y = "", x = "Wt% Oxide/ppm") +
  theme(aspect.ratio = 1) +
  facet_wrap(~ oxide, ncol = 4, scales = "free", )
raw.sd.plot




wr.summary.table <- wr.summary.table[, c(full.plot.order)]
wr.summary.table <- data.frame( mass = c( rep( sample.mass, times = nrow(wr.summary.table))),
                                wr.summary.table )



