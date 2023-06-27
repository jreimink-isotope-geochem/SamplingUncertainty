
options(shiny.sanitize.errors = FALSE)
library(shiny)
library(rhandsontable)
library(ggplot2)
library(DT)
library(dplyr)
library(viridis)
library(tidyr)                            

# Define UI
ui <- fluidPage(
  
  titlePanel(""),
  sidebarPanel(
    tags$div(
      style = "background-color: #f2f2f2;",
      h4(style = "text-align: center;", "Geochemical Uncertainty in Geological Sampling"),
      p(style = "text-align: center;", 
        HTML("This app was written by <br> 
             Drs. Jesse Reimink and Tom Chacko <br> 
             and is largely based on the sampling error algorithms of <br> 
             Stanley (2003) but also incorporate the effects of closure on <br>
             uncertainty estimates associated with whole-rock compositions")),
      tags$hr(style = "border: 1px solid black; margin-top: 10px;"),
      h5("Instructions:"),
      p(HTML("<ol>
      <li>Enter data in the <b>Input Data</b> tab below</li>
      <li>Enter the sample name that will be saved to output tables and plots</li>
      <li>Select the <b>Type of Calculation</b> below. <i>Range of Sample Sizes</i> models 
      the rock over a range of sample sizes (masses), 
      while the <i>Defined Sample Size</i> uses a user-defined sample size</li>
      <li>Press the <b>Save Table</b> button - this must be done for the program to run
      the calculation for the desired rock composition</li>
      <li>Press <b>Run Modeling</b> and wait for program to finish</li>
      <li>View and download output files in the tabs in the output section</li>")),
    ),
    tags$hr(style = "border: 1px solid black; margin-top: 10px;"),      
    radioButtons(
      inputId ="calc_type",
      label = 'Type of Calculation',
      choiceNames = c("Range of Sample Sizes", "Defined Sample Size"),
      choiceValues = c('Range', 'Single'),
      selected = "Range",
      inline = TRUE ),
    #deactivate input sample size field when range of samples is selected
    conditionalPanel(
      condition = "input.calc_type == 'Single'",
      numericInput( "sampleSize", "Sample Size (in grams):", value = 625, min = 1 )
    ),
    br(),
    # uiOutput("saveButton"),
    actionButton("saveBtn", "Save Table"),
    br(),
    actionButton("runModelBtn", "Run Modeling"),),
  
  mainPanel(
    tabsetPanel(
      # Tab 1: Output table
      
      
      tabPanel("Input Data", 
               verbatimTextOutput("sumOutput"),
               br(),
               textInput( "sampleName", "Sample Name:", value = "sample name here" ),
               br(),
               p("Dimensions in cm, densities in kg/cm3"),
               rHandsontableOutput("MineralParameters") ),
      
      # Tab 2: Output table - 
      tabPanel("Results Table", 
               h4("Whole-rock Means"), dataTableOutput("output_means"),
               downloadButton("dwnldMeanTablebtn", "Download Means Table" ),
               br(),
               h4("Whole-rock Relative Standard Deviations"), dataTableOutput("output_rsds"),
               downloadButton("dwnldRSDbtn", "Download RSD Table" ),
               br(),
               h4("Whole-rock Standard Deviations"), dataTableOutput("output_sds"),
               downloadButton("dwnldSDTablebtn", "Download StDev Table" ) ),
      
      
      # tabPanel("Oxide Plot", plotOutput("OxidePlotChangeable"),
      #          br(),
      #          selectInput(
      #            inputId = "column",
      #            label = "Select an oxide:",
      #            choices = c( "SiO2", "TiO2", "Al2O3", "FeO", "MnO", "MgO", "CaO", "Na2O", "K2O")
      #          ) ),
      
      #Tab 3: Output plot - Discordance Plot (lower and Upper int)
      tabPanel("Modeled Oxide RSDs", plotOutput("RSDsPlot"),
               br(),
               downloadButton( "dwnldRSDplotbtn", "Download Plot" ) ),
      
      # Tab 4: Output plot - Download StDev plot
      tabPanel("Standard Deviations", plotOutput("SDsPlot"),
               br(),
               downloadButton("dwnlSDplotbtn", "Download Plot" ) ),
    )
  ),
  # CSS styling
  tags$head(
    tags$style(HTML("
      #runModelBtn {
        background-color: #213e47;
        color: white;
      }
      
      #runModelBtn:hover {
        background-color: #3A8997;
      }
      
            #saveBtn {
        background-color: #f3d098;
        color: #000000;
      }
      
      #saveBtn:hover {
        background-color: #645153;
        color: #ffffff
      }
      
    ") ) )
)

# Define server logic
server <- function(input, output) {
  
  # Read the data from CSV
  min_parameters <- reactiveVal(read.csv("min_parameters.csv", stringsAsFactors = FALSE))
  ## oxide to plot
  data.oxide <- reactiveVal(NULL)
  
  # Render the min_density table
  output$MineralParameters <- renderRHandsontable({
    rhandsontable(min_parameters(), rowHeaders = NULL, width = "100%" ) %>%
      hot_col(col = "Distribution", type = "dropdown", source = c("h", "b", "p")) %>%
      hot_cols(list(editable = TRUE)) %>%
      hot_col("Modes", renderer = "
        function(instance, td, row, col, prop, value, cellProperties) {
          Handsontable.renderers.TextRenderer.apply(this, arguments);
          if (parseFloat(value) !== 0) {
            $(td).css('background-color', 'lightyellow');
          }
        }
      ")
  })
  

  
  # Update the data when the table is edited
  observeEvent(input$MineralParameters$changes$changes, {
    min_parameters(hot_to_r(input$MineralParameters))
  })
  
  # Auto-sort the data based on the 'Name' column
  observeEvent(min_parameters(), {
    sorted_data <- min_parameters() %>%
      arrange( desc(Modes))
    
    min_parameters(sorted_data)
  })
  
  # Render the Save Table button
  output$saveButton <- renderUI({
    color <- if (values$tableChanged) "green" else "red"
    actionButton("saveBtn", "Save Table", style = sprintf("color: white; background-color: %s;", color))
  })

  observeEvent(input$saveBtn, {
    newData <- hot_to_r(input$MineralParameters)
    sample.size.data <- data.frame( sample.size = input$sampleSize )
    # Replace "path/to/file.csv" with the actual path to your CSV file
    write.csv(newData, "Shiny_OutPut.csv", row.names = FALSE)
    write.csv(newData, paste(input$sampleName,"_min_parameters.csv", sep=""), row.names = FALSE)
    write.csv(sample.size.data, "sample_size.csv", row.names = FALSE)
    showModal(modalDialog(
      "Changes saved successfully!",
      title = "Success"
    ))
  })
  
  ## Sum up the mode data and display
  sumModes <- reactive({
    input$MineralParameters %>%
      hot_to_r() %>%
      dplyr::select(Modes) %>%
      summarise("Total Mode" = sum(Modes))
  })
  
  # Display the sum of the Modes column
  output$sumOutput <- renderPrint({
    sumModes()
  })
  
  # do stuff with the sample size value
  observeEvent(input$sampleSize, {
    sampleSize <- input$sampleSize
    # Do something with the sampleSize value
    # ...
  })
  

  
  # Run the reduction with an if statement depending on whether it is pertaining to many sample sizes or
  #. just a few
  observeEvent(input$runModelBtn, {
    showModal(
      modalDialog(
        "The Modeling Process is Running",
        easyClose = TRUE
      )
    )

    # Set data1 as an environment variable for the sourced script
    env <- new.env()
    
    # Source the R script in the created environment
      if (input$calc_type == "Single") {
        source("Sampling_Uncertainty_Modeling_ForShinyOneSampleSize.R", local = env)
        
      } else if (input$calc_type == "Range") {
        # Load file for "Range of Samples" option
        # Source the different file or perform necessary operations
        source("Sampling_Uncertainty_Modeling_ForShinyRangeOfSamples.R", local = env)
    }
    
    
    # Retrieve the result from the environment
    if (input$calc_type == "Single") {

      wr.summary.table <- env$wr.summary.table
      wr.comp.model.comb <- env$wr.comp.model.comb
      
    } else if (input$calc_type == "Range") {
      # Load file for "Range of Samples" option
      # Source the different file or perform necessary operations
      rsds <- env$summary.rsds
      sds <- env$summary.raw.sds
      summary.means <- env$summary.means
      wr.comp.model.comb <- env$wr.comp.model.comb
    }
    
    
    output$RSDsPlot <- renderPlot({
      env$raw.rsd.plot 
    })   
    output$SDsPlot <- renderPlot({
      env$raw.sd.plot
    })
    

    
    # Render the data.table
    output$output_means <- renderDataTable({
      if (input$calc_type == "Single") {
        datatable(wr.summary.table)
      } else if (input$calc_type == "Range") {
        datatable(summary.means)
      }
    })
    
    ## Render other tables if Range is selected
    output$output_rsds <- renderDataTable({
      if (input$calc_type == "Single") {
        NULL
      } else if (input$calc_type == "Range") {
        datatable(rsds)
      }
    })

    output$output_sds <- renderDataTable({
      if (input$calc_type == "Single") {
        NULL
      } else if (input$calc_type == "Range") {
        datatable(sds)
      }
    })
    
    observeEvent( input$runModelBtn, {
      showModal(
        modalDialog(
          if (input$calc_type == "Single") {
            print("Single Sample Size Calculation Completed")
          } else if (input$calc_type == "Range") {
            print("Range of Sample Sizes Calculation Completed")
          }
          
        )
      )
    })

    
    # Render the data.table
    output$myTable <- renderDataTable({
      
      if (input$calc_type == "Single") {
        datatable( wr.summary.table )
        
      } else if (input$calc_type == "Range") {
        # Load file for "Range of Samples" option
        # Source the different file or perform necessary operations
        datatable( rsds )
      }
    })

    ## download data tables
    output$dwnldMeanTablebtn <- downloadHandler(
      filename = function() {
        paste(input$sampleName,"Model_Means_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        if (input$calc_type == "Single") {
          write.csv( wr.summary.table, file, row.names = TRUE )
        } else if (input$calc_type == "Range") {
          write.csv( summary.means, file, row.names = TRUE )
        }
      }
    )
    
    output$dwnldRSDbtn <- downloadHandler(
      filename = function() {
        paste(input$sampleName,"Model_RSDs_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        if (input$calc_type == "Single") {
          NULL
        } else if ( input$calc_type == "Range") {
          write.csv( rsds, file, row.names = TRUE )
        }
      }
    )
    
    output$dwnldSDTablebtn <- downloadHandler(
      filename = function() {
        paste(input$sampleName,"Model_SDs_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        if (input$calc_type == "Single") {
          NULL
        } else if ( input$calc_type == "Range") {
          write.csv( sds, file, row.names = TRUE )
        }
      }
    )
    
    ## download the mass vs RSD for each oxide plot
    output$dwnlSDplotbtn <- downloadHandler(
      filename = function() {
        paste(input$sampleName,"Model_SDs_", Sys.Date(), ".pdf", sep = "") # Set the filename for the downloaded file
      },
      content = function(file) {
        # Retrieve the saved plot and save it to the file
        plot <- env$raw.sd.plot
        # Set the dimensions of the exported PDF
        width <- 14  # Specify the desired width in inches
        height <- 14  # Specify the desired height in inches
        
        # Save the plot to the file as a PDF with the specified dimensions
        pdf(file, width = width, height = height)
        print(plot)
        dev.off()
      }
    )
    
    
    output$dwnldRSDplotbtn <- downloadHandler(
      filename = function() {
        paste(input$sampleName,"Model_RSDs_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        # Retrieve the saved plot and save it to the file
        plot <- env$raw.rsd.plot
        # Set the dimensions of the exported PDF
        width <- 14  # Specify the desired width in inches
        height <- 8  # Specify the desired height in inches
        
        # Save the plot to the file as a PDF with the specified dimensions
        pdf(file, width = width, height = height)
        print(plot)
        dev.off()
      }
    )
  }
  )
  
  observeEvent( input$column, {
    output$OxidePlotChangeable <- renderDataTable({
      ggplot( wr.comp.model.comb, aes( x = wr.comp.model.comb$input$column ) ) +
        geom_histogram()
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server)

