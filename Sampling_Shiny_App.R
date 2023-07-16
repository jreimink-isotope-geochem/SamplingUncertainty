
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
  h1( id="big-heading", "Geochemical Uncertainty in Geological Sampling"),
  tags$style( HTML( "#big-heading{ color: #213e47;
                    text-align:center}")),
  sidebarPanel(
    tags$div(
      style = "background-color: #f2f2f2;",
      p(style = "text-align: left;", 
        HTML("This app was written by Jesse Reimink and Tom Chacko 
        and is largely based on the sampling error algorithms of 
        Stanley (2003) but also incorporates the effects of closure on 
        uncertainty estimates associated with whole-rock compositions
             ")),
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
    downloadButton("dwnldTableBtn", "Download Table"),
    br(),
    actionButton("runModelBtn", "Run Modeling"),),
  
  mainPanel(
    tabsetPanel(

      # Tab 1: Input table -       
      tabPanel( "Instructions", 
                tags$div(
                  tags$hr(style = "border: 1px solid black; margin-top: 10px;"),
                h5("Instructions:"),
                tags$hr(style = "border: 1px solid black; margin-top: 10px;"),
                p(HTML("<ol>
      <li>To enter the input data directly in the app, click on the <b>Input Data</b> tab above.
      Note that this is a one-time data entry where input will be lost when the app is closed. 
      Alternatively, select the <b>Download Table</b> button to the left and enter and 
      save input data in an Excel (csv format) table. Do NOT change the formatting of this table,
      save it as a .csv file, and re-upload it. To re-upload, select <b>Upload File</b> button 
      on right to upload your .csv input file. Mineral mode inputs are in volume % and 
      grain dimensions in centimeters. Distribution types: h = hypergeometric; 
      b = binomial; p = Poisson. Mineral compositions (from EPMA data) are in 
      weight % oxide.</li>
      <li>In the <b>Input Data</b> tab window, enter the sample name, which will be saved to 
      output tables and plots.</li>
      <li>Select the <b>Type of Calculation</b> (left side panel): <i>Range of Sample Sizes</i> 
      models the rock over a range of sample sizes (masses) from 625 to 20,000 grams</li>
      <li>Press the <b>Save Table</b> button (left panel) - this must be done for the program to
      run the calculation for the desired rock composition. Additionally, it is recommended that you 
      download the final table before running by using the <b>Download Table</b> button and saving the 
      file for your records. </li>
      <li>Press <b>Run Modeling</b> (left panel) and wait for program to finish.</li>
      <li>View and download output table and plots using the output tabs above. Output table (<i>Results</i>) gives the average major-, 
                       minor- and selected trace-element composition derived from the 
                       10,000 simulations of the input rock modal composition. The table 
                       also reports, for a given sample mass and mineral grain sizes, 
                       the predicted sampling error standard deviation and relative 
                       standard deviation (RSD) for each oxide or element.</li>" ) ) ) ),
      # Tab 2: Input table - 
      tabPanel("Input Data", 
               textInput( "sampleName", "Sample Name:", value = "sample name here" ),
               br(),
               selectizeInput(
                 inputId = "inputOxides",
                 label = "Select and Order Oxides",
                 choices = c("SiO2", "TiO2", "Al2O3", "Fe2O3", "FeO", "FeOt", "MnO", "MgO", "CaO",
                             "Na2O", "K2O", "P2O5", "H2O", "CO2", "Li2O", "ThO2", "BaO", "ZrO2", 
                             "La2O3", "Ce2O3", "Pr2O3", "Nd2O3", "Sm2O3", "Eu2O3", "Gd2O3", "Tb2O3",
                             "Dy2O3", "Ho2O3", "Er2O3", "Tm2O3", "Yb2O3", "Lu2O3", "Y2O3" ),
                 selected = c("SiO2", "TiO2", "Al2O3", "FeOt", "MnO", "MgO", "CaO",
                              "Na2O", "K2O", "P2O5" ),
                 multiple = TRUE
               ),
               selectizeInput(
                 inputId = "inputTraces",
                 label = "Select and Order Trace Elements",
                 choices = c("Zr", "Th", "Ba", "Li", "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", 
                             "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Y" ),
                 selected = c( "Zr", "Th", "Ba" ),
                 multiple = TRUE
               ),
               br(),
               verbatimTextOutput("sumOutput"),
               br(),
               
               radioButtons(
                 inputId ="input_type",
                 label = 'Data Input Type',
                 choiceNames = c("Edit File", "Upload File"),
                 choiceValues = c('Edit', 'Upload'),
                 selected = "Edit",
                 inline = TRUE ),
               #deactivate input sample size field when range of samples is selected
               conditionalPanel(
                 condition = "input.input_type == 'Upload'",
                 fileInput("fileInput", "Upload CSV File"),
                 dataTableOutput( "tableOutput")
               ),
               p("Dimensions in cm, densities in kg/cm3"),
               rHandsontableOutput("MineralParameters") ),
      
      # Tab 2: Output table - 
      tabPanel("Results Table", 
               h4("Whole-rock Means"), dataTableOutput("output_means"),
               downloadButton("dwnldMeanTablebtn", "Download Means Table" ) ),
               # br(),
               # h4("Whole-rock Relative Standard Deviations"), dataTableOutput("output_rsds"),
               # downloadButton("dwnldRSDbtn", "Download RSD Table" ),
               # br(),
               # h4("Whole-rock Standard Deviations"), dataTableOutput("output_sds"),
               # downloadButton("dwnldSDTablebtn", "Download StDev Table" ) ),
      
      
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
server <- function(input, output, session) {
  
  # Read the data from CSV
  min_parameters <- reactiveVal(read.csv("min_parameters_input.csv", stringsAsFactors = FALSE))

  ## created the selected columns for export
  selectedOxides <- reactive({
    input$inputOxides
  })
  selectedTraces <- reactive({
    input$inputTraces
  })
  
  
  
  
  uploadedData <- reactive({
    req(input$fileInput)
    inFile <- input$fileInput
    if (!is.null(inFile)) {
      # Read the uploaded file
      input.data <- read.csv(inFile$datapath, stringsAsFactors = FALSE)
      return(input.data)
    }
  })
  

  
  
  # output$MineralParameters <- renderRHandsontable({
  #   if (input$dataOption = "Edit") {
  #     rhandsontable(min_parameters(), rowHeaders = NULL, width = "100%") %>%
  #       hot_col(col = "Distribution", type = "dropdown", source = c("h", "b", "p")) %>%
  #       hot_cols(list(editable = TRUE)) %>%
  #       hot_col("Modes", renderer = "
  #       function(instance, td, row, col, prop, value, cellProperties) {
  #         Handsontable.renderers.TextRenderer.apply(this, arguments);
  #         if (parseFloat(value) !== 0) {
  #           $(td).css('background-color', 'lightyellow');
  #         }
  #       }
  #     ")
  #   } else {
  #     if (!is.null(input$fileInput ) ) {
  #       inFile <- input$fileInput
  #       data <- read.csv(inFile$datapath, stringsAsFactors = FALSE)
  #       rhandsontable(data, rowHeaders = NULL, width = "100%") %>%
  #         hot_col(col = "Distribution", type = "dropdown", source = c("h", "b", "p")) %>%
  #         hot_cols(list(editable = TRUE)) %>%
  #         hot_col("Modes", renderer = "
  #         function(instance, td, row, col, prop, value, cellProperties) {
  #           Handsontable.renderers.TextRenderer.apply(this, arguments);
  #           if (parseFloat(value) !== 0) {
  #             $(td).css('background-color', 'lightyellow');
  #           }
  #         }
  #       ")
  #     }
  #   }
  # })
  observe( 
    if( input$input_type == "Upload" ) {
      output$MineralParameters <- renderRHandsontable({
        rhandsontable(uploadedData(), rowHeaders = NULL, width = "100%" ) %>%
          hot_cols( columnSorting = T, manualColumnResize = T ) %>%
          hot_col(col = "Distribution", type = "dropdown", source = c("h", "b", "p")) %>%
          hot_cols(list(editable = TRUE), fixedColumnsLeft = 1, columnSorting = T, manualColumnResize = T ) %>%
          hot_rows(fixedRowsTop = 1) %>%
          hot_col("Modes", renderer = "
        function(instance, td, row, col, prop, value, cellProperties) {
          Handsontable.renderers.TextRenderer.apply(this, arguments);
          if (parseFloat(value) !== 0) {
            $(td).css('background-color', 'lightyellow');
          }
        }
      ")
      })
    }
    else {
      output$MineralParameters <- renderRHandsontable({
        rhandsontable(min_parameters(), rowHeaders = NULL ) %>%
          # hot_cols( colwidths = c( "100px", "300px", "100px") ) %>%
          # hot_cols( ) %>%
          hot_col(col = "Distribution", type = "dropdown", source = c("h", "b", "p")) %>%
          hot_cols(list(editable = TRUE), fixedColumnsLeft = 1, columnSorting = T, manualColumnResize = T ) %>%
          hot_col("Modes", renderer = "
        function(instance, td, row, col, prop, value, cellProperties) {
          Handsontable.renderers.TextRenderer.apply(this, arguments);
          if (parseFloat(value) !== 0) {
            $(td).css('background-color', 'lightyellow');
          }
        }
      ")
      })
    }
    )


  observeEvent(input$saveBtn, {
    if( input$input_type == "Upload" ) {
      newData <- hot_to_r(input$MineralParameters)
    }
    else {
      newData <- hot_to_r(input$MineralParameters)
    }
    sample.size.data <- data.frame( sample.size = input$sampleSize )
    # Replace "path/to/file.csv" with the actual path to your CSV file
    write.csv(newData, "Shiny_OutPut.csv", row.names = FALSE)
    write.csv(newData, paste(input$sampleName,"_min_parameters.csv", sep=""), row.names = FALSE)
    write.csv(sample.size.data, "sample_size.csv", row.names = FALSE)
    write.csv(selectedOxides(), "oxideschosen.csv", row.names = FALSE)
    write.csv(selectedTraces(), "traceschosen.csv", row.names = FALSE)
    showModal(modalDialog(
      "Changes saved successfully!",
      title = "Success"
    ))
  })
  
  # Update the data when the table is edited
  observeEvent(input$MineralParameters$changes$changes, {
    min_parameters(hot_to_r(input$MineralParameters))
  })
  
  # Auto-sort the data based on the 'Modes' column
  observeEvent(
    min_parameters(), {
    sorted_data <- min_parameters() %>%
      arrange( desc( Modes ) )
    min_parameters(sorted_data)
  })
  
  # Render the Save Table button
  output$saveButton <- renderUI({
    color <- if (values$tableChanged) "green" else "red"
    actionButton("saveBtn", "Save Table", style = sprintf("color: white; background-color: %s;", color))
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
  
  ## download input data tables when button pushed
  output$dwnldTableBtn <- downloadHandler(
    newData <- hot_to_r(input$MineralParameters),
    filename = function() {
      paste(input$sampleName,"Input_Data_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv( hot_to_r(input$MineralParameters), file, row.names = FALSE )
    }
  )
  
  
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
      summary.means <- env$wr.summary.table.export
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
        datatable(summary.means, options = list(pageLength = 30, info = FALSE))
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
    
    # output$dwnldRSDbtn <- downloadHandler(
    #   filename = function() {
    #     paste(input$sampleName,"Model_RSDs_", Sys.Date(), ".csv", sep = "")
    #   },
    #   content = function(file) {
    #     if (input$calc_type == "Single") {
    #       NULL
    #     } else if ( input$calc_type == "Range") {
    #       write.csv( rsds, file, row.names = TRUE )
    #     }
    #   }
    # )
    # 
    # output$dwnldSDTablebtn <- downloadHandler(
    #   filename = function() {
    #     paste(input$sampleName,"Model_SDs_", Sys.Date(), ".csv", sep = "")
    #   },
    #   content = function(file) {
    #     if (input$calc_type == "Single") {
    #       NULL
    #     } else if ( input$calc_type == "Range") {
    #       write.csv( sds, file, row.names = TRUE )
    #     }
    #   }
    # )
    

    
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

