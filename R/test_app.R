# Reading in Necessary Libraries ------------------------------------------

library(shiny)
library(bslib)


# Creating a User Interface -----------------------------------------------
# QC-control
qc_min = 1
qc_max = 50
qc_value = 30

score_min = 1
score_max = 50
score_value = 30

ui <- navbarPage(
  "sgRNA Analyzer",
  
  # Input panel
  tabPanel("Input Gene ID",
           card(
             card_header = "Gene Input",
             "Please enter valid Ensembl Gene ID below and press submit to",
             "find sgRNAs.",
             textInput(
               "geneID",
               label = "Ensembl Gene ID",
               value = "ENSG00000107485"
             ),
             actionButton("submit", 
                          label = "Submit",
                          style = 'width:150px'
             ),
           )
  ),
  
  # Quality Control Plot panel
  tabPanel("Quality Control",
           card(
             card_header = "Quality Control",
             plotOutput(outputId = "qcplot"),
             sliderInput(
               inputId = "binsqc",
               label = "Number of bins:",
               min = qc_min,
               max = qc_max,
               value = qc_value
             )
           )
  ),
  
  # sgRNA score plot
  tabPanel("sgRNA Score Plot",
           card(
             card_header = "sgRNA Efficacy",
             plotOutput( outputId = "scoreplot"),
             sliderInput(
               inputId = "binsscore",
               label = "Number of bins:",
               min = score_min,
               max = score_max,
               value = score_value
             )
           )
  )
) 


server <- function(input, output) {
  
  output$qcplot <- renderPlot({
    x    <- faithful$waiting
    bins <- seq(min(x), max(x), length.out = input$binsqc + 1)
    
    hist(x, breaks = bins, col = "#007bc2", border = "orange",
         xlab = "Waiting time to next eruption (in mins)",
         main = "Histogram of waiting times")
    
  })
  
  output$scoreplot <- renderPlot({
    x2 <- faithful$eruptions
    bins2 <- seq(min(x2), max(x2), length.out = input$binsscore + 1)
    
    hist(x2, breaks = bins2, col = "green", border = "orange",
         xlab = "Eruptions",
         main = "Histogram of eruptions")
    
  })
  
}

shinyApp(ui = ui, server = server)
