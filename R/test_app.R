library(shiny)
library(bslib)

# Define UI for app that draws a histogram ----
ui <- page_fluid(
  # App title ----
  title = "Welcome to the ...thing",
  # Sidebar panel for inputs ----
  card(plotOutput(outputId = "distPlot"),
       sliderInput(
         inputId = "bins",
         label = "Number of bins:",
         min = 5,
         max = 50,
         value = 30
       )
       ),
    
  card(plotOutput(outputId = "distPlot2"),
       sliderInput(
         inputId = "bins2",
         label = "Number of bins:",
         min = 0,
         max = 25,
         value = 10))
  # Output: Histogram ----
  
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  output$distPlot <- renderPlot({
    
    x    <- faithful$waiting
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    hist(x, breaks = bins, col = "#007bc2", border = "orange",
         xlab = "Waiting time to next eruption (in mins)",
         main = "Histogram of waiting times")
    
  })
  
  output$distPlot2 <- renderPlot({
    x2 <- faithful$eruptions
    bins2 <- seq(min(x2), max(x2), length.out = input$bins2 + 1)
    
    hist(x2, breaks = bins2, col = "green", border = "orange",
         xlab = "Eruptions",
         main = "Histogram of eruptions")
    
  })
  
}

shinyApp(ui = ui, server = server)
