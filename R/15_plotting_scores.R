# 15_plotting_scores.R

# Purpose:
# Explore score plots and tables for interactive rshiny interface

# Notes:
# Uses test dataset predictions instead of generated guides for input genes

# Input:
# - Linear model
# - XGB model
# - variables and testing dataset index

# Reading in necessary libraries ------------------------------------------

library(shiny)
library(xgboost)
library(ggplot2)
library(plotly)
library(DT)

# Reading in Models -------------------------------------------------------

xgb_model <- xgb.load(
  "~/GitHub/Advanced-Bioinformatics/models/xgb_model.ubj")
linear_model <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/models/linear_model.rds")


# Reading in data ---------------------------------------------------------

variables <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/variables.rds")
test_index <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/data/processed/test_index.rds")



# Preparing data ----------------------------------------------------------

pos <- as.matrix(variables[, grep("pos", colnames(variables))])
additional_data <- as.matrix(variables[, c("GC", "strand", "chromosome", 
                                           "region_type", "gene_overlap_count", 
                                           "H3K27ac_peak_overlap")
                                       ]
                             )
X <- cbind(pos, additional_data)


df <- data.frame(
  sgRNA_id = variables[test_index, "UID"],
  guide = variables[test_index, "Guide"],
  linear_predictions = predict(linear_model, variables[test_index,]),
  xgb_predictions = predict(xgb_model, X[test_index,],)
)


# Interactive Plot --------------------------------------------------------

# User Interface
ui <- navbarPage(title = "sgRNA analysis",
                 
  # Welcome Panel
  tabPanel("Welcome",
           h2("Welcome Panel"),
           helpText("This is the helptext for the welcome panel.")
           ),
  # Score Panel
  tabPanel("Score",
               fluidPage( 
                 sidebarLayout(
                   sidebarPanel(
                     h4("Settings Panel"),
                     helpText(
                       "Choose a model for ranking sgRNAs (linear, XGBoost or ",
                       "combination). For combinated model, please select ", 
                       "weight of the linear model (0-1). The number of ",
                       "selected sgRNAs ranked according to the selected model",
                       " are highlighted in the plot and listed in the table ",
                       "below."
                       ),
                     selectInput( "select",
                                  "Rank candidate guides by:",
                                  choices = list(
                                    "Linear Model" = "linear_predictions", 
                                    "XGBoost Model" = "xgb_predictions",
                                    "Combination of Models" = "combination"),
                                  selected = "combination"),
                     sliderInput("linear_weight",
                                 "Proportion Linear Model:",
                                 min = 0, max = 1, value = 0.5, step = 0.05),
                     br(),
                     numericInput("top_n", "Number of sgRNAs:",
                                  value = 50, min = 1, max = 2000),
                     sliderInput("max_difference",
                                 "Maximal difference accepted:",
                                 min = 0, max = 100, value = 25, post = "%"),
                     helpText("The maximum percentage difference between the ",
                              "linear and XGBoost model predictions can be set ", 
                              "to flag guides with large variability in model ",
                              "scores.")
                     ),
                   
                   mainPanel(
                     h4("XGBoost vs. Linear Model Predicted Scores"),
                     plotly::plotlyOutput("plot"),
                     br(),
                     h4("Top Candidate sgRNAs"),
                     DT::DTOutput("top_sgrna_table")
                     )
                   )
                 )
           )
  )

# Defining Server Logic
server <- function (input, output) {
  
  # Creating Dataframe for Score Evaluation
  evaluated_df <- reactive ({
    d <- df
    weight_linear <- input$linear_weight
    weight_xgb <- 1 - weight_linear
    
    # Computing combined model score
    d$combination <- (weight_linear*d$linear_predictions) + 
      (weight_xgb*d$xgb_predictions)
    
    # Computing difference between model scores
    d$difference <- abs(100 * (d$linear_predictions - d$xgb_predictions)/
      (d$linear_predictions + d$xgb_predictions)/2)
    
    # Ordering sgRNAs for selected model
    d[order(d[[input$select]]), ]
    }
  )
  
  # Creating Interactive Plot of sgRNA scores
  output$plot <- plotly::renderPlotly({
    d <- evaluated_df()
    top_d <- head(d, min(input$top_n, nrow(d)))
    
    p <- ggplot(
      data = d, 
      aes(
        x = linear_predictions,
        y = xgb_predictions,
        colour = combination,
        text = paste0("<b>sgRNA: </b>", sgRNA_id, "<br>",
                      "<b>Guide: </b>", guide, "<br>", 
                      "<b>Combination: </b>", round(combination, 4), "<br>", 
                      "<b>Linear: </b>", round(linear_predictions, 4), "<br>",
                      "<b>XGBoost: </b>", round(xgb_predictions, 4), "<br>",
                      "<b>Difference: </b>", round(difference, 2), "<br>")
                )) + 
      geom_point(alpha = 0.3) +
      geom_point(
        data = top_d,
        aes(x = linear_predictions,
            y = xgb_predictions),
            colour = "maroon") +
      labs( x = "Linear Model Prediction",
            y = "XGBoost Model Prediction")+
      theme_light()

    plotly::ggplotly(p, tooltip = "text")
    
  })
  
  # Creating Data Table for Candidate sgRNAs
  output$top_sgrna_table <- DT::renderDT({
    d <- evaluated_df()
    d <- head(d, min(input$top_n, nrow(d)))
    
    score_table <- 
      data.frame(
        "ID" = d$sgRNA_id, 
        "Guide" = d$guide,
        "Linear Score" = round(d$linear_predictions, 4),
        "XGBoost Score" = round(d$xgb_predictions, 4),
        "Combined Score" = round(d$combination, 4),
        "Difference %" = round(d$difference, 2),
        check.names = FALSE)
    
    # Creating Column for Flagging guides with high score variability
    score_table[["Flagged"]] <- ifelse(
      score_table[["Difference %"]] > input$max_difference,
      "YES",
      "NO"
    )
  
  # Formatting datatable
  dt <- DT::datatable(
    score_table,
    rownames = FALSE,
    options = list(
      pageLength=5,
      lengthMenu= c(5,10,25,50),
      scrollY = "300px",
      scrollX=TRUE
    )
  )
  
  # Adding color to Flagged column
  DT::formatStyle(
    dt,
    "Flagged",
    backgroundColor = DT::styleEqual(
      c("NO", "YES"),
      c("#EAF5EC", "#FFF3E0")
    ),
    color = DT::styleEqual(
      c("NO", "YES"),
      c("#216E39","#9A5A08")
    )
  )
  })
  
}

# Running App
shinyApp(ui = ui, server = server)
