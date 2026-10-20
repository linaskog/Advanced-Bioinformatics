# app_score_panel.R

score_panel <- tabPanel(
  "Score",
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
                     "Rank Candidate Guides By:",
                     choices = list(
                       "Linear Model" = "LinearScore", 
                       "XGBoost Model" = "XGBScore",
                       "Combination of Models" = "Combination"),
                     selected = "Combination"),
        sliderInput("linear_weight",
                    "Proportion Linear Model:",
                    min = 0, max = 1, value = 0.5, step = 0.05),
        br(),
        numericInput("top_n", "Number of sgRNAs:",
                     value = 50, min = 1, max = 2000),
        sliderInput("max_difference",
                    "Maximal Difference Accepted (%):",
                    min = 0, max = 100, value = 25, post = "%"),
        helpText("The maximum percentage difference between the ",
                 "linear and XGBoost model predictions can be set ", 
                 "to flag guides with large variability in model ",
                 "scores.")
      ),
      
      mainPanel(
        h4("XGBoost versus Linear Model Predicted Scores"),
        plotly::plotlyOutput("plot"),
        br(),
        h4("Top Candidate sgRNAs"),
        DT::DTOutput("top_sgrna_table")
      )
    )
    )
  ) 
             
