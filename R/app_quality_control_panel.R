# app_quality_control_panel.R
# quality_control_panel object for Rshiny interface ("R/test_app.R")

quality_control_panel <- 
  tabPanel("Quality Control",
         card(
           card_header = "Quality Control",
           textOutput(outputId = "qc_guide_count"),
           textOutput(outputId = "qc_mean_gc"),
           
           sliderInput(
             inputId = "gc_range",
             label = "GC-content reference range (%)",
             min = 0,
             max = 100,
             value = c(40,70),
             step = 5,
             post = "%"
           ),
           textOutput(outputId = "qc_gc_range_summary"),
           
           plotOutput(outputId = "qcplot"),
           helpText(
             "The shaded region represents the user-selected GC-content range. ",
             "GC content is one sequence-quality indicator, but it does not ",
             "directly measure sgRNA editing efficiency."
           ),
           tags$hr(),
           h4("Candidate sgRNA Quality Table"),
           DT::DTOutput(outputId = "guide_quality_table"),
           helpText(
             "Sequence QC flag: 'No flags' means the guide passes the selected ",
             "GC-content range and contains no detected poly-T (TTTT) or ",
             "other homopolymer runs of four or more nucleotides. ",
             "'Review' indicates at least one of these criteria was flagged. ",
             "These checks do not measure editing efficiency or off-target specificity."
           )
           
         )
)