# app_input_panel.R
# Encodes input_panel object for rshiny interface ("R/test_app.R")

input_panel <-
  tabPanel("Input Gene ID",
           card(
             card_header = "Gene Input",
             "Please enter valid Ensembl Gene ID below and press submit to",
             "find 20 nucleotide long sgRNAs.",
             br(),
             br(),
             textInput(
               "geneID",
               label = "Ensembl Gene ID",
               value = "ENSG00000107485"
             ),
             br(),
             actionButton(
               "submit",
               label = "Submit",
               style = 'width:150px'
             ),
             br(),
             textOutput(
               "output_text"
               )
             )
           )