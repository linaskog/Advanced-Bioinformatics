# input_panel

input_panel <-
  tabPanel("Input Gene ID",
           card(
             card_header = "Gene Input",
             "Please enter valid Ensembl Gene ID below and press submit to",
             "find sgRNAs.",
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