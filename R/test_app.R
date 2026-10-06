# Reading in Datasets -----------------------------------------------------

# Human reference genome hg19
library(BSgenome.Hsapiens.UCSC.hg19)

# Human genome annotated with Ensembl IDs
library(EnsDb.Hsapiens.v75)

# Reading in Necessary Libraries ------------------------------------------

library(shiny)
library(bslib)
library(GenomicRanges)
library(GenomeInfoDb)
library(Biostrings)
library(xgboost)
library(stringr)


ensembl <- EnsDb.Hsapiens.v75
genome <- BSgenome.Hsapiens.UCSC.hg19


# Reading in Models -------------------------------------------------------

xgb_model <- xgb.load("xgb_model.ubj")

# Finding sgRNAs ----------------------------------------------------------

finding_sgRNAs <- function(gene_sequence, PAM = "NGG") {
  # Locate PAM sequences in gene of interest
  pam_sequences <- matchPattern(
    PAM,
    gene_sequence,
    fixed = FALSE
  )
  
  # Locating PAM start positions for valid (20 nt) sgRNAs
  pam_start <- start(pam_sequences)
  pam_start <- pam_start[pam_start > 20]
  
  # Creating IRange object for each sgRNA
  guide_ranges <- IRanges(
    start = pam_start - 20,
    end = pam_start -1
  )
  
  # Creating data frame for Guide sequences extracted using IRanges
  data.frame( "Guide" = as.character(extractAt(
    gene_sequence,
    guide_ranges))
  )
}

# Creating a User Interface -----------------------------------------------
# QC-control
qc_min = 1
qc_max = 20
qc_value = 10

score_min = 1
score_max = 30
score_value = 15

ui <- navbarPage(
  "sgRNA Analyzer",
  
  # Input panel
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
  
  submit_reactive <- eventReactive ( input$submit, {
    
    gene <- genes(
      ensembl,
      filter = ~gene_id == input$geneID
    )
    # Changing style to match datasets
    seqlevelsStyle(gene) <- "UCSC"
    # Extracting GOI from reference genome --> DNAString Class
    gene_sequence <- getSeq(genome, gene)[[1]]
    # Creating reverese complement of GOI
    gene_sequence_reverse <- reverseComplement(gene_sequence)
    
    # Creating data frame with all forward and reverse sgRNAs
    guides <- rbind(finding_sgRNAs(gene_sequence), 
                    finding_sgRNAs(gene_sequence_reverse))
    
    # Adding identification column to guide data frame
    guides$ID <- 1:nrow(guides)
    
    guides 
    })
    
  output$output_text <- renderText({
    guides <- submit_reactive()
    paste(input$geneID, "loaded,", nrow(guides), "sgRNAs found")
    })
  
  output$qcplot <- renderPlot({
    guides <- submit_reactive()
    guides$GC <- (str_count(guides$Guide, "G") + 
                    str_count(guides$Guide, "C")
    )/str_length(guides$Guide)
    
    bins <- seq(min(guides$GC), max(guides$GC), length.out = input$binsqc + 1)
    
    hist(guides$GC, breaks = bins, col = "lightgreen", border = "orange",
         xlab = "GC-ratio",
         main = "GC-ratio Distribution")
    
  })
  
  output$scoreplot <- renderPlot({
    guides <- submit_reactive()
    
    guides$GC <- (str_count(guides$Guide, "G") + 
                   str_count(guides$Guide, "C")
    )/str_length(guides$Guide)
    
    sequence_matrix <- str_split_fixed(
      guides$Guide, 
      pattern = "", 
      n = 20)
    
    X <- do.call(
      cbind,
      lapply(1:20, function(i) {
        model.matrix(~ sequence_matrix[, i] - 1)
      })
    )
    
    colnames(X) <- paste0(
      rep(paste0("pos", 1:20), each = 4),
      "_",
      rep(c("A", "C", "G", "T"), 20)
    )
    
    X <- cbind(X, GC = guides$GC)
    
    guides$predictions <- predict(xgb_model, X)
    
    bins2 <- seq(min(guides$predictions), 
                 max(guides$predictions), 
                 length.out = input$binsscore + 1 )
    
    hist(guides$predictions, breaks = bins2, col = "lightblue", border = "grey",
         xlab = "XGB Score",
         main = "Distribtution of XGB Scores")
    
  })
  
}

shinyApp(ui = ui, server = server)
