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
  ## display total nr of candidate sgRNAs and mean GC
  output$qc_guide_count <- renderText({
    guides <- submit_reactive()

    paste("Total candidate sgRNAs: ", nrow(guides))
  })

  output$qc_mean_gc <- renderText({
    guides <- submit_reactive()

    gc <- (
      str_count(guides$Guide, "G")+
      str_count(guides$Guide, "C")
    ) / str_length(guides$Guide)
    paste0("Mean GC content: ", round(mean(gc) *100, 1),"%")
  })

  #sgRNA quality table
  output$guide_quality_table <- DT::renderDT({
    guides <- submit_reactive()

    #calculate GC content as precentage
    gc_precent<- (
      str_count(guides$Guide, "G")+
      str_count(guides$Guide, "C")
    )/ str_length(guides$Guide)*100

    #quality table
    quality_table<- data.frame(
      "Guide ID"=guides$ID,
      "Guide sequence"= guides$Guide,
      "GC (%)"= round(gc_precent,1),
      check.names = FALSE
    )

    #flag potential poly-T transcription termination motifs
    quality_table[["poly-T (TTTT)"]] <- ifelse(
      stringr::str_detect(guides$Guide, "TTTT"),
      "YES",
      "NO"
    )

    #connect the table to the slider
    #check whether each guide is within the selected GC range

    quality_table[["Within GC range"]]<- ifelse(
      quality_table[["GC (%)"]] >= input$gc_range[1]&
      quality_table[["GC (%)"]] <= input$gc_range[2],
      "YES",
      "NO"
    )

    #Detect homopolymer runs of A, C, or G (4+ nucleotides)
    other_runs<- stringr::str_extract(
      guides$Guide,
      "A{4,}|C{4,}|G{4,}"
    )
    quality_table[["Other homopolymer (4+)"]] <- ifelse(
      is.na(other_runs),
      "NO",
      other_runs
    )

    # combine basic sequence QC checks
    gc_outside <- (
      quality_table[["GC (%)"]] < input$gc_range[1] |
        quality_table[["GC (%)"]] > input$gc_range[2]
    )

    poly_t_warning <- stringr::str_detect(
      guides$Guide, "TTTT"
    )

    homopolymer_warning <- stringr::str_detect(
      guides$Guide, "A{4,}|C{4,}|G{4,}"
    )

    # create an overall sequence QC flag
    quality_table[["Sequence QC flag"]] <- ifelse(
      gc_outside | poly_t_warning | homopolymer_warning,
      "Review",
      "No flags"
    )
    ## Create an interactive sgRNA quality table

    dt<- DT::datatable(
      quality_table,
      rownames = FALSE,
      options = list(
        pageLength=5,
        lengthMenu= c(5,10,25,50),
        scrollY = "300px",
        scrollX=TRUE
      )
    )

    #highlight seq QC results
    DT::formatStyle(
      dt,
      "Sequence QC flag",
      backgroundColor = DT::styleEqual(
        c("No flags", "Review"),
        c("#EAF5EC", "#FFF3E0")
      ),
      color = DT::styleEqual(
        c("No flags", "Review"),
        c("#216E39","#9A5A08")
      )
    )

   })


  ## number and precentage of guides within selected GC range
  output$qc_gc_range_summary<- renderText({
    guides <-submit_reactive()

    gc<- (
      str_count(guides$Guide, "G")+
      str_count(guides$Guide, "C")
    )/ str_length(guides$Guide)

    gc_min<- input$gc_range[1]/100
    gc_max<- input$gc_range[2]/100

    n_within <- sum(gc>= gc_min & gc<= gc_max)
    precentage <- n_within/ length(gc) *100

    paste0(
      "Guides within selected range: ",
      n_within, " of ", length(gc),
      " (", round(precentage, 1), "%)"
    )
  })

  output$qcplot <- renderPlot({
    guides <- submit_reactive()
    guides$GC <- (str_count(guides$Guide, "G") +
                    str_count(guides$Guide, "C")
    )/str_length(guides$Guide)

    bins <- seq(min(guides$GC), max(guides$GC), length.out = 11)

    hist(guides$GC, breaks = bins, col = "#7FA8CC", border = "white",
         xlab = "GC-ratio",
         main = "GC-ratio Distribution")
    ## make his_plot respond to the GC slider
    rect(xleft = input$gc_range[1]/100,
         ybottom = 0,
         xright = input$gc_range[2]/100,
         ytop = par("usr")[4],
         col = adjustcolor("#F1CA81", alpha.f = 0.34),
         border = NA

         )
  #add vertical lines to the selected region
    abline(
      v= input$gc_range/100,
      col="#A76114",
      lty=2,
      lwd=2
    )
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
