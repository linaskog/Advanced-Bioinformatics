# Reading in Necessary Libraries ------------------------------------------

library(shiny)
library(bslib)
library(xgboost)

# takes gene ID as input, returns variables
source("~/GitHub/Advanced-Bioinformatics/R/app_get_variables.R")
source("~/GitHub/Advanced-Bioinformatics/R/app_input_panel.R")
source("~/GitHub/Advanced-Bioinformatics/R/app_quality_control_panel.R")
source("~/GitHub/Advanced-Bioinformatics/R/app_score_panel.R")

# Reading in Models -------------------------------------------------------

xgb_model <- xgb.load(
  "~/GitHub/Advanced-Bioinformatics/models/xgb_model.ubj")
linear_model <- readRDS(
  "~/GitHub/Advanced-Bioinformatics/models/linear_model.rds")

# Creating a User Interface -----------------------------------------------

# Creating a User Interface -----------------------------------------------

ui <- navbarPage(
  "sgRNA Analyzer",
  # Input Panel (R/app_input_panel.R)
  input_panel,
  # Quality Control Plot Panel (R/app_quality_control_panel.R)
  quality_control_panel,
  # Score Plot Panel (R/app_score_panel.R)
  score_panel
)


server <- function(input, output) {
  submit_reactive <- eventReactive ( input$submit, {
    guides <- get_variables(geneID = input$geneID)
  })

  output$output_text <- renderText({
    guides <- submit_reactive()
    paste0(
      attr(guides,"gene_symbol"),
      " (", attr(guides, "ensembl_id"), ") loaded, ",
      nrow(guides), " sgRNAs found"
    )
    })
  ## display total nr of candidate sgRNAs and mean GC
  output$qc_guide_count <- renderText({
    guides <- submit_reactive()
    paste("Total candidate sgRNAs: ", nrow(guides))
  })

  output$qc_mean_gc <- renderText({
    guides <- submit_reactive()
    paste0("Mean GC content: ", round(mean(guides$GC) *100, 1),"%")
  })

  #sgRNA quality table
  output$guide_quality_table <- DT::renderDT({
    guides <- submit_reactive()

    #calculate GC content as precentage
    gc_precent<- guides$GC*100

    #quality table
    quality_table<- data.frame(
      "Guide ID"=guides$UID,
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

    gc<- guides$GC

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

  evaluated_scores <- reactive({
    guides <- submit_reactive()
    
    pos <- as.matrix(guides[, grep("pos", colnames(guides))])
    additional_data <- as.matrix(guides[, c("GC", "strand", "chromosome", 
                                            "region_type", "gene_overlap_count", 
                                            "H3K27ac_peak_overlap")
    ]
    )
    xgb_data <- cbind(pos, additional_data)
    
    scores <- data.frame("UID" = guides$UID,
                         "Guide" = guides$Guide,
                         "LinearScore" = predict(linear_model, guides),
                         "XGBScore" = predict(xgb_model, xgb_data),
                         check.names = FALSE)
    
    weight_linear <- input$linear_weight
    weight_xgb <- 1 - weight_linear
    
    # Computing combined model score
    scores$Combination <- (weight_linear*scores$LinearScore) + 
      (weight_xgb*scores$XGBScore)
    
    # Computing difference between model scores
    scores$Difference <- abs(100 * (scores$LinearScore - scores$XGBScore)/
                               (scores$LinearScore + scores$XGBScore)/2)
    
    # Ordering sgRNAs for selected model
    scores[order(scores[[input$select]]), ]
  })
  
  output$plot <- plotly::renderPlotly({
    scores <- evaluated_scores()
    top_scores <- head(scores, min(input$top_n, nrow(scores)))
    
    p <- ggplot(
      data = scores, 
      aes(
        x = LinearScore,
        y = XGBScore,
        colour = Combination,
        text = paste0("<b>sgRNA: </b>", UID, "<br>",
                      "<b>Guide: </b>", Guide, "<br>", 
                      "<b>Combination: </b>", round(Combination, 4), "<br>", 
                      "<b>Linear: </b>", round(LinearScore, 4), "<br>",
                      "<b>XGBoost: </b>", round(XGBScore, 4), "<br>",
                      "<b>Difference: </b>", round(Difference, 2), "<br>")
      )) + 
      geom_point(alpha = 0.3) +
      geom_point(
        data = top_scores,
        aes(x = LinearScore,
            y = XGBScore),
        colour = "maroon") +
      labs( x = "Linear Model Prediction",
            y = "XGBoost Model Prediction")+
      theme_light()
    
    plotly::ggplotly(p, tooltip = "text")
  })
  output$top_sgrna_table <- DT::renderDT({
    scores <- evaluated_scores()
    scores <- head(scores, min(input$top_n, nrow(scores)))
    
    score_table <- 
      data.frame(
        "ID" = scores$UID, 
        "Guide" = scores$Guide,
        "Linear Score" = round(scores$LinearScore, 4),
        "XGBoost Score" = round(scores$XGBScore, 4),
        "Combined Score" = round(scores$Combination, 4),
        "Difference %" = round(scores$Difference, 2),
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

shinyApp(ui = ui, server = server)
