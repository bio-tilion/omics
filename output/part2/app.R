library(shiny)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(DT)

load(file = "dge.RData")

df.plot.gene$gene_label <- ifelse(df.plot.gene$gene_name == "", df.plot.gene$ensembl_id, df.plot.gene$gene_name)

ui <- fluidPage(
  titlePanel("Interactive plots"),
  
  # first row for plots
  fluidRow(
    # Description
    column(
      3, 
      p("This interactive page generates two exploratory plots: a volcano plot showing the results of the differential gene expression analysis, and a boxplot showing the normalized counts for a single gene across all samples in the dataset."),
      p("At the bottom of the page there is a table showing the results of the analysis that can be used for searching for specific genes. The drop-down menu can be used to filter for ", em("Significance"), "category. The search text bar can be used to filter for any matching pattern."),
      p("Since the focus of the coursework was", em("SRRM4,"), "by default the boxplot shows its expression and the corresponding point in the volcano plot is highlighted and labelled. By clicking on a row of the table, you are selecting/deselecting the corresponding gene, which in turn updates the selected points in the volcano plot, as well as the boxplot. In case you selected multiple genes, the boxplot will always show the last gene in your selection, but the volcano plot will have all of those genes highlighted. If the selection is empty, the default", em("SRRM4"), "is shown again."),
    ),
    
    # Volcano plot
    column(6, plotOutput("volcano")),
    
    # Boxplot
    column(3, plotOutput("boxplot")) 
  ),
  
  # second row for table
  fluidRow(
    column(
      # width
      3,
      # list to choose from
      selectInput(
        "sig",
        "Significance",
        c("All", levels(df.plot.gene$sig))
      )
    )
  ),
  # new row for the table
  fluidRow(
    column(
      12, DT::dataTableOutput("table")
    )
  )
)

server <- function(input, output) {
  # Dataframe
  output$table <- DT::renderDataTable({
    df.shiny <- df.plot.gene
    if (input$sig != "All") {
      df.shiny <- df.plot.gene[df.plot.gene$sig == input$sig,]
    }
    df.shiny
  })
  
  # function to get the gene selection
  df.click <- function(input){
    if (length(input)){
      t <- input
      df.plot.gene[t, ]
    } else {
      df.plot.gene[df.plot.gene$gene_name == "SRRM4",]
    }
  }
  
  # raective function for gene selection
  get.df <- reactive({
    df.click(input$table_rows_selected)
  })
  
  # Volcano plot
  output$volcano <- renderPlot({
    p <- ggplot(df.plot.gene, aes(log2FoldChange, -log10(padj), col=sig)) +
      geom_point(size=0.5) + 
      scale_color_manual(values = c("red", "blue", "gray50", "gray80")) +
      theme_bw(base_size = 16)
    p$labels$colour <- c("Significance")

    df.temp <- get.df()
    p <- p + geom_point(data = df.temp, size=1.5, col="black", shape = 23) +
      scale_fill_manual(values = c("red", "blue", "gray50", "gray80")) +
      geom_text_repel(
        data = df.temp,
        aes(label = gene_label),
        col = "black"
    )
    p
  })
  
  # Boxplot
  output$boxplot <- renderPlot({
    # get gene name and ensembl id
    df.temp <- get.df()
    gene_ensembl <- tail(df.temp, 1)$ensembl_id
    gene_title <- tail(df.temp, 1)$gene_label
    
    gene_data <- plotCounts(dds_DGE, gene_ensembl, "diagnosis", returnData = TRUE)
    
    q <- ggplot(gene_data, aes(diagnosis, count)) +
      geom_boxplot() +
      geom_jitter(aes(col = diagnosis), alpha = 0.4, width = 0.1) +
      labs(
        title = gene_title,
        x = "Diagnosis",
        y = "Normalized count"
      ) +
      guides(col = "none") +
      theme_bw(base_size = 16) +
      theme(plot.title = element_text(hjust = 0.5))
    q
  })
}

shinyApp(ui = ui, server = server)
