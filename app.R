#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(colourpicker)
library(shinyjs)
library(shinyWidgets) # for downloadBttn

# test data frame for Volcano plot
#input_volcano <- data.frame(geneID=paste0(LETTERS, 1:2000), 
#                            pvalue=c(sample(seq(0,0.05, 0.000001), 200, replace=T), sample(seq(0,1, 0.0001), 1800, replace=T)),
#                            fc=c(sample(seq(-6.2, 8.3, 0.1), 800, replace=T), sample(seq(-0.3, 0.3, 0.1), 1200, replace=T)),
#                            stringsAsFactors = F)


# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Volcano plot"),
   
   # Sidebar with a slider input for the p-value threshold
   sidebarLayout(
      sidebarPanel(
        fileInput("file_volcano", "Upload text or csv file", accept = c("text/txt", "text/csv")),
        
        uiOutput("volcano_columns1"),
        uiOutput("volcano_columns2"),
        # P-value threshold for horizontal bar and points color
        sliderInput("volcano_pval_threshold",
                     "P-value threshold:",
                     min = 0, max = 0.1,
                     value = 0.05, step=0.005),
        # Fold change threshold for vertical bars and points color
         sliderInput("volcano_fc_threshold",
                     "fold change threshold:",
                     min = -10, max = 10,
                     value = c(-2,2), step=0.1),
        # Plot title and title size
        textInput("volcano_title", label = h4("Plot title"), value = "Volcano plot"),
        sliderInput("volcano_size_title",
                    "Title font size:",
                    min = 10, max = 60,
                    value = 2, step=2),
        textInput("volcano_xlab", label = h4("y-axis label"), value = "log2(fold-change)"),
        sliderInput("volcano_size_xlab",
                    "x-axis font size:",
                    min = 6, max = 40,
                    value = 6, step=2),
        textInput("volcano_ylab", label = h4("y-axis label"), value = "-log10(p-value)"),
        sliderInput("volcano_size_ylab",
                    "y-axis font size:",
                    min = 6, max = 40,
                    value = 6, step=2),
        colourpicker::colourInput("volcano_color_up", "Color for up-regulated genes", "red", allowTransparent = TRUE),
        colourpicker::colourInput("volcano_color_dw", "Color for down-regulated genes", "green", allowTransparent = TRUE),
        colourpicker::colourInput("volcano_color_none", "Color for unchanged genes", "black", allowTransparent = TRUE)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        actionBttn("doVolcano", "Display Vocano plot", icon("bar-chart-o"), style = "jelly", size = "sm"),
         plotOutput("volcanoPlot"),
         downloadBttn("downloadVolcano", "Save Volcano plot", size = "sm", style = "jelly")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  input_volcano <- reactive({
    inFile <- input$file_volcano
    if (is.null(inFile))
      return(NULL)
    read.table(inFile$datapath, sep="\t", header=T, as.is = T)
    
  })
  shinyjs::enable("doVolcano")
  
  observe({
    if(!is.null(input$file_volcano)){
      volcano1 <- input_volcano()
      
      output$volcano_columns1 <- renderUI({
        columns_volcano <- colnames(volcano1[unlist(lapply(volcano1, is.numeric))])
        selectInput(inputId="volcano_pval_col", "P-value column", choices=columns_volcano, selected=input$volcano_pval_col)
      })
      
      output$volcano_columns2 <- renderUI({
        columns_volcano <- colnames(volcano1[unlist(lapply(volcano1, is.numeric))])
        selectInput(inputId="volcano_fc_col", "Fold change column", choices=columns_volcano, selected=input$volcano_fc_col)
      })  
    }
  })
  
   #output$volcanoPlot <- renderPlot({
    drawVolcano <- reactive({
     # if no input file, do not output anything
     if (is.null(input$file_volcano))
       return(NULL)
      
      volcano1 <- input_volcano()
     
     # volcano1 <- input_volcano()
     # 
     # output$volcano_columns1 <- renderUI({
     #   columns_volcano <- colnames(volcano1[unlist(lapply(volcano1, is.numeric))])
     #   selectInput(inputId="volcano_pval_col", "P-value column", choices=columns_volcano, selected=input$volcano_pval_col)
     # })
     # 
     # output$volcano_columns2 <- renderUI({
     #   columns_volcano <- colnames(volcano1[unlist(lapply(volcano1, is.numeric))])
     #   selectInput(inputId="volcano_fc_col", "Fold change column", choices=columns_volcano, selected=input$volcano_fc_col)
     # })  
     
     if(is.null(input$volcano_pval_threshold) | is.null(input$volcano_fc_threshold))
       return(NULL)
     # rerieve limits for drawing lines and setting colours
     pval_threshold  <- input$volcano_pval_threshold
     volcano_fc_threshold_up <- input$volcano_fc_threshold[2]
     volcano_fc_threshold_dw <- input$volcano_fc_threshold[1]

      # Convert p-value to -log10 p-value
     volcano1[input$volcano_pval_col] <- -log10(volcano1[input$volcano_pval_col])

     
      ## If the user chooses the option to color points according to down/up-regulated, add a column to the data frame
       volcano1$select <- "Not regulated"
       volcano1$select[!is.na(volcano1[input$volcano_pval_col]) & volcano1[input$volcano_pval_col] >= -log10(pval_threshold) &
                         volcano1[input$volcano_fc_col] >= volcano_fc_threshold_up] <- "Up-regulated"
       volcano1$select[!is.na(volcano1[input$volcano_pval_col]) & volcano1[input$volcano_pval_col] >= -log10(pval_threshold) &
                         volcano1[input$volcano_fc_col] <= volcano_fc_threshold_dw] <- "Down-regulated"
       
        ggplot(volcano1, aes_string(x=input$volcano_fc_col, y=input$volcano_pval_col, color="select")) + 
          geom_point() + 
          theme_bw() + # remove grey background
          ggtitle(input$volcano_title) + # set plot title
          geom_vline(xintercept = c(volcano_fc_threshold_dw, volcano_fc_threshold_up), col="red") + # draw 2 vertical lines for positive and negative fold change thresholds
          geom_hline(yintercept = -log10(pval_threshold), col="red") + # draw 1 horizontal line for p-value
          scale_color_manual(values=c(input$volcano_color_dw, input$volcano_color_none, input$volcano_color_up)) + # set user-chosen colors for none regulated, up-regulated and down-regulated genes
          theme(plot.title = element_text(size = input$volcano_size_title, face = "bold", hjust=0.5), # size of volcano plot title
                axis.title.x=element_text(size=input$volcano_size_xlab), # size of x-axis label
                axis.title.y=element_text(size=input$volcano_size_ylab)) + # size of y-axis label
          labs(x=input$volcano_xlab, y=input$volcano_ylab) # user-chosen labels
   })
    
    observeEvent(input$doVolcano, {
      output$volcanoPlot <- renderPlot({
        drawVolcano()
      })
    })
    
   # export Volcano plot
   output$downloadVolcano <- downloadHandler(
     filename = function(){paste("volcano", input$download_type_heat, sep = ".")},
     content = function(file){
       p <- drawVolcano()
       export(p, file)
     }
   )
   
}

# Run the application 
shinyApp(ui = ui, server = server)

