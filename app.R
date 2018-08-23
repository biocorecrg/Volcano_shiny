### Shiny app for building Volcano plots
# Author: Sarah Bonnin - sarah.bonnin@crg.eu


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


# Define UI for application that draws a volcano plot
ui <- fluidPage(
   
   # Application title
   titlePanel("Volcano plot"),
   
   # Sidebar with options
   sidebarLayout(
      sidebarPanel(
        # Input file can be in text or csv format
        fileInput("file_volcano", "Upload text or csv file", accept = c("text/txt", "text/csv")),
        
        # Once the user file is loaded, the numeric columns are retrieved and put as options for the user to select the appropriate ones
        uiOutput("volcano_columns1"), # P-value column
        uiOutput("volcano_columns2"), # Fold change column
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
        # Set plot title and title size
        textInput("volcano_title", label = h4("Plot title"), value = "Volcano plot"),
        sliderInput("volcano_size_title",
                    "Title font size:",
                    min = 10, max = 60,
                    value = 2, step=2),
        textInput("volcano_xlab", label = h4("y-axis label"), value = "log2(fold-change)"),
        # X and y-axis labels and label size
        sliderInput("volcano_size_xlab",
                    "x-axis font size:",
                    min = 6, max = 40,
                    value = 6, step=2),
        textInput("volcano_ylab", label = h4("y-axis label"), value = "-log10(p-value)"),
        sliderInput("volcano_size_ylab",
                    "y-axis font size:",
                    min = 6, max = 40,
                    value = 6, step=2),
        # Colour chosen for down-, up- and non-regulated genes (selected according to above criteria)
        colourpicker::colourInput("volcano_color_up", "Color for up-regulated genes", value="red", allowTransparent = TRUE),
        colourpicker::colourInput("volcano_color_dw", "Color for down-regulated genes", value="green", allowTransparent = TRUE),
        colourpicker::colourInput("volcano_color_none", "Color for unchanged genes", value="black", allowTransparent = TRUE)
      ),
       # Show the plot when the button "Display Volcano plot" is clicked
      mainPanel(
        actionBttn("doVolcano", "Display Vocano plot", icon("bar-chart-o"), style = "jelly", size = "sm"),
         plotOutput("volcanoPlot"),
        conditionalPanel( # if the user clicked on "display Volcano plot", the download button also appears
          condition = "input.doVolcano == true",
          downloadBttn("downloadVolcano", "Save Volcano plot", size = "sm", style = "jelly")
        )
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # Reactive function to read in user input file
  input_volcano <- reactive({
    inFile <- input$file_volcano
    if (is.null(inFile))
      return(NULL)
    read.table(inFile$datapath, sep="\t", header=T, as.is = T)
    
  })
  shinyjs::enable("doVolcano")
  
  # Observe if input file was read and output numeric columns choices
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
      
      # Retrieve user input / data frame
      volcano1 <- input_volcano()
     
     if(is.null(input$volcano_pval_threshold) | is.null(input$volcano_fc_threshold))
       return(NULL)
     # rerieve limits for drawing lines and setting colours for up-, down- and non-regulated genes
     pval_threshold  <- input$volcano_pval_threshold
     volcano_fc_threshold_up <- input$volcano_fc_threshold[2]
     volcano_fc_threshold_dw <- input$volcano_fc_threshold[1]

      # Convert p-value to -log10 p-value
      volcano1[input$volcano_pval_col] <- -log10(volcano1[input$volcano_pval_col])
     
      # If the user chooses the option to color points according to down/up-regulated, add a column to the data frame
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
      # If the user clicks on the "Display Vocano" button, draw volcano plot
      output$volcanoPlot <- renderPlot({
        drawVolcano()
      })
      # The option to export Volcano plot appears at the same time as the plot itself
      output$downloadVolcano <- downloadHandler(
        filename = function(){paste(Sys.Date(), "_Volcano.png", sep ="")},
        content = function(file){
          png(file) # picking up correct file name in Chrome only...
          p <- drawVolcano()
          print(p)
          dev.off()
          #ggsave(filename=file, plot=p)
        }
      )
      
    })
    
   
}

# Run the application 
shinyApp(ui = ui, server = server)
# runApp(launch.browser=TRUE)

