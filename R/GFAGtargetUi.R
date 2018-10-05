#' @title GFAGtargetUi
#' @description Launch Launch GFAGtargetUi Graphical User Interface (GUI) in
#' local machine or default browser
#' @importFrom dplyr mutate filter summarize arrange
#' @importFrom ggpubr ggarrange
#' @importFrom reshape2 melt
#' @importFrom DT dataTableOutput renderDataTable
#' @importFrom data.table data.table
#' @importFrom shiny actionButton actionLink addResourcePath column
#'         conditionalPanel downloadButton downloadHandler
#'         eventReactive fileInput fluidPage helpText isolate
#'         mainPanel need numericInput observe observeEvent
#'         outputOptions plotOutput radioButtons
#'         reactive reactiveValues renderPlot renderUI runApp
#'         selectInput shinyApp shinyServer shinyUI sidebarLayout
#'         sidebarPanel sliderInput stopApp tabPanel tabsetPanel
#'         textInput textOutput titlePanel uiOutput tags HTML
#'         h4 img icon updateTabsetPanel updateTextInput validate
#'         wellPanel checkboxInput br checkboxGroupInput a strong
#'         renderPrint fluidRow showNotification brushedPoints
#' @import ggsignif RColorBrewer colorRamps gridExtra knitr stringr stringi
#' @import varhandle DBI ggplot2 ggrepel testthat
#' @details This shiny (GUI) enables the creation of a plot that shows the
#' expression of the targets in a selected path between control and case.
#' See the vignette for more information.
#' @usage GFAGtargetUi(browser)
#' @param browser is a logical variable necessary to run the app.
#' When browser=TRUE the app is launched in your default web browser.
#' When browser=FALSE the app is launched in your local machine.
#' @return return a shiny Graphical User Interface (GUI).
#' @export
#' @author Giordano Bruno Sanches Seco <giordano.bruno@unesp.br>
#' @examples
#' #Creating input files
#' data(DiffHs)
#' data(ExpressionHs)
#' data(GeneFunctionHs)
#' data(ResultAnalysisHs)
#' \dontrun{
#' #Input file 1
#' write.table(DiffHs,file = "DiffHs.txt",sep = "\t", col.names = TRUE,
#' row.names = FALSE, quote = FALSE)
#' #Input file 2
#' write.table(ExpressionHs,file = "ExpressionHs.txt",sep = "\t",
#' col.names = TRUE, row.names = FALSE, quote = FALSE)
#' #Input file 3
#' write.table(GeneFunctionHs,file = "GeneFunctionHs.txt",sep = "\t",
#' col.names = TRUE, row.names = FALSE, quote = FALSE)
#' #Input file 4
#' write.table(ResultAnalysisHs,file = "ResultAnalysisHs.txt",sep = "\t",
#' col.names = TRUE, row.names = FALSE, quote = FALSE)
#' #Grahphical analysis
#' GFAGtargetUi(browser = TRUE) # Launch the app in your default web browser.
#' GFAGtargetUi(browser = FALSE) # Launch the app in your local machine.
#' }
GFAGtargetUi <-function(browser){
    tryCatch({
        if (browser){
            shiny::runApp(paste0(system.file("shiny",
                                package="ADAM"),
                                "/app2.R"),launch.browser = TRUE)
        } else {
            shiny::runApp(paste0(system.file("shiny",
                                package="ADAM"),"/app2.R"))
        }
    },
    error = function(e) {print("Parameter must be TRUE or FALSE")
    }
    )}
