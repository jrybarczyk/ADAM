#' @title GFAGpathUi
#' @description Launch GFAGpathUi Graphical User Interface (GUI) in local
#' machine or default browser.
#' @details This shiny (GUI) enables the visualization and plotting of
#' P.value distribution of the Paths in the GFAG output file. The plots
#' were optimized for a maximun of 50 Paths, you can still plot more than 50
#' paths but the graphs tend to get rather long.
#' @usage GFAGpathUi(browser)
#' @param browser is a logical variable necessary to run the app. When x=TRUE
#' the app is launched in your default web browser. When x=FALSE the app is
#' launched in your local machine.
#' @return return a shiny Graphical User Interface (GUI).
#' @export
#' @author Giordano Bruno Sanches Seco <giordano.bruno@unesp.br>
#' @examples
#' #Creating input files
#' data(DiffAedes)
#' data(ExpressionAedes)
#' data(GeneFunctionAedes)
#' data(ResultAnalysisAedes)
#' \dontrun{
#' #Input file 1
#' write.table(DiffAedes,file = "DiffAedes.txt",sep = "\t", col.names = TRUE,
#' row.names = FALSE, quote = FALSE)
#' #Input file 2
#' write.table(ExpressionAedes,file = "ExpressionAedes.txt",sep = "\t",
#' col.names = TRUE, row.names = FALSE, quote = FALSE)
#' #Input file 3
#' write.table(GeneFunctionAedes,file = "GeneFunctionAedes.txt",sep = "\t",
#' col.names = TRUE, row.names = FALSE, quote = FALSE)
#' #Input file 4
#' write.table(ResultAnalysisAedes,file = "ResultAnalysisAedes.txt",sep = "\t",
#' col.names = TRUE, row.names = FALSE, quote = FALSE)
#' #Grahphical analysis
#' GFAGpathUi(browser = TRUE) # Launch the app in your default web browser.
#' GFAGpathUi(browser = FALSE) # Launch the app in your local machine.
#' }

GFAGpathUi <-function(browser){
    tryCatch({
        if (browser){
            shiny::runApp(paste0(system.file("shiny",
                                            package="ADAM"),
                                            "/app.R"),
                                            launch.browser = TRUE)
        } else {
            shiny::runApp(paste0(system.file("shiny",
                                            package="ADAM"),
                                            "/app.R"))
        }
    },
    error = function(e) {print("Parameter must be TRUE or FALSE")
    }
    )}

