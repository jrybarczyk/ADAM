suppressMessages(require(shiny))
suppressMessages(require(shinydashboard))
suppressMessages(require(shinyjs))
suppressMessages(require(stringr))
suppressMessages(require(stringi))
suppressMessages(require(varhandle))
suppressMessages(require(DBI))
suppressMessages(require(GO.db))
suppressMessages(require(ggplot2))
suppressMessages(require(ggrepel))
suppressMessages(require(ggpubr))
suppressMessages(require(ggsignif))
suppressMessages(require(reshape2))
suppressMessages(require(RColorBrewer))
suppressMessages(require(colorRamps))
suppressMessages(require(plyr))
suppressMessages(require(dplyr))
suppressMessages(require(DT))
suppressMessages(require(data.table))
suppressMessages(require(gridExtra))


#xx <- as.list(GOTERM)



ui <- fluidPage(
    shinyjs::useShinyjs(),
    br(),
    sidebarLayout(
        sidebarPanel(
            tabsetPanel(id='tabs',
                tabPanel('File Upload & Case Selection',value = 1,
                    br(),
                    p('GFAG Data',style='font-size: 25px; font-style: bold;
                        color: black;'),
                    uiOutput('GFAGui'),
                    tags$hr(style='border: 1px #8d8e90 solid;'),
                    p('Expression Data',style='font-size: 25px;
                        font-style: bold; color: black;'),
                    uiOutput('NEXui'),
                    tags$hr(style='border: 1px #8d8e90 solid;'),
                    p('Path-to-Target Data',style='font-size: 25px;
                        font-style: bold; color: black;'),
                    uiOutput('Pathui1'),
                    tags$hr(style='border: 1px #8d8e90 solid;'),
                    p('Diferential Expression Data',style='font-size: 25px;
                        font-style: bold;color: black;'),
                    uiOutput('DEAui'),
                    tags$hr(style='border: 1px #8d8e90 solid;'),
                    actionButton(inputId = 'Reset',
                        label = 'Reset Selected files/columns'),
                    helpText('Click in the buttom above to reset
                        the selected files and erase updated columns')
                ),
                tabPanel('GFAG Data Filter',value = 2,
                    br(),
                    uiOutput('NN1'),
                    br(),
                    uiOutput('NN2'),
                    br(),
                    uiOutput('BT1'),
                    br(),
                    uiOutput('BT2'),
                    br(),
                    uiOutput('BT3'),
                    br(),
                    uiOutput('BT4'),
                    tags$hr(),
                    uiOutput('golist'),
                    #tags$hr(),
                    #uiOutput('golist2'),
                    uiOutput('BT5'),
                    helpText('Click in the buttom above to update the value of
                        the selected GFAG/Path.'),
                    tags$hr(),
                    uiOutput('Plot_to_start'),
                    helpText('Click in the buttom above to start the plotting
                        phase. Note that this buttom will only be active if
                        currently there is a GFAG/Path selected.')),
                tabPanel('Plot',value = 3, br(),
                    uiOutput('qval_sld'),
                    br(),
                    uiOutput('logFC_sld'),
                    hr(),
                    uiOutput('height'),
                    uiOutput('width'),
                    br(),
                    uiOutput('plotfont1'),
                    br(),
                    uiOutput('fontsize1'),
                    br(),
                    uiOutput('fontsize2'),
                    br(),
                    uiOutput('fontsize3'),
                    br(),
                    uiOutput('fontsize4'),
                    br(),
                    uiOutput('fontface1'),
                    br(),
                    uiOutput('fontface2'),
                    br(),
                    uiOutput('fontface3'),
                    hr(),
                    uiOutput('logfcbt'),
                    br(),
                    uiOutput('qvalbt'),
                    br(),
                    uiOutput('qvlgbt'),
                    br(),
                    uiOutput('nfbt'),
                    tags$hr(),
                    downloadButton(outputId = 'fplot',
                        label = 'Download Plot')))),
        mainPanel(conditionalPanel(condition = 'input.tabs==1',
            tabsetPanel(
                tabPanel('Summary',
                    hr(),
                    textOutput('GFAG_text11'),
                    textOutput('GFAG_text21'),
                    textOutput('GFAG_text22'),
                    textOutput('GFAG_text31'),
                    hr(),
                    textOutput('NEX_text11'),
                    textOutput('NEX_text21'),
                    textOutput('NEX_text31'),
                    textOutput('NEX_text32'),
                    textOutput('NEX_text41'),
                    textOutput('NEX_text51'),
                    hr(),
                    textOutput('Path_text11'),
                    textOutput('Path_text21'),
                    textOutput('Path_text31'),
                    textOutput('Path_text40'),
                    textOutput('Path_text41'),
                    textOutput('Path_text42'),
                    textOutput('Path_text43'),
                    textOutput('Path_text44'),
                    textOutput('Path_text51'),
                    hr(),
                    textOutput('DEA_text11'),
                    textOutput('DEA_text21'),
                    textOutput('DEA_text31'),
                    textOutput('DEA_text41'),
                    textOutput('DEA_text51'),
                    textOutput('DEA_text52'),
                    textOutput('DEA_text53'),
                    textOutput('DEA_text54'),
                    hr()),
                tabPanel('Tables',
                    h3('GFAG Data'),
                    DT::dataTableOutput('GFAG_table'),
                    tags$hr(style='border: 1px #8d8e90 solid;'),
                    h3('Expression Data'),
                    DT::dataTableOutput('NEX_table'),
                    tags$hr(style='border: 1px #8d8e90 solid;'),
                    h3('Path-to-Target Data'),
                    DT::dataTableOutput('Path_table'),
                    tags$hr(style='border: 1px #8d8e90 solid;'),
                    h3('Differential Expression Analisys Data'),
                    DT::dataTableOutput('DEA_table'),
                    tags$hr(style='border: 1px #8d8e90 solid;')))),
        conditionalPanel(condition = 'input.tabs==2',
            h3('GFAG Data:'),
            textOutput('IDAT4'),
            DT::dataTableOutput('ftable2')),
        conditionalPanel(condition = 'input.tabs==3',
            h3('Your Plot:'),
            uiOutput('crplot')))

    )
)

server <- function(input, output, session) {

    shinyjs::useShinyjs()
    options(shiny.maxRequestSize=30*1024^2)
    rv<-reactiveValues()


    shinyjs::disable('upGene1')
    shinyjs::disable('upcase')
    shinyjs::disable('upcontrol')
    shinyjs::disable('upgfag')
    shinyjs::disable('upGene2')
    shinyjs::disable('uplog')
    shinyjs::disable('uppval')
    shinyjs::disable('ViewH')
    shinyjs::disable('ViewN')
    shinyjs::disable('ViewNH')
    shinyjs::disable('ViewR')
    shinyjs::disable('upfiltered')
    shinyjs::disable('plot_start')
    shinyjs::disable('upPath')
    shinyjs::disable('upGene')
    shinyjs::disable('fplot')


    observeEvent(input$NEXdata,{
        rv$Case_index<-NULL
        rv$Control_index<-NULL
        file1 <- input$NEXdata
        if(is.null(file1)){
            rv$NEX_data<-NULL
            return()
        }
            else{
                a<-read.delim(file=file1$datapath,
                    sep='\t',
                    header = TRUE,
                    stringsAsFactors = TRUE)
      if ( (is.null(rv$GFAG_data) & is.null(a)) | (!is.null(rv$GFAG_data)
            & is.null(a)) | ((!(is.null(rv$GFAG_data)) & !is.null(a)) &
            (length(colnames(a))>=3) &
            (TRUE %in% grepl(rv$Controlid,colnames(a)))
            & (TRUE %in% grepl(rv$Caseid,colnames(a)))
            & (TRUE %in% grepl('GO:00',a[,grep('ID',colnames(a),value = T)])))
            | (((!(is.null(rv$GFAG_data))) & (!is.null(a)) ) &
            (!(length(colnames(a))>=3) |
            !(TRUE %in% grepl(rv$Controlid,colnames(a))) |
            !(TRUE %in% grepl(rv$Caseid,colnames(a)))))) {
                rv$NEX_data1<-a
                rv$NEX_data2<-NULL
      }
                else{
                    if (((!(is.null(rv$GFAG_data)) & !is.null(a)) &
                        (length(colnames(a))>=3) &
                        (TRUE %in% grepl(rv$Controlid,colnames(a))) &
                        (TRUE %in% grepl(rv$Caseid,colnames(a))) &
                        !(TRUE %in% grepl('GO:00',a[,grep('ID',
                            colnames(a),value = T)])))) {
                            rv$NEX_data1<-a
                            rv$NEX_data2<-a
                    }
                }
            }

    })


    output$NEX_text1<-renderText({
        if((is.null(rv$GFAG_data2)) & is.null(rv$NEX_data1)){
            return('Warning: You need to load a GFAG output file first!')
        }
    })


    output$NEX_text11<-renderText({
        if((is.null(rv$GFAG_data2)) & is.null(rv$NEX_data1)){
            return(paste0('The uploaded expression file contains:
                No expression file was uploaded!'))
        }
    })


    output$NEX_text2<-renderText({
        if(!(is.null(rv$GFAG_data2)) & is.null(rv$NEX_data1)){
            return('Warning: No expression file was uploaded!')
        }
   })


    output$NEX_text21<-renderText({
        if(!(is.null(rv$GFAG_data2)) & is.null(rv$NEX_data1)){
            return(paste0('The uploaded expression file contains:
                No expression file was uploaded!'))
       }
    })


    output$NEX_text3<-renderText({
        if (((!(is.null(rv$GFAG_data)) & !is.null(rv$NEX_data1))
            & (length(colnames(rv$NEX_data1))>=3) &
            (TRUE %in% grepl(rv$Controlid,colnames(rv$NEX_data1))) &
            (TRUE %in% grepl(rv$Caseid,colnames(rv$NEX_data1))) &
            !(TRUE %in% grepl('GO:00',rv$NEX_data1[,grep('ID',
                colnames(rv$NEX_data1),value = T)])))) {
                return('Expression file sucessfully uploaded!')
        }
    })


    output$NEX_text31<-renderText({
        req(input$gfagdata)
        req(input$NEXdata)
        if (((!(is.null(rv$GFAG_data)) & !is.null(rv$NEX_data1)) &
            (length(colnames(rv$NEX_data1))>=3) &
            (TRUE %in% grepl(rv$Controlid,colnames(rv$NEX_data1))) &
            (TRUE %in% grepl(rv$Caseid,colnames(rv$NEX_data1))) &
            !(TRUE %in% grepl('GO:00',rv$NEX_data1[,grep(
                 'ID',colnames(rv$NEX_data1),value = T)])))) {
                return(paste0('The uploaded expression file contains:
                    A total of ',length(rownames(rv$NEX_data2)),' entries.'))
        }
    })


    output$NEX_text32<-renderText({
        req(input$gfagdata)
        req(input$NEXdata)
        if (((!(is.null(rv$GFAG_data)) & !is.null(rv$NEX_data1)) &
            (length(colnames(rv$NEX_data1))>=3) &
            (TRUE %in% grepl(rv$Controlid,colnames(rv$NEX_data1))) &
            (TRUE %in% grepl(rv$Caseid,colnames(rv$NEX_data1))) &
            !(TRUE %in% grepl('GO:00',rv$NEX_data1[,grep('ID',
                colnames(rv$NEX_data1),value = T)])))) {
                return(paste0('The uploaded expression file contains: ',
                    length(unique(rv$NEX_data2[,1])),' unique targets.'))
        }
    })


    output$NEX_text4<-renderText({
        req(input$gfagdata)
        req(input$NEXdata)
        if (((!(is.null(rv$GFAG_data)) & !is.null(rv$NEX_data1)) &
            (length(colnames(rv$NEX_data1))>=3) &
            (TRUE %in% grepl(rv$Controlid,colnames(rv$NEX_data1)))
            & (TRUE %in% grepl(rv$Caseid,colnames(rv$NEX_data1))) &
            (TRUE %in% grepl('GO:00',rv$NEX_data1[,grep('ID',
                colnames(rv$NEX_data1),value = T)])))) {
                return('Error: Failed to upload diferential expression file!
                    The selected file may not contain the columns of the
                    same Case/Control contemplated in the GFAG output file
                    or does not have the minimum amount of required
                    columns=3 (Gene/Probe,Control,Case). Please remember that
                    the GFAG output file should not be altered!')
        }
    })


    output$NEX_text41<-renderText({
        req(input$gfagdata)
        req(input$NEXdata)
        if (((!(is.null(rv$GFAG_data)) & !is.null(rv$NEX_data1)) &
            (length(colnames(rv$NEX_data1))>=3) &
            (TRUE %in% grepl(rv$Controlid,colnames(rv$NEX_data1))) &
            (TRUE %in% grepl(rv$Caseid,colnames(rv$NEX_data1))) &
            (TRUE %in% grepl('GO:00',rv$NEX_data1[,grep('ID',
                colnames(rv$NEX_data1),value = T)])))) {
                return(paste0('The uploaded expression file contains:
                    No expression file was uploaded!'))
        }
    })


    output$NEX_text5<-renderText({
        req(input$gfagdata)
        req(input$NEXdata)
        if ((((!(is.null(rv$GFAG_data))) & (!is.null(rv$NEX_data1))) &
            (!(length(colnames(rv$NEX_data1))>=3) |
            !(TRUE %in% grepl(rv$Controlid,colnames(rv$NEX_data1))) |
            !(TRUE %in% grepl(rv$Caseid,colnames(rv$NEX_data1)))))) {
                return('Error: Failed to upload diferential expression file!
                    The selected file may not contain the columns of the same
                    Case/Control contemplated in the GFAG output file or does
                    not have the minimum amount of required columns=3
                    (Gene/Probe,Control,Case). Please remember that the GFAG
                    output file should not be altered!')
        }
    })


    output$NEX_text51<-renderText({
        req(input$gfagdata)
        req(input$NEXdata)
        if ((((!(is.null(rv$GFAG_data))) & (!is.null(rv$NEX_data1))) &
            (!(length(colnames(rv$NEX_data1))>=3) |
            !(TRUE %in% grepl(rv$Controlid,colnames(rv$NEX_data1))) |
            !(TRUE %in% grepl(rv$Caseid,colnames(rv$NEX_data1)))))) {
                return(paste0('The uploaded expression file contains:
                    No expression file was uploaded!'))
        }
    })


    output$NEX_table<-DT::renderDataTable({
        if(!(is.null(rv$NEX_data2))){
            rv$NEX_data2
        }
            else{
                return()
            }
    })


    output$NEXui<-renderUI({
        if ((is.null(rv$GFAG_data2)) & is.null(rv$NEX_data1)) {
            shinyjs::disabled(wellPanel(
                style = "background-color: #c6c6c6;border: 1px #c6c6c6 solid;",
                textOutput('NEX_text1'),
                tags$head(tags$style("#NEX_text1{color: blue;font-size: 15px;
                    font-style: italic;}"))))
        }
            else{
                if ((!(is.null(rv$GFAG_data2))) & is.null(rv$NEX_data1)) {
                    wellPanel(
                        style = "background-color: #ffffff;border:
                             1px #86e6fc solid;",
                        fileInput(label = 'Upload Expression Data:',
                                  inputId = 'NEXdata'),
                    textOutput('NEX_text2'),
                    tags$head(tags$style("#NEX_text2{color: blue;
                        font-size: 15px; font-style: italic;}")))
            }
                else{
                    if ((((!(is.null(rv$GFAG_data))) &
                        (!is.null(rv$NEX_data1))) &
                        (!(length(colnames(rv$NEX_data1))>=3) |
                        !(TRUE %in% grepl(rv$Controlid,colnames(rv$NEX_data1)))
                        | !(TRUE %in% grepl(rv$Caseid,
                            colnames(rv$NEX_data1))))))  {
                                wellPanel(style = "background-color: #ffffff;
                                    border: 1px #86e6fc solid;",
                                        fileInput(
                                            label =
                                                'Upload Expression Data:',
                                            inputId = 'NEXdata'),
                                            textOutput('NEX_text5'),
                                            tags$head(
                                                tags$style(
                                                    "#NEX_text5{color: red;
                                                        font-size: 15px;
                                                        font-style: italic;
                                }")))
                    }
                        else{
                            if((((!(is.null(rv$GFAG_data))) &
                                (!is.null(rv$NEX_data1))) &
                                (length(colnames(rv$NEX_data1))>=3) &
                                (TRUE %in% grepl(rv$Controlid,
                                    colnames(rv$NEX_data1))) &
                                (TRUE %in% grepl(rv$Caseid,
                                    colnames(rv$NEX_data1))) &
                                (! (TRUE %in% grepl('GO:00',
                                    rv$NEX_data1[,grep('ID',
                                        colnames(rv$NEX_data1),
                                        value = T)]))))) {
                                        shinyjs::disabled(wellPanel(
                                            style = "background-color:
                                                #c6c6c6;border:
                                                1px #c6c6c6 solid;",
                                            textOutput('NEX_text3'),
                                            tags$head(tags$style(
                                                "#NEX_text3{color: green;
                                                font-size: 15px;
                                                font-style: italic;
                                        }"))))
                                }
                                    else{
                                        if(((!(is.null(rv$GFAG_data)) &
                                             !is.null(rv$NEX_data1)) &
                                            (length(colnames(rv$NEX_data1))>=3)
                                            & (TRUE %in% grepl(rv$Controlid,
                                                colnames(rv$NEX_data1))) &
                                            (TRUE %in% grepl(rv$Caseid,
                                                colnames(rv$NEX_data1))) &
                                            (TRUE %in% grepl('GO:00',
                                                rv$NEX_data1[,grep('ID',
                                                colnames(rv$NEX_data1),
                                                value = T)])))) {
                                                    wellPanel(style =
                                                    "background-color:
                                                    #ffffff;border:
                                                    1px #86e6fc solid;",
                                                    fileInput(label =
                                                        'Upload Expression
                                                        Data:',
                                                    inputId = 'NEXdata'),
                                                    textOutput('NEX_text4'),
                                                    tags$head(tags$style(
                                                        "#NEX_text4{color:
                                                        red; font-size: 15px;
                                                        font-style: italic;
                                                    }")))
                                        }
                                    }
                        }
                }
            }
    })


    output$selectGene1<-renderUI({
        if (is.null(rv$NEX_data2)) {
            shinyjs::disabled(selectInput(label = 'Select a Target column:',
                inputId = 'Genesel1',choices = c()))
        }
            else{
                if(!(is.null(rv$NEX_data2))){
                    shinyjs::enable('upGene1')
                    selectInput(label = 'Select a Target column:',
                        inputId = 'Genesel1',
                        choices = colnames(rv$NEX_data2))
                }
            }
    })


    observeEvent(input$upcase,{
        rv$Gene1_index<-grep(paste0('^',input$Genesel1,'$'),
            colnames(rv$NEX_data2))
    })


    output$selectcase<-renderUI({
        if (is.null(rv$NEX_data2)) {
            shinyjs::disabled(selectInput(label = 'Select a case:',
                inputId = 'Casesel',choices = c()))
        }
            else{
                if(!(is.null(rv$NEX_data2))){
                    shinyjs::enable('upcase')
                    selectInput(label = 'Select a case:',inputId = 'Casesel',
                        choices = colnames(rv$NEX_data2))
                }
            }
    })


    observeEvent(input$upcase,{
        rv$Case_index<-grep(paste0('^',input$Casesel,'$'),
            colnames(rv$NEX_data2))
      })


    output$selectcontrol<-renderUI({
        if (is.null(rv$NEX_data2)) {
            shinyjs::disabled(selectInput(label = 'Select a control:',
                inputId = 'Controlsel',choices = c()))
        }
            else{
                if(!(is.null(rv$NEX_data2))){
                    shinyjs::enable('upcontrol')
                    selectInput(label = 'Select a control:',
                        inputId = 'Controlsel',
                        choices = colnames(rv$NEX_data2))
                }
            }
    })


    observeEvent(input$upcontrol,{
        rv$Control_index<-grep(paste0('^',input$Controlsel,'$'),
            colnames(rv$NEX_data2))
    })


    observeEvent(input$DEAdata,{
        rv$logsel_index<-NULL
        rv$pvalsel_index<-NULL
        rv$colselsum<-0
        rv$genekey<-0
        rv$logkey<-0
        rv$pvalkey<-0
        file2 <- input$DEAdata
        if(is.null(file2)){
            rv$DEA_data<-NULL
        }
            else{
                c<-read.delim(file=file2$datapath,
                    sep=input$sepDEA,
                    header = input$headerDEA,
                    stringsAsFactors = input$stringAsFactorsDEA)
                d<-c
                rv$DEA_data1<-c
                rv$DEA_data2<-d
            }

    })


    observeEvent(input$headerDEA,{
        rv$logsel_index<-NULL
        rv$pvalsel_index<-NULL
        rv$colselsum<-0
        rv$genekey<-0
        rv$logkey<-0
        rv$pvalkey<-0
        file2 <- input$DEAdata
        if(is.null(file2)){
            rv$DEA_data<-NULL
        }
            else{
                c<-read.delim(file=file2$datapath,
                    sep=input$sepDEA,
                    header = input$headerDEA,
                    stringsAsFactors = input$stringAsFactorsDEA)
                    d<-c
                    rv$DEA_data2<-d
           }

    })


    observeEvent(input$stringAsFactorsDEA,{
        rv$logsel_index<-NULL
        rv$pvalsel_index<-NULL
        rv$colselsum<-0
        rv$genekey<-0
        rv$logkey<-0
        rv$pvalkey<-0
        file2 <- input$DEAdata
        if(is.null(file2)){
        rv$DEA_data<-NULL
        }
            else{
                c<-read.delim(file=file2$datapath,
                    sep=input$sepDEA,
                    header = input$headerDEA,
                    stringsAsFactors = input$stringAsFactorsDEA)
                    d<-c
                    rv$DEA_data2<-d
            }
    })


    observeEvent(input$sepDEA,{
        rv$logsel_index<-NULL
        rv$pvalsel_index<-NULL
        rv$colselsum<-0
        rv$genekey<-0
        rv$logkey<-0
        rv$pvalkey<-0
        file2 <- input$DEAdata
        if(is.null(file2)){
            rv$DEA_data1<-NULL
        }
            else{
                c<-read.delim(file=file2$datapath,
                    sep=input$sepDEA,
                    header = input$headerDEA,
                    stringsAsFactors = input$stringAsFactorsDEA)
                    d<-c
                    rv$DEA_data2<-d
            }

    })


    DEA_data<-reactive({
        file2 <- input$DEAdata
        if(is.null(file2)){
            return()
        }
            else{
                a<-read.delim(file=file2$datapath,
                    sep=input$sepDEA,
                    header = input$headerDEA,
                    stringsAsFactors = input$stringAsFactorsDEA)
                    b<-a
                    b
            }
    })


    output$DEA_text1<-renderText({
        if(((((is.null(rv$GFAG_data)) & (is.null(rv$NEX_data1))) &
            is.null(rv$Path_data1))) & is.null(rv$DEA_data1)){
            return('Warning: You need to load a GFAG output file first!')
        }
            else{
                NULL
            }
    })


    output$DEA_text11<-renderText({
        if(((((is.null(rv$GFAG_data)) &
            (is.null(rv$NEX_data1))) & is.null(rv$Path_data1))) &
            is.null(rv$DEA_data1)){
            return('The uploaded differential expression file contains:
                no diferential expression file was uploded!')
        }
            else{
                NULL
            }
    })


    output$DEA_text2<-renderText({
        if(((((!is.null(rv$GFAG_data)) &
            (is.null(rv$NEX_data1))) & is.null(rv$Path_data1))) &
            is.null(rv$DEA_data1)){
            return('Warning: You need to load a expression file first!')
        }
            else{
                NULL
            }
    })


    output$DEA_text21<-renderText({
        if(((((!is.null(rv$GFAG_data)) & (is.null(rv$NEX_data1))) &
        is.null(rv$Path_data1))) &
        is.null(rv$DEA_data1)){
            return('The uploaded differential expression file contains:
                no diferential expression file was uploded!')
        }
            else{
                NULL
            }
    })


    output$DEA_text3<-renderText({
        if(((((!is.null(rv$GFAG_data)) &
            !(is.null(rv$NEX_data1))) & is.null(rv$Path_data1))) &
            is.null(rv$DEA_data1)){
            return('Warning: You need to load a path-to-Target file first!')
        }
            else{
                NULL
            }
    })


    output$DEA_text31<-renderText({
        if(((((!is.null(rv$GFAG_data)) &
            !(is.null(rv$NEX_data1))) & is.null(rv$Path_data1))) &
            is.null(rv$DEA_data1)){
            return('The uploaded differential expression file contains:
                no diferential expression file was uploded!')
        }
            else{
                NULL
            }
    })


    output$DEA_text4<-renderText({
        if(((((!is.null(rv$GFAG_data)) &
            !(is.null(rv$NEX_data1))) & !is.null(rv$Path_data1))) &
            is.null(rv$DEA_data1)){
            return('Warning: No diferential expression file was uploaded!')
        }
            else{
                NULL
            }
    })


    output$DEA_text41<-renderText({
        if(((((!is.null(rv$GFAG_data)) & !(is.null(rv$NEX_data1))) &
            !is.null(rv$Path_data1))) & is.null(rv$DEA_data1)){
            return('The uploaded differential expression file contains:
                no diferential expression file was uploded!')
        }
            else{
                NULL
            }
    })


    output$DEA_text5<-renderText({
        if(((((!is.null(rv$GFAG_data)) &
            !(is.null(rv$NEX_data1))) & !is.null(rv$Path_data1))) &
            !is.null(rv$DEA_data1) &
            (rv$genekey==0 | rv$logkey==0 | rv$pvalkey==0)){
            return('Diferential expression file sucessfully uploaded!
                Please proced with the selection of table parameters
                and columns!')
        }
            else{
                NULL
            }
    })


    output$DEA_text51<-renderText({
        req(input$gfagdata)
        req(input$NEXdata)
        req(input$pathgene)
        req(input$DEAdata)
        if(((((!is.null(rv$GFAG_data)) & !(is.null(rv$NEX_data1))) &
            !is.null(rv$Path_data1))) & !is.null(rv$DEA_data1)){
            return(paste0('The uploaded differential expression file contains:
                ',length(rv$DEA_data1[,1]), ' entries!'))
        }
            else{
                NULL
            }
    })


    output$DEA_text52<-renderText({
        req(input$gfagdata)
        req(input$NEXdata)
        req(input$pathgene)
        req(input$DEAdata)
        if(((((!is.null(rv$GFAG_data)) & !(is.null(rv$NEX_data1))) &
            !is.null(rv$Path_data1))) & !is.null(rv$DEA_data1)){
            return(paste0('The uploaded differential expression file contains:
                ',length(unique(rv$DEA_data1[,1])), ' unique targets!'))
        }
            else{
                NULL
            }
    })


    output$DEA_text53<-renderText({
        req(input$gfagdata)
        req(input$NEXdata)
        req(input$pathgene)
        req(input$DEAdata)
        if(((((!is.null(rv$GFAG_data)) & !(is.null(rv$NEX_data1))) &
            !is.null(rv$Path_data1))) & !is.null(rv$DEA_data1)){
            if (length(intersect(unique(rv$NEX_data2[,1]),rv$DEA_data2[,1]))==
                length(unique(rv$DEA_data2[,1]))) {
                return(paste0('The uploaded differential expression
                    file contains: ',length(intersect(unique(rv$NEX_data2[,1]),
                    rv$DEA_data2[,1])),' targets in common with expression
                    data. All targets in the differential expression file are
                    contemplated in the expression file!'))
            }
                else{
                    if ((length(intersect(unique(rv$NEX_data2[,1]),
                        rv$DEA_data2[,1]))-length(unique(rv$DEA_data2[,1])))<0)
                        {
                        return(paste0('The uploaded differential expression
                            file contains: ',
                            length(intersect(unique(rv$NEX_data2[,1]),
                            rv$DEA_data2[,1])),
                            ' targets in common with expression data.
                            There are ',
                            abs(length(intersect(unique(rv$NEX_data2[,1]),
                                rv$DEA_data2[,1]))-
                                length(unique(rv$DEA_data2[,1]))),
                            ' targets in the differential expression file not
                            contemplated in your expression file!'))
                        }
                            else{
                                if ((length(intersect(unique(rv$NEX_data2[,1]),
                                    rv$DEA_data2[,1]))-
                                    length(unique(rv$DEA_data2[,1])))>0) {
                                    return(paste0('The uploaded differential
                                        expression file contains: ',
                                    length(intersect(unique(rv$NEX_data2[,1]),
                                        rv$DEA_data2[,1])),
                                    ' targets in common with expression data.
                                    There are ',
                                    abs(length(intersect(unique(
                                    rv$NEX_data2[,1]),rv$DEA_data2[,1]))-
                                    length(unique(rv$DEA_data2[,1]))),
                                    ' targets in the expression file not
                                    contemplated in your differential
                                    expression file!'))
                                }
                            }
                }
            }
                else{
                    NULL
                }
    })


    output$DEA_text54<-renderText({
        req(input$gfagdata)
        req(input$NEXdata)
        req(input$pathgene)
        req(input$DEAdata)
        if(((((!is.null(rv$GFAG_data)) & !(is.null(rv$NEX_data1))) &
            !is.null(rv$Path_data1))) & !is.null(rv$DEA_data1)){
            if (length(intersect(unique(rv$Path_data2[,1]),
                unique(rv$DEA_data2[,1])))==length(unique(rv$DEA_data2[,1]))) {
                return(paste0('The uploaded differential expression
                    file contains: ',
                length(intersect(unique(rv$Path_data2[,1]),rv$DEA_data2[,1])),
                ' targets in common with Path-to-Target data. All targets in the
                differential expression file are contemplated in the
                Path-to-Target file!'))
            }
                else{
                    if ((length(intersect(unique(rv$Path_data2[,1]),
                        unique(rv$DEA_data2[,1])))-
                        length(unique(rv$DEA_data2[,1])))<0) {
                        return(paste0('The uploaded differential expression
                            file contains: ',
                        length(intersect(unique(rv$Path_data2[,1]),
                            rv$DEA_data2[,1])),
                        ' targets in common with Path-to-Target data.
                        There are ',
                        abs(length(intersect(unique(rv$Path_data2[,1]),
                            unique(rv$DEA_data2[,1])))-
                            length(unique(rv$DEA_data2[,1]))),
                        ' targets in the differential expression file
                        not contemplated in your Path-to-Target file!'))
                    }
                        else{
                            if ((length(intersect(unique(rv$Path_data2[,1]),
                                unique(rv$DEA_data2[,1])))-
                                length(unique(rv$DEA_data2[,1])))>0) {
                                return(paste0('The uploaded differential
                                    expression file contains: ',
                                length(intersect(unique(rv$Path_data2[,1]),
                                    rv$DEA_data2[,1])),
                                ' targets in common with Path-to-Target data.
                                There are ',
                                abs(length(intersect(unique(rv$Path_data2[,1]),
                                    unique(rv$DEA_data2[,1])))-
                                    length(unique(rv$DEA_data2[,1]))),
                                ' targets in the Path-to-Target file not
                                contemplated in your differential
                                expression file!'))
                            }
                        }
                }
            }
            else{
               NULL
            }
    })


    output$DEA_text6<-renderText({
        if(((((!is.null(rv$GFAG_data)) & !(is.null(rv$NEX_data1))) &
            !is.null(rv$Path_data1))) & !is.null(rv$DEA_data1) &
            (rv$genekey==1 & rv$logkey==1 & rv$pvalkey==1 )){
            return('Diferential expression file sucessfully uploaded;
            columns sucessfully selected!')
        }
            else{
                NULL
            }
    })


    output$DEAui<-renderUI({
        if(((((is.null(rv$GFAG_data)) & (is.null(rv$NEX_data1))) &
            is.null(rv$Path_data1))) & is.null(rv$DEA_data1)){
                shinyjs::disabled(
                    wellPanel(style = "background-color: #c6c6c6;border:
                              1px #c6c6c6 solid;",
                        textOutput('DEA_text1'),
                        textOutput('DEA_text2'),
                        textOutput('DEA_text3'),
                        textOutput('DEA_text4'),
                        #textOutput('DEA_text5'),
                        #textOutput('DEA_text6'),
                        tags$style(type="text/css",
                        ".shiny-output-error { visibility: hidden; }",
                        ".shiny-output-error:before { visibility: hidden; }"),
                        tags$head(tags$style("#DEA_text1{color: blue;
                            font-size: 15px; font-style: italic;}"))))
        }
            else{
                if(((((!is.null(rv$GFAG_data)) & (is.null(rv$NEX_data1))) &
                    is.null(rv$Path_data1))) & is.null(rv$DEA_data1)){
                    shinyjs::disabled(
                        wellPanel(
                            style = "background-color: #c6c6c6;border:
                                1px #c6c6c6 solid;",
                            textOutput('DEA_text1'),
                            textOutput('DEA_text2'),
                            textOutput('DEA_text3'),
                            textOutput('DEA_text4'),
                            #textOutput('DEA_text5'),
                            #textOutput('DEA_text6'),
                            tags$style(type="text/css",
                            ".shiny-output-error { visibility: hidden; }",
                            ".shiny-output-error:before { visibility: hidden;
                            }"),
                            tags$head(tags$style("#DEA_text2{color: blue;
                                font-size: 15px; font-style: italic;}"))))
                }
                    else{
                        if(((((!is.null(rv$GFAG_data)) &
                            !(is.null(rv$NEX_data1))) &
                            is.null(rv$Path_data1))) & is.null(rv$DEA_data1)){
                            shinyjs::disabled(
                                wellPanel(
                                    style = "background-color: #c6c6c6;border:
                                        1px #c6c6c6 solid;",
                                    textOutput('DEA_text1'),
                                    textOutput('DEA_text2'),
                                    textOutput('DEA_text3'),
                                    textOutput('DEA_text4'),
                                    #textOutput('DEA_text5'),
                                    #textOutput('DEA_text6'),
                                    tags$style(type="text/css",
                                        ".shiny-output-error { visibility:
                                        hidden; }",
                                        ".shiny-output-error:before
                                        { visibility: hidden; }"),
                                        tags$head(
                                            tags$style("#DEA_text3{color: blue;
                                                font-size: 15px;
                                                font-style: italic;}"))))
                        }
                            else{
                                if(((((!is.null(rv$GFAG_data)) &
                                    !(is.null(rv$NEX_data1))) &
                                    !is.null(rv$Path_data1))) &
                                    is.null(rv$DEA_data1)){
                                    tagList(
                                        wellPanel(style="background-color:
                                            #ffffff;border: 1px #86e6fc solid;
                                            ",
                                        fileInput(label = 'Upload Differential
                                            Expression Analysis Data for a
                                            choosen Case/Control:',
                                            inputId = 'DEAdata'),
                                        textOutput('DEA_text1'),
                                        textOutput('DEA_text2'),
                                        textOutput('DEA_text3'),
                                        textOutput('DEA_text4'),
                                        # textOutput('DEA_text5'),
                                        #  textOutput('DEA_text6'),
                                        tags$style(type="text/css",
                                            ".shiny-output-error {
                                            visibility: hidden; }",
                                        ".shiny-output-error:before
                                        { visibility: hidden; }"),
                                        tags$head(tags$style(
                                            "#DEA_text4{color: blue;
                                            font-size: 15px;
                                            font-style: italic;}")),
                                        tags$hr(),
                                        helpText("Select the read.table
                                            parameters below:"),
                                        checkboxInput(inputId = 'headerDEA',
                                            label = 'Header',
                                            value = TRUE),
                                        checkboxInput(
                                            inputId = "stringAsFactorsDEA",
                                            "stringAsFactors",
                                            FALSE),
                                        helpText('Select the separator type:'),
                                        radioButtons(inputId = 'sepDEA',
                                            label = 'Separator',
                                            choices = c(Comma=',',
                                                Semicolon=';',
                                                Tab='\t',
                                                Space=''),
                                            selected = '\t'),
                                        tags$hr(),
                                        p('Column Selection:',
                                            style='font-size: 22px;
                                            font-style: bold;'),
                                        p('Upload a diferential expression
                                            file first!',
                                            style='font-style: bold;')))
                                    }
                                    else{
                                        if(((((!is.null(rv$GFAG_data)) &
                                            !(is.null(rv$NEX_data1))) &
                                            !is.null(rv$Path_data1))) &
                                            !(is.null(rv$DEA_data1)) &
                                            (rv$genekey==0 |
                                            rv$logkey==0 | rv$pvalkey==0 )){
                                            tagList(
                                                wellPanel(
                                                style="background-color:
                                                #c6c6c6;border:
                                                1px #86e6fc solid;",
                                                fileInput(
                                                label = 'Upload Differential
                                                Expression Analysis Data for a
                                                choosen Case/Control:',
                                                inputId = 'DEAdata'),
                                                textOutput('DEA_text1'),
                                                textOutput('DEA_text2'),
                                                textOutput('DEA_text3'),
                                                textOutput('DEA_text4'),
                                                textOutput('DEA_text5'),
                                                textOutput('DEA_text6'),
                                                tags$style(type="text/css",
                                                    ".shiny-output-error
                                                    {visibility: hidden; }",
                                                    ".shiny-output-error:before
                                                    {visibility: hidden; }"),
                                                tags$head(tags$style(
                                                    "#DEA_text5{color: green;
                                                    font-size: 15px;
                                                    font-style: italic;}")),
                                                tags$hr(),
                                                helpText("Select the read.table
                                                    parameters below:"),
                                                checkboxInput(inputId =
                                                'headerDEA',
                                                label = 'Header',
                                                value = TRUE),
                                                checkboxInput(inputId =
                                                    "stringAsFactorsDEA",
                                                    "stringAsFactors", FALSE),
                                                helpText('Select the
                                                    separator type:'),
                                                radioButtons(inputId =
                                                    'sepDEA',
                                                    label = 'Separator',
                                                    choices = c(Comma=',',
                                                        Semicolon=';',
                                                        Tab='\t',
                                                        Space=''),
                                                    selected = '\t'),
                                                tags$hr(),
                                                p('Column Selection:',
                                                    style='font-size: 22px;
                                                    font-style: bold;'),
                                                hr(),
                                                uiOutput('selectGene2'),
                                                actionButton(
                                                    inputId = 'upGene2',
                                                    label = 'Update selected
                                                    gene column'),
                                                hr(),
                                                uiOutput('selectlog'),
                                                actionButton(inputId = 'uplog',
                                                    label = 'Update Selected
                                                    logFC column'),
                                                hr(),
                                                uiOutput('selectpval'),
                                                actionButton(
                                                    inputId = 'uppval',
                                                    label = 'Update Selected
                                                    P-value column')))
                                            }
                                                else{
                                                    if(((((
                                                        !is.null(
                                                            rv$GFAG_data)) &
                                                        !(is.null(
                                                            rv$NEX_data1))) &
                                                        !is.null(
                                                            rv$Path_data1))) &
                                                        !is.null(rv$DEA_data1)
                                                        & ((rv$genekey==1
                                                        & rv$logkey==1 &
                                                        rv$pvalkey==1 ))){
                                                        shinyjs::disabled(
                                                        wellPanel(
                                                        style =
                                                        "background-color:
                                                        #c6c6c6;
                                                        border: 1px
                                                        #c6c6c6 solid;",
                                                        textOutput(
                                                            'DEA_text1'),
                                                        textOutput(
                                                            'DEA_text2'),
                                                        textOutput(
                                                            'DEA_text3'),
                                                        textOutput(
                                                            'DEA_text4'),
                                                        textOutput(
                                                            'DEA_text5'),
                                                        textOutput(
                                                            'DEA_text6'),
                                                        tags$style(
                                                            type="text/css",
                                                        ".shiny-output-error
                                                        {visibility: hidden;
                                                        }",
                                                        ".shiny-output-error:
                                                        before { visibility:
                                                        hidden; }"),
                                                        tags$head(tags$style(
                                                        "#DEA_text6{color:
                                                        green;
                                                        font-size: 15px;
                                                        font-style:
                                                        italic;}"))))
                                                    }
                                                }
                                    }
                            }
                    }
            }
    })


    output$DEA_table<-DT::renderDataTable({
        if(!(is.null(rv$DEA_data2))){
            rv$DEA_data2
        }
            else{
                return()
            }
    })


    output$selectGene2<-renderUI({
        if (is.null(rv$DEA_data2)) {
            shinyjs::disabled(selectInput(label = 'Select a Target column:',
                inputId = 'Genesel2',choices = c()))
        }
            else{
                if(!(is.null(rv$DEA_data2))){
                    shinyjs::enable('upGene2')
                    selectInput(label = 'Select a Target column:',
                        inputId = 'Genesel2',
                        choices = colnames(rv$DEA_data2))
                }
                    
            }
    })
    

    observeEvent(input$upGene2,{
        rv$Gene2_index<-grep(paste0('^',input$Genesel2,'$'),
                             colnames(rv$DEA_data2))
        rv$colselsum <- rv$colselsum + 1
        rv$genekey<-1
        rv$gkey2<-'ok'
    })


    output$selectlog<-renderUI({
        if (is.null(rv$DEA_data2)) {
        shinyjs::disabled(selectInput(label = 'Select the LogFC column:',
            inputId = 'logsel',choices = c()))
        }
            else{
                if(!(is.null(rv$DEA_data2))){
                    shinyjs::enable('uplog')
                    selectInput(label = 'Select a LogFC column:',
                        inputId = 'logsel',
                        choices = colnames(rv$DEA_data2))
                }
            }
    })


    output$logtest1<-renderText({
        if (is.null(rv$logkey)) {
        'hello'
        }
    })


    observeEvent(input$uplog,{
        rv$logsel_index<-grep(paste0('^',input$logsel,'$'),
            colnames(rv$DEA_data2))
        rv$colselsum <- rv$colselsum + 1
        rv$logkey<-1
    })


    output$selectpval<-renderUI({
        if (is.null(rv$DEA_data2)) {
        shinyjs::disabled(selectInput(label = 'Select the P-value column:',
            inputId = 'pvalsel',choices = c()))
        }
            else{
                if(!(is.null(rv$DEA_data2))){
                    shinyjs::enable('uppval')
                    selectInput(label = 'Select the P-value column:',
                    inputId = 'pvalsel',
                    choices = colnames(rv$DEA_data2))
                }
            }
    })


    observeEvent(input$uppval,{
        rv$pvalsel_index<-grep(paste0('^',input$pvalsel,'$'),
            colnames(rv$DEA_data2))
        rv$colselsum <- rv$colselsum + 1
        rv$pvalkey<-1
    })

    ###########################################################################



    observeEvent(input$pathgene,{
        rv$path_index<-NULL
        rv$gene_index<-NULL
        file3 <- input$pathgene
        if(is.null(file3)){
            rv$path_data1<-NULL
        }
            else{
                e<-read.delim(file=file3$datapath,
                sep='\t',
                header = TRUE,
                stringsAsFactors = TRUE)
                f<-e
                if (is.null(e)) {
                    rv$Path_data1<-e
                }
                    else{
                        rv$viasum<-sum(unique(e[,2]) %in%
                        unique(rv$GFAG_data2$ID))
                        rv$elemsum<-sum(unique(e[,1])
                        %in% unique(rv$NEX_data2$gene))
                        rv$Path_data1<-e
                        rv$Path_data2<-NULL
                        if (!(is.null(e)) & (length(colnames(e))==2)) {
                            if (rv$viasum>0 & rv$elemsum>0) {
                                rv$Path_data1<-e
                                rv$Path_data2<-e
                            }
                        }
                    }
            }
    })


    output$Path_text1<-renderText({
        if((((is.null(rv$GFAG_data2)) & (is.null(rv$NEX_data2))) &
            is.null(rv$Path_data1))){
            return('Warning: You need to load a GFAG output file first!')
        }
            else{
                NULL
            }
    })


    output$Path_text11<-renderText({
        if((((is.null(rv$GFAG_data2)) & (is.null(rv$NEX_data2))) &
            is.null(rv$Path_data1))){
            return('The uploaded Path-to-Target file contains:
                No Path-to-Target file was uploaded!')
        }
            else{
                NULL
            }
    })


    output$Path_text2<-renderText({
        if((!(is.null(rv$GFAG_data2)) & (is.null(rv$NEX_data2))) &
            is.null(rv$Path_data1)){
            return('Warning: You need to load a expression file first!')
        }
            else{
                NULL
            }
    })


    output$Path_text21<-renderText({
        if((!(is.null(rv$GFAG_data2)) & (is.null(rv$NEX_data2))) &
            is.null(rv$Path_data1)){
            return('The uploaded Path-to-Target file contains:
                No Path-to-Target file was uploaded!')
        }
            else{
                NULL
            }
    })


    output$Path_text3<-renderText({
        if((!(is.null(rv$GFAG_data2)) & !(is.null(rv$NEX_data2))) &
            is.null(rv$Path_data1)) {
            return('Warning: No Path-to-Target data was uploaded!')
        }
            else{
                NULL
            }
    })


    output$Path_text31<-renderText({
        if((!(is.null(rv$GFAG_data2)) & !(is.null(rv$NEX_data2))) &
            is.null(rv$Path_data1)) {
            return('The uploaded Path-to-Target file contains:
                No Path-to-Target file was uploaded!')
        }
            else{
                NULL
            }
    })


    output$Path_text4<-renderText({
        if((!(is.null(rv$GFAG_data2)) & !(is.null(rv$NEX_data2))) &
            !(is.null(rv$Path_data1)) & (length(colnames(rv$Path_data1))==2)
            & rv$viasum>0 & rv$elemsum>0) {
            return('Path-to-Target data sucessfully uploaded!')
        }
            else{
                NULL
            }

    })


    output$Path_text40<-renderText({
        req(input$gfagdata)
        req(input$NEXdata)
        req(input$pathgene)
        if((!(is.null(rv$GFAG_data2)) & !(is.null(rv$NEX_data2))) &
            !(is.null(rv$Path_data1)) &
            (length(colnames(rv$Path_data1))==2) &
            rv$viasum>0 & rv$elemsum>0) {
            return(paste0('The uploaded Path-to-Target file contains:
                A total of ',
                length((rv$Path_data2[,1])),' entries.'))
            }
                else{
                    NULL
                }
    })


    output$Path_text41<-renderText({
        req(input$gfagdata)
        req(input$NEXdata)
        req(input$pathgene)
        if((!(is.null(rv$GFAG_data2)) & !(is.null(rv$NEX_data2))) &
            !(is.null(rv$Path_data1)) &
            (length(colnames(rv$Path_data1))==2) &
            rv$viasum>0 & rv$elemsum>0) {
            return(paste0('The uploaded Path-to-Target file contains: ',
                length(unique(rv$Path_data2[,1])),' unique targets'))
        }
            else{
                NULL
            }
    })


    output$Path_text42<-renderText({
        req(input$gfagdata)
        req(input$NEXdata)
        req(input$pathgene)
        if((!(is.null(rv$GFAG_data2)) & !(is.null(rv$NEX_data2))) &
            !(is.null(rv$Path_data1)) &
            (length(colnames(rv$Path_data1))==2) &
            rv$viasum>0 & rv$elemsum>0) {
            return(paste0('The uploaded Path-to-Target file contains: ',
                length(unique(rv$Path_data2[,2])),' unique paths'))
        }
            else{
                NULL
            }
    })


    output$Path_text43<-renderText({
        req(input$gfagdata)
        req(input$NEXdata)
        req(input$pathgene)
        if((!(is.null(rv$GFAG_data2)) & !(is.null(rv$NEX_data2))) &
            !(is.null(rv$Path_data1)) & (length(colnames(rv$Path_data1))==2) &
            rv$viasum>0 & rv$elemsum>0) {
            if(length(intersect(unique(rv$Path_data2[,2]),
                rv$GFAG_data2[,1]))==length(rv$GFAG_data2[,1])){
                return(paste0('The uploaded Path-to-Target file contains: ',
                    length(intersect(unique(rv$Path_data2[,2]),
                        rv$GFAG_data2[,1])),
                    ' paths in common with GFAG data. All of your GFAGs
                    are contemplated in the Path-to-Target data.'))
            }
                else{
                    if ((length(intersect(unique(rv$Path_data2[,2]),
                        rv$GFAG_data2[,1]))-
                        length(unique(rv$Path_data2[,2])))<0) {
                        return(paste0('The uploaded Path-to-Target
                            file contains: ',
                            length(intersect(unique(rv$Path_data2[,2]),
                                rv$GFAG_data2[,1])),
                            ' paths in common with GFAG data. There are ',
                            abs(length(intersect(unique(rv$Path_data2[,2]),
                                rv$GFAG_data2[,1]))-
                                length(unique(rv$Path_data2[,2]))),
                            ' GFAGs in your Path-to-Target file not
                            contemplated in your output file!'))
                    }
                        else{
                            if ((length(intersect(unique(rv$Path_data2[,2]),
                                rv$GFAG_data2[,1]))-
                                length(unique(rv$Path_data2[,2])))>0) {
                                return(paste0('The uploaded Path-to-Target
                                    file contains: ',
                                    length(intersect(unique(rv$Path_data2[,2]),
                                        rv$GFAG_data2[,1])),
                                    ' paths in common with GFAG data.
                                    There are ',
                                    abs(length(intersect(unique(
                                        rv$Path_data2[,2]),
                                        rv$GFAG_data2[,1]))-
                                        length(unique(rv$Path_data2[,2]))),
                                    ' GFAGs in your output file not
                                    contemplated in your
                                    Path-to-Target file!'))

                            }
                        }
                }
        }
            else{
                NULL
            }
    })


    output$Path_text44<-renderText({
        req(input$gfagdata)
        req(input$NEXdata)
        req(input$pathgene)
        if((!(is.null(rv$GFAG_data2)) & !(is.null(rv$NEX_data2))) &
            !(is.null(rv$Path_data1)) & (length(colnames(rv$Path_data1))==2)
            & rv$viasum>0 & rv$elemsum>0) {
                if (length(intersect(unique(rv$Path_data2[,1]),
                    rv$NEX_data2[,1]))==length(unique(rv$Path_data2[,1]))) {
                    return(paste0('The uploaded Path-to-Target file
                        contains: ',
                    length(intersect(unique(rv$Path_data2[,1]),
                    rv$NEX_data2[,1])),
                    ' targets in common with expression data.'))
                }
                    else{
                        if ((length(intersect(unique(rv$Path_data2[,1]),
                            rv$NEX_data2[,1]))-
                            length(unique(rv$Path_data2[,1])))<0) {
                            return(paste0('The uploaded Path-to-Target file
                                contains: ',
                            length(intersect(unique(rv$Path_data2[,1]),
                                rv$NEX_data2[,1])),
                            ' targets in common with expression data.
                            There are ',
                            abs(length(intersect(unique(rv$Path_data2[,1]),
                                rv$NEX_data2[,1]))-
                                length(unique(rv$Path_data2[,1]))),
                            ' targets in your Path-to-Target file not
                            contemplated in your expression file!'))
                        }
                            else{
                                if ((length(intersect(unique(
                                    rv$Path_data2[,1]),rv$NEX_data2[,1]))-
                                    length(unique(rv$Path_data2[,1])))>0) {
                                    return(paste0('The uploaded Path-to-Target
                                        file contains: ',
                                    length(intersect(unique(rv$Path_data2[,1]),
                                        rv$NEX_data2[,1])),
                                    ' targets in common with expression data.
                                    There are ',
                                    abs(length(intersect(unique(
                                        rv$Path_data2[,1]),rv$NEX_data2[,1]))-
                                    length(unique(rv$Path_data2[,1]))),
                                    ' targets in your expression file not
                                    contemplated in your Path-to-Target file!'))
                                }
                            }
                    }
        }
            else{
                NULL
            }
    })


    output$Path_text5<-renderText({
        if((!(is.null(rv$GFAG_data2)) & !(is.null(rv$NEX_data2))) &
            !(is.null(rv$Path_data1)) &
            ((length(colnames(rv$Path_data1)))!=2) |
            (rv$viasum==0 | rv$elemsum==0)) {
            return('Error: Incorrect Path-to-Target file format!')
        }
            else{
                NULL
            }
    })


    output$Path_text51<-renderText({
        req(input$gfagdata)
        req(input$NEXdata)
        req(input$pathgene)
        if((!(is.null(rv$GFAG_data2)) & !(is.null(rv$NEX_data2))) &
            !(is.null(rv$Path_data1)) &
            ((length(colnames(rv$Path_data1)))!=2) |
            (rv$viasum==0 | rv$elemsum==0)) {
            return('The uploaded Path-to-Target file contains:
                No Path-to-Target file was uploaded!')
        }
            else{
                NULL
           }
    })


    output$Pathui1<-renderUI ( {
        if ((((is.null(rv$GFAG_data2)) & (is.null(rv$NEX_data2))) &
            is.null(rv$Path_data1))) {
            shinyjs::disabled(wellPanel(style = "background-color:
                #c6c6c6;border: 1px #c6c6c6 solid;",
            textOutput('Path_text1'),
            tags$head(tags$style("#Path_text1{color: blue; font-size: 15px;
                font-style: italic;}"))))
        }
            else{
                if((!(is.null(rv$GFAG_data2)) &
                    (is.null(rv$NEX_data2))) &
                    is.null(rv$Path_data1)) {
                    shinyjs::disabled(
                        wellPanel(style = "background-color: #c6c6c6;border:
                            1px #c6c6c6 solid;",
                        textOutput('Path_text2'),
                        tags$head(tags$style("#Path_text2{color: blue;
                            font-size: 15px; font-style: italic;}"))))
                }
                    else{
                        if((!(is.null(rv$GFAG_data2)) &
                            !(is.null(rv$NEX_data2))) &
                            is.null(rv$Path_data1)) {
                            wellPanel(
                                style="background-color: #ffffff;border:
                                    1px #86e6fc solid;",
                                fileInput(inputId = 'pathgene',
                                    label = 'Select a file with a pathway
                                    (Gene ontology, KEGG, etc) to gene
                                    relationship:'),
                                textOutput('Path_text3'),
                                tags$head(tags$style("#Path_text3{color: blue;
                                    font-size: 15px; font-style: italic;}")))
                        }
                            else{
                                if((!(is.null(rv$GFAG_data2)) &
                                    !(is.null(rv$NEX_data2))) &
                                    !(is.null(rv$Path_data1)) &
                                    ((length(colnames(rv$Path_data1))!=2)
                                    | (rv$viasum==0 | rv$elemsum==0))) {
                                    wellPanel(style="background-color: #ffffff;
                                        border: 1px #86e6fc solid;",
                                    fileInput(inputId = 'pathgene',
                                        label = 'Select a file with a pathway
                                        (Gene ontology, KEGG, etc) to gene
                                        relationship:'),
                                    textOutput('Path_text5'),
                                    tags$head(tags$style(
                                        "#Path_text5{color: red; font-size:
                                        15px; font-style: italic;}")))
                            }
                                else{
                                    if((!(is.null(rv$GFAG_data2)) &
                                        !(is.null(rv$NEX_data2))) &
                                        !(is.null(rv$Path_data1)) &
                                        (length(colnames(rv$Path_data1))==2) &
                                        rv$viasum>0 & rv$elemsum>0) {
                                        shinyjs::disabled(
                                            wellPanel(
                                            style = "background-color:
                                            #c6c6c6;border:
                                            1px #c6c6c6 solid;",
                                        textOutput('Path_text4'),
                                        tags$head(tags$style("#Path_text4{
                                        color: green; font-size: 15px;
                                        font-style: italic;}"))))
                                    }
                                }
                            }
                    }

            }
    })


    output$Path_table<-DT::renderDataTable({
        if(!(is.null(rv$Path_data2))){
            rv$Path_data2
        }
            else{
                return()
            }
    })


    output$selectPath<-renderUI({
        if (is.null(rv$Path_data1)) {
            shinyjs::disabled(selectInput(label = 'Select the Pathway column:',
                inputId = 'Pathsel',choices = c()))
        }
            else{
                if(!(is.null(rv$Path_data1))){
                    shinyjs::enable('upPath')
                    selectInput(label = 'Select a Pathway column:',
                        inputId = 'Pathsel',choices = colnames(rv$Path_data1))
                }
            }
    })


    observeEvent(input$upPath,{
        rv$path_index<-grep(paste0('^',input$Pathsel,'$'),
            colnames(rv$Path_data1))
    })


    output$selectGene<-renderUI({
        if (is.null(rv$Path_data1)) {
            shinyjs::disabled(selectInput(label = 'Select the Gene column:',
                inputId = 'Genesel',choices = c()))
        }
            else{
                if(!(is.null(rv$Path_data1))){
                    shinyjs::enable('upGene')
                    selectInput(label = 'Select the Gene column:',
                    inputId = 'Genesel',
                    choices = colnames(rv$Path_data1))
                }
            }
    })


    observeEvent(input$upGene,{
        rv$gene_index<-grep(paste0('^',input$Genesel,'$'),
            colnames(rv$Path_data1))
    })


  #############################################################################



    observeEvent(input$gfagdata,{
        file3 <- input$gfagdata
        if(is.null(file3)){
            rv$GFAG_data<-NULL
            rv$genekey<-0
            rv$logkey<-0
            rv$pvalkey<-0
        }
            else{
                rv$GFAG_data0<-read.delim(file=file3$datapath,
                    sep='\t',
                    header = TRUE,
                    stringsAsFactors = TRUE)
                if ((!(is.null(rv$GFAG_data0)) &
                    grepl('H_',colnames(rv$GFAG_data0)[5])) |
                    (!(is.null(rv$GFAG_data0)) &
                    (length(colnames(rv$GFAG_data0))>=16 &
                    length(colnames(rv$GFAG_data0))<=22))) {
                    rv$GFAG_data<-read.delim(file=file3$datapath,
                        sep='\t',
                        header = TRUE,
                        stringsAsFactors = TRUE)
                    rv$GFAG_data2<-read.delim(file=file3$datapath,
                        sep='\t',
                        header = TRUE,
                        stringsAsFactors = TRUE)
                    rv$Controlid<-str_split(
                        colnames(rv$GFAG_data0)[5],'_')[[1]][2]
                    rv$Caseid<-str_split(
                        colnames(rv$GFAG_data0)[6],'_')[[1]][2]
                }
                    else{
                        if((!(is.null(rv$GFAG_data0)) &
                            !(grepl('H_',colnames(rv$GFAG_data0)[5])))
                            | !(!(is.null(rv$GFAG_data0)) &
                            (length(colnames(rv$GFAG_data0))>=16 &
                            length(colnames(rv$GFAG_data0))<=22))) {
                            rv$GFAG_data<-read.delim(file=file3$datapath,
                            sep='\t', header = TRUE, stringsAsFactors = TRUE)
                            rv$GFAG_data2<-NULL
                            rv$Controlid<-NULL
                            rv$Caseid<-NULL
                        }
                    }
            }
    })


    output$GFAG_text1<-renderText({
        if((is.null(rv$GFAG_data))){
            return('Warning: No GFAG output file was uploaded!')
        }
            else{
                NULL
            }
    })


    output$GFAG_text11<-renderText({
        if((is.null(rv$GFAG_data))){
            return(paste0('The uploaded GFAG file contains: ',
                'No GFAG output file was uploaded!'))
        }
            else{
                NULL
            }
    })


    output$GFAG_text2<-renderText({
        req(input$gfagdata)
        req(rv$GFAG_data)
        if ((!(is.null(rv$GFAG_data)) &
            grepl('H_',colnames(rv$GFAG_data)[5])) |
            (!(is.null(rv$GFAG_data)) &
            (length(colnames(rv$GFAG_data))>=16 &
            length(colnames(rv$GFAG_data))<=22))) {
            return('GFAG output file sucessfully uploaded!')
        }
            else{
                NULL
            }
    })


    output$GFAG_text21<-renderText({
        req(input$gfagdata)
        req(rv$GFAG_data)
        if ((!(is.null(rv$GFAG_data)) &
            grepl('H_',colnames(rv$GFAG_data)[5])) |
            (!(is.null(rv$GFAG_data)) &
            (length(colnames(rv$GFAG_data))>=16 &
            length(colnames(rv$GFAG_data))<=22))) {
            return(paste0('The uploaded GFAG file contains: A total of ',
                length(rownames(rv$GFAG_data2)),' entries.'))
        }
            else{
                NULL
            }
    })


    output$GFAG_text22<-renderText({
        req(input$gfagdata)
        req(rv$GFAG_data)
        if ((!(is.null(rv$GFAG_data)) &
            grepl('H_',colnames(rv$GFAG_data)[5])) |
            (!(is.null(rv$GFAG_data)) &
            (length(colnames(rv$GFAG_data))>=16 &
            length(colnames(rv$GFAG_data))<=22))) {
            return(paste0('The uploaded GFAG file contains: ',
                length(unique(rv$GFAG_data2[,1])),' unique GFAGs.'))
        }
            else{
                NULL
            }
    })


    output$GFAG_text3<-renderText({
        if ((!(is.null(rv$GFAG_data)) &
            !(grepl('H_',colnames(rv$GFAG_data)[5]))) |
            !(!(is.null(rv$GFAG_data)) &
            (length(colnames(rv$GFAG_data))>=16 &
            length(colnames(rv$GFAG_data))<=22))) {
            return('Error: Failed to upload GFAG output file!
                The selected file is not in the correct format,
                possible issues may include incorrect number of
                columns or incorrect column naming.
                Please remember that the GFAG output file should
                not be altered!')
        }
    })


    output$GFAG_text31<-renderText({
        req(input$gfagdata)
        req(input$NEXdata)
        req(input$pathgene)
        req(input$DEA_data)
        if ((!(is.null(rv$GFAG_data)) &
            !(grepl('H_',colnames(rv$GFAG_data)[5]))) |
            !(!(is.null(rv$GFAG_data)) &
            (length(colnames(rv$GFAG_data))>=16 &
            length(colnames(rv$GFAG_data))<=22))) {
            return(paste0('The uploaded GFAG file contains: ',
                'No GFAG output file was uploaded!'))
        }
    })


    output$GFAGui<-renderUI({
        if (is.null(rv$GFAG_data)) {
            tagList(wellPanel(
                style="background-color: #ffffff; border: 1px
                    #86e6fc solid;",
                fileInput(label = 'Upload a GFAG output file for
                    choosen Case/Control:',
                inputId = 'gfagdata'),
                textOutput('GFAG_text1'),
                tags$head(tags$style("#GFAG_text1{color: blue; font-size: 15px;
                    font-style: italic; }"))))
        }
            else{
                if ((!(is.null(rv$GFAG_data)) &
                    grepl('H_',colnames(rv$GFAG_data)[5])) |
                    (!(is.null(rv$GFAG_data)) &
                    (length(colnames(rv$GFAG_data))>=16 &
                    length(colnames(rv$GFAG_data))<=22))) {
                    tagList((wellPanel(
                        style = "background-color:
                        #c6c6c6;border: 1px #c6c6c6 solid;",
                        textOutput('GFAG_text2'),
                        tags$head(tags$style("#GFAG_text2{color: green;
                        font-size: 15px; font-style: italic;}")))))
                }
                    else{
                        if((!(is.null(rv$GFAG_data)) &
                            !(grepl('H_',colnames(rv$GFAG_data)[5]))) |
                            !(!(is.null(rv$GFAG_data)) &
                            (length(colnames(rv$GFAG_data))>=16 &
                            length(colnames(rv$GFAG_data))<=22))) {
                            tagList(wellPanel(style="background-color:
                                #ffffff;border: 1px #ffffff solid;",
                                fileInput(label = 'Upload a GFAG output
                                    file for choosen Case/Control:',
                                    inputId = 'gfagdata'),
                                textOutput('GFAG_text3'),
                                tags$head(tags$style("#GFAG_text3{color: red;
                                     font-size: 15px;
                                     font-style: italic;}"))))
                        }
                    }
            }
    })


    output$GFAG_table<-DT::renderDataTable({
        if(!(is.null(rv$GFAG_data2))){
            rv$GFAG_data2
        }
            else{
                return()
            }
    })


    observeEvent(input$Spcsel,{
       rv$Species<-input$Species
    })


    observeEvent(input$Reset,{
        shinyjs::disable('upcase')
        shinyjs::disable('upcontrol')
        shinyjs::disable('upgfag')
        shinyjs::disable('uplog')
        shinyjs::disable('uppval')
        shinyjs::disable('ViewH')
        shinyjs::disable('ViewN')
        shinyjs::disable('ViewNH')
        shinyjs::disable('ViewR')
        shinyjs::disable('upfiltered')
        shinyjs::disable('plot_start')
        rv$Case_index<-NULL
        rv$Control_index<-NULL
        rv$logsel_index<-NULL
        rv$pvalsel_index<-NULL
        rv$NEX_data1<-NULL
        rv$NEX_data2<-NULL
        rv$DEA_data1<-NULL
        rv$DEA_data2<-NULL
        rv$GFAG_data<-NULL
        rv$GFAG_data2<-NULL
        rv$Path_data1<-NULL
        rv$Path_data2<-NULL
        rv$CACO_cols<-NULL
        rv$QVAL_cols<-NULL
        rv$LOGFC_cols<-NULL
        rv$gg<-NULL
        rv$filtered_dt<-NULL
        rv$filtered_gfag<-NULL
        rv$genekey<-0
        rv$logkey<-0
        rv$pvalkey<-0
        rv$startfilter<-0
        rv$filtereddt<-0
        rv$ftable<-NULL
        shinyjs::reset('NEXdata')
        shinyjs::reset('pathgene')
        shinyjs::reset('DEAdata')
        shinyjs::reset('gfagdata')
    })


    output$GFAG_filter<-renderUI({
        if (!(is.null(rv$Case_index)) & !(is.null(rv$Control_index)) &
            !(is.null(rv$logsel_index)) & !(is.null(rv$pvalsel_index)) &
            !(is.null(rv$GFAG_data))) {
            actionButton(inputId = 'Gofilter',
                label = 'Proced to GFAG selection')
        }
            else{
                shinyjs::disabled(actionButton(inputId = 'Gofilter',
                    label = 'Proced to GFAG selection'))
            }
    })



  ############### Second tab ####################


    output$NN1<-renderUI({
        if (!(is.null(rv$Gene2_index)) & !(is.null(rv$logsel_index)) &
            !(is.null(rv$pvalsel_index)) & !(is.null(rv$GFAG_data2))) {
            numericInput(inputId = 'qvalN',
                label = 'Select a q-value cutoff for activity:',
                value = 0.05)
        }
            else{
                shinyjs::disabled(numericInput(inputId = 'qvalN',
                    label = 'Select a q-value cutoff for activity:',
                    value = 0.05))
            }
    })


    output$NN2<-renderUI({
        if (!(is.null(rv$Gene2_index)) & !(is.null(rv$logsel_index)) &
            !(is.null(rv$pvalsel_index)) & !(is.null(rv$GFAG_data2))) {
            numericInput(inputId = 'qvalH',
                label = 'Select a q-value cutoff for diversity:',
                value = 0.05)
        }
            else{
                shinyjs::disabled(numericInput(inputId = 'qvalH',
                    label = 'Select a q-value cutoff for diversity:',
                    value = 0.05))
            }
    })


    output$BT1<-renderUI({
        if (!(is.null(rv$Gene2_index)) & !(is.null(rv$logsel_index)) &
            !(is.null(rv$pvalsel_index)) & !(is.null(rv$GFAG_data2))) {
            actionButton(inputId = 'ViewN',label = 'Filter Activity')
        }
            else{
                shinyjs::disabled(actionButton(inputId = 'ViewN',
                    label = 'Filter Activity'))
            }
    })


    output$BT2<-renderUI({
        if (!(is.null(rv$Gene2_index)) & !(is.null(rv$logsel_index)) &
            !(is.null(rv$pvalsel_index)) & !(is.null(rv$GFAG_data2))) {
            actionButton(inputId = 'ViewH',label = 'Filter Diversity')
        }
            else{
                shinyjs::disabled(actionButton(inputId = 'ViewH',
                    label = 'Filter Diversity'))
            }
    })


    output$BT3<-renderUI({
        if (!(is.null(rv$Gene2_index)) & !(is.null(rv$logsel_index)) &
            !(is.null(rv$pvalsel_index)) & !(is.null(rv$GFAG_data2))) {
            actionButton(inputId = 'ViewNH',
                label = 'Filter Activity & Diversity')
        }
            else{
                shinyjs::disabled(actionButton(inputId = 'ViewNH',
                    label = 'Filter Activity & Diversity'))
            }
    })


    output$BT4<-renderUI({
        if (!(is.null(rv$Gene2_index)) & !(is.null(rv$logsel_index)) &
            !(is.null(rv$pvalsel_index)) & !(is.null(rv$GFAG_data2))) {
            actionButton(inputId = 'ViewR',label = 'No Filter')
        }
            else{
                shinyjs::disabled(actionButton(inputId = 'ViewR',
                    label = 'No Filter'))
            }
    })


    output$BT5<-renderUI({
        if (!(is.null(rv$Gene2_index)) & !(is.null(rv$logsel_index)) &
            !(is.null(rv$pvalsel_index)) & !(is.null(rv$GFAG_data2))) {
            actionButton('upfiltered','Update selected GFAG/Path')
        }
            else{
            shinyjs::disabled(actionButton('upfiltered',
                'Update selected GFAG/Path'))
            }
    })


    output$IDAT4<-renderText({
        if(!(is.null(rv$GFAG_data2))){
            return()
        }
            else{
                return('No GFAG data file was uploaded')
            }
    })


    GFAG_table2<-reactive({
        if(!(is.null(rv$GFAG_data2))){
            rv$GFAG_data2
        }
            else{
               return()
           }
    })


    GFAG_tableN<-reactive({
        if(!(is.null(rv$GFAG_data2))){
            dt<-rv$GFAG_data2
            dt<-as.data.frame(dt)
            dt1<-subset(dt,dt$qValue_n<input$qvalN)
            dt1
        }
            else{
                return()
           }
    })


    GFAG_tableH<-reactive({
        if(!(is.null(rv$GFAG_data2))){
            dt<-rv$GFAG_data2
            dt<-as.data.frame(dt)
            dt1<-subset(dt,dt$qValue_h<input$qvalH)
            dt1
        }
            else{
                return()
            }
    })


    GFAG_tableNH<-reactive({
        if(!(is.null(rv$GFAG_data2))){
            dt<-rv$GFAG_data2
            dt<-as.data.frame(dt)
            dt1<-subset(dt,dt$qValue_h<input$qvalH & dt$qValue_n<input$qvalN)
            dt1
        }
            else{
                return()
            }
    })


    rv<-reactiveValues()


    observeEvent(input$ViewN,{
        tt<-GFAG_tableN()
        rv$ftable<-tt
    })


    observeEvent(input$ViewH,{
        tt<-GFAG_tableH()
        rv$ftable<-tt
    })


    observeEvent(input$ViewNH,{
        tt<-GFAG_tableNH()
        rv$ftable<-tt
    })


    observeEvent(input$ViewR,{
        tt<-GFAG_table2()
        rv$ftable<-tt
    })


    output$ftable2<-DT::renderDataTable({
        req(rv$GFAG_data)
        req(rv$genekey)
        req(rv$logkey)
        req(rv$pvalkey)
        if (!(is.null(rv$GFAG_data)) & (rv$genekey==1 & rv$logkey==1 &
        rv$pvalkey==1)) {
            rv$ftable
        }else{
            NULL
        }

    })


    output$ftable3<-renderTable({
        rv$DEA_data2[1:10,]
    })


    output$ftable4<-renderTable({
        rv$DEA_data3[1:10,]
    })


    output$ftable5<-renderTable({
        rv$filtered_dt[1:10,]
    })


    output$ftable6<-renderTable({
        rv$mydt1[1:10,]
    })


    output$ftable7<-renderTable({
        rv$mydt2[1:10,]
    })


    output$ftable8<-renderTable({
        rv$mydt3[1:10,]
    })


    output$ftable9<-renderTable({
        rv$mydt4[1:10,]
    })


    output$ftable10<-renderTable({
        rv$mydt5[1:10,]
    })


    output$ftable11<-renderTable({
        rv$DEA_data3[1:10,]
    })


    output$ftable12<-renderTable({
        rv$data.m2
    })


    output$golist<-renderUI({
        if(!(is.null(rv$Gene2_index)) & !(is.null(rv$logsel_index)) &
           !(is.null(rv$pvalsel_index)) & !(is.null(rv$GFAG_data2)) 
           & !(is.null(rv$ftable))){
            selectInput(inputId = 'gol',
                label = 'Select a GFAG/Path of interest:',
                choices = rv$ftable[,1])
        }
            else{
                shinyjs::disabled(selectInput(inputId = 'gol',
                    label = 'Select a GFAG/Path of interest:',
                    choices = c()))
            }
    })


    output$golist2<-renderUI({
        if(!(is.null(rv$GFAG_data))){
            textInput(inputId = 'gol2',
                label = 'Write down a GFAG/Path of interest:')
        }
            else{
                shinyjs::disabled( textInput(inputId = 'gol2',
                    label = 'Write down a GFAG/Path of interest:'))
            }
    })


    observeEvent(input$gol,{
        rv$filtered_gfag<-NULL
        rv$CACO_cols<-NULL
        rv$QVAL_cols<-NULL
        rv$LOGFC_cols<-NULL
        rv$gg<-NULL
        rv$filtered_dt<-NULL
        rv$startfilter<-0
        rv$filtereddt<-0
    })


    observeEvent(input$upfiltered,{
        rv$filtered_gfag<-input$gol
    })


    rv$startfilter<-0


    output$Plot_to_start<-renderUI({
        if (!(is.null(rv$filtered_gfag))) {
            actionButton('Plot_start','Start Plotting')
        }
            else{
                shinyjs::disabled(actionButton('Plot_start',
                    'Start Plotting'))
            }
    })

    output$fdt_test1<-renderText({colnames((rv$DEA_data3))})

    observeEvent(input$Plot_start,{
        rv$startfilter<-rv$startfilter+1
        path_data<-rv$Path_data2
        colnames(path_data)[2]<-'GO.ID'
        colnames(path_data)[1]<-'Probe_Gene_ID'

        NEX_data3<-rv$NEX_data2[,c(1,as.numeric(grep(paste0('^',
            rv$Controlid,'$'),colnames(rv$NEX_data2))),as.numeric(
                grep(paste0('^',rv$Caseid,'$'),colnames(rv$NEX_data2))))]
        colnames(NEX_data3)[1]<-'Probe_Gene_ID'
        colnames(NEX_data3)[2]<-'Control'
        colnames(NEX_data3)[3]<-'Case'

        mydt1<-merge(path_data,NEX_data3,by.x='Probe_Gene_ID',all.x=T)

        mydt1<-mydt1[!(is.na(mydt1$GO.ID)),]

        DEA_data3<-rv$DEA_data2[,c(as.numeric(rv$Gene2_index),
            as.numeric(rv$logsel_index),as.numeric(rv$pvalsel_index))]
        
        colnames(DEA_data3)<-NA
        
        colnames(DEA_data3)[grep(colnames(
        rv$DEA_data2)[as.numeric(rv$Gene2_index)],
        colnames(rv$DEA_data2))]<-'Probe_Gene_ID'
        
        colnames(DEA_data3)[grep(colnames(
        rv$DEA_data2)[as.numeric(rv$logsel_index)],
        colnames(rv$DEA_data2))]<-'logFC'

        colnames(DEA_data3)[grep(colnames(
        rv$DEA_data2)[as.numeric(rv$pvalsel_index)],
        colnames(rv$DEA_data2))]<-'q.value'

        rv$DEA_data3<-DEA_data3
        mydt3<-merge(mydt1,DEA_data3,by.x='Probe_Gene_ID',all.x=T)
        mydt3<-mydt3[!(is.na(mydt3$GO.ID)),]

        rv$mydt5<-subset(mydt3,mydt3$GO.ID %in% rv$filtered_gfag)
        rv$mydt5<-rv$mydt5[!(is.na(rv$mydt5$GO.ID)),]

        rv$filtered_dt<-rv$mydt5

        data.m <- melt(rv$filtered_dt,
            id.vars=c('Probe_Gene_ID','logFC','GO.ID','q.value'))
        rv$data.m2<-data.m%>%
        mutate(G_logFC=ifelse(data.m$logFC > -1 & data.m$logFC < 1,
            '-1<logFC<1',ifelse(logFC <= -1,
            'logFC<=-1','logFC>=1')))%>%
            mutate(G_qvalue=ifelse(
                format(data.m$q.value,scientific = F) <= 0.05 &
                format(data.m$q.value,scientific = F) > 0.01,'*',
                ifelse(format(data.m$q.value,scientific = F) <= 0.01 &
                format(data.m$q.value,scientific = F) > 0.001,'**',
                ifelse(format(data.m$q.value,scientific = F) <= 0.001 ,
                '***','NS'))))
        sg<-rep(sign(rv$data.m2$value[rv$data.m2$variable=='Case']-
            rv$data.m2$value[rv$data.m2$variable=='Control']),2)
        rv$data.m2$Sign<-sg
    })

    output$headers1<-renderUI({
        if (!(is.null(rv$filtered_gfag)) & rv$startfilter>0) {
            radioButtons(inputId = 'selectCACO',
                label = 'Select a method for filtering Case & Control
                    values of repeated targets:',
                choices = c('None'='none',
                    'Case: Maximum Absolute Value'='max.ca.value',
                    'Case: Minimum Absolute Value'='min.ca.value',
                    'Control: Maximum Absolute Value'='max.co.value',
                    'Control: Minimum Absolute Value'='min.co.value',
                    'Variation: Maximum Absolute Value'='max.var',
                    'Variation: Minimum Absolute Value'='min.var',
                    'Take the median'='median.value'),
                selected = 'none')
        }
            else{
                if (rv$startfilter==0) {
                    shinyjs::disabled(radioButtons(inputId = 'selectCACO',
                        label = 'Select a method for filtering Case &
                            Control values of repeated targets:',
                        choices = c('None'='none',
                            'Case: Maximum Absolute Value'='max.ca.value',
                            'Case: Minimum Absolute Value'='min.ca.value',
                            'Control: Maximum Absolute Value'='max.co.value',
                            'Control: Minimum Absolute Value'='min.co.value',
                            'Variation: Maximum Absolute Value'='max.var',
                            'Variation: Minimum Absolute Value'='min.var',
                            'Take the median'='median.value'),
                        selected = 'none'))
                }
            }
    })


    rv$filtereddt<-0
    output$FilteredDT_Bt1<-renderUI({
        if (!(is.null(rv$filtered_gfag)) & rv$startfilter>0 &
            input$selectCACO!='none') {
            actionButton(inputId = 'filtered_dt',
                label = 'Create Filtered Dataset')
        }
            else{
                shinyjs::disabled(actionButton(inputId = 'filtered_dt',
                    label = 'Create Filtered Dataset'))
            }
    })


    actionButton(inputId = 'filtered_dt',label = 'Create Filtered Dataset')


    output$qval_sld<-renderUI({
        if (rv$startfilter>0) {
            sliderInput(inputId = 'sld_qval',
                label = 'Choose an inteval for q.values:',
                min = 0,
                max = 1,
                value = c(0,1))
        }
            else{
                shinyjs::disabled(sliderInput(inputId = 'sld_qval',
                    label = 'Choose an inteval for q.values:',
                    min = 0,
                    max = 1,
                    value = c(0,1)))
            }
    })


    output$logFC_sld<-renderUI({
        if (rv$startfilter>0) {
            sliderInput(inputId = 'sld_logfc',
                label = 'Choose an inteval for LogFC::',
                min = -10,
                max = 10,
                value = c(-2,2))
        }
            else{
            shinyjs::disabled(sliderInput(inputId = 'sld_logfc',
                label = 'Choose an inteval for LogFC::',
                min = -10,
                max = 10,
                value = c(-2,2)))
           }
    })


    output$height<-renderUI({
        if ( rv$startfilter>0) {
            textInput('height', 'Plot Height', value="500")
        }
            else{
                shinyjs::disabled(
                    textInput('height', 'Plot Height', value="500"))
            }
    })


    output$width<-renderUI({
        if ( rv$startfilter>0) {
            textInput('width', 'Width', value="100%")
        }
        else{
            shinyjs::disabled(textInput('width', 'Width', value="100%"))
        }
    })


    output$plotfont1<-renderUI({
        if ( rv$startfilter>0) {
            selectInput(inputId = 'Plot_font',
                'Select an .eps font type:',
                choices = names(postscriptFonts()),
                selected = "Bookman")
        }
            else{
                shinyjs::disabled(selectInput(
                    inputId = 'Plot_font','Select an .eps font type:',
                    choices = names(postscriptFonts()),selected = "Bookman"))
            }
    })


    output$fontsize1<-renderUI({
        if (rv$startfilter>0) {
            numericInput(inputId = 'Font_size1',
                label = 'Enter with a font size for plot elements:',
                min = 1,
                max = 20,
                value = 3 )
        }
            else{
                shinyjs::disabled(numericInput(inputId = 'Font_size1',
                    label = 'Enter with a font size for plot elements:',
                    min = 1,
                    max = 20,
                    value = 3 ))
            }
    })


    output$fontsize2<-renderUI({
        if (rv$startfilter>0) {
            numericInput(inputId = 'Font_size2',
                label = 'Enter with a font size for the plot title:',
                min = 1,
                max = 20,
                value = 12)
        }
            else{
                shinyjs::disabled(numericInput(inputId = 'Font_size2',
                    label = 'Enter with a font size for the plot title:',
                    min = 1,
                    max = 20,
                    value = 12))
            }
  })


    output$fontsize3<-renderUI({
        if (rv$startfilter>0) {
            numericInput(inputId = 'Font_size3',
                label = 'Enter with a font size for the axis text:',
                min = 1,
                max = 20,
                value = 10)
        }
            else{
                shinyjs::disabled(numericInput(inputId = 'Font_size3',
                    label = 'Enter with a font size for the axis text:',
                    min = 1,
                    max = 20,
                    value = 10))
            }
    })


    output$fontsize4<-renderUI({
        if (rv$startfilter>0) {
            numericInput(inputId = 'Font_size4',
                label = 'Enter with a font size for the axis title:',
                min = 1,
                max = 20,
                value = 10)
        }
            else{
                shinyjs::disabled(numericInput(inputId = 'Font_size4',
                    label = 'Enter with a font size for the axis title:',
                    min = 1,
                    max = 20,
                    value = 10))
            }
    })


    output$fontface1<-renderUI({
        if (rv$startfilter>0) {
            selectInput(inputId = 'font_face1',
                label = 'Select a font face for axis text:',
                choices = c('bold','italic','bold.italic','plain'),
                selected = 'plain')
        }
            else{
                shinyjs::disabled(selectInput(inputId = 'font_face1',
                    label = 'Select a font face for axis text:',
                    choices = c('bold','italic','bold.italic','plain'),
                    selected = 'plain'))
            }
    })


    output$fontface2<-renderUI({
        if (rv$startfilter>0) {
            selectInput(inputId = 'font_face2',
                label = 'Select a font face for axis title:',
                choices = c('bold','italic','bold.italic','plain'),
                selected = 'plain')
        }
            else{
                shinyjs::disabled(selectInput(inputId = 'font_face2',
                    label = 'Select a font face for axis title:',
                    choices = c('bold','italic','bold.italic','plain'),
                    selected = 'plain'))
            }
    })


    output$fontface3<-renderUI({
        if (rv$startfilter>0) {
            selectInput(inputId = 'font_face3',
                label = 'Select a font face for plot title:',
                choices = c('bold','italic','bold.italic','plain'),
                selected = 'plain')
        }
            else{
                shinyjs::disabled(selectInput(inputId = 'font_face3',
                    label = 'Select a font face for plot title:',
                    choices = c('bold','italic','bold.italic','plain'),
                    selected = 'plain'))
            }
    })


    output$logfcbt<-renderUI({
        if (rv$startfilter>0) {
            actionButton(inputId = 'plot_logfc',
                label = 'Plot with filtered LogFC values')
        }
            else{
                shinyjs::disabled(actionButton(inputId = 'plot_logfc',
                    label = 'Plot with filtered LogFC values'))
            }
    })


    output$qvalbt<-renderUI({
        if (rv$startfilter>0) {
            actionButton(inputId = 'plot_qval',
                label = 'Plot with filtered q.values')
        }
            else{
                shinyjs::disabled(actionButton(inputId = 'plot_qval',
                    label = 'Plot with filtered q.values'))
            }
    })


    output$qvlgbt<-renderUI({
        if (rv$startfilter>0) {
            actionButton(inputId = 'plot_qvlg',
                label = 'Plot with filtered q.values and LogFC values')
        }
            else{
                shinyjs::disabled(actionButton(inputId = 'plot_qvlg',
                    label = 'Plot with filtered q.values and LogFC values'))
            }
    })


    output$nfbt<-renderUI({
        if (rv$startfilter>0) {
            actionButton(inputId = 'plot_nf',label ='Plot with no filter')
        }
            else{
                shinyjs::disabled(actionButton(inputId = 'plot_nf',
                    label ='Plot with no filter'))
            }
    })


    output$ggtext<-renderText({
        rv$ggt
    })


    output$ggtext2<-renderText({
        rv$ggt2
    })


    output$fonttext<-renderText({
        as.character(input$Font_size1)
     })


    suppressWarnings(observeEvent(input$plot_nf,{
        if (rv$startfilter>0) {
            data.m <- melt(rv$filtered_dt,
                id.vars=c('Probe_Gene_ID','logFC','GO.ID','q.value'))
            data.m2<-data.m%>%
            mutate(G_logFC=ifelse(data.m$logFC > -1 & data.m$logFC < 1,
                '-1<logFC<1',ifelse(logFC <= -1,'logFC<=-1','logFC>=1')))%>%
                mutate(G_qvalue=ifelse(format(data.m$q.value,scientific = F) <=
                    0.05 & format(data.m$q.value,scientific = F) > 0.01,'*',
                    ifelse(format(data.m$q.value,scientific = F) <= 0.01 &
                    format(data.m$q.value,scientific = F) > 0.001,'**',
                    ifelse(format(data.m$q.value,scientific = F) <= 0.001 ,
                    '***','NS'))))
            sg<-rep(sign(data.m2$value[data.m2$variable=='Case']
                -data.m2$value[data.m2$variable=='Control']),2)
            data.m2$Sign<-sg
            gene_median<-(data.m2$value[data.m2$variable=='Case']+
                data.m2$value[data.m2$variable=='Control'])/2
            data.m2$gene_median<-gene_median
            dtca<-subset(data.m2,data.m2$variable=='Case')
            dtco<-subset(data.m2,data.m2$variable=='Control')
            data.m2<-rbind(dtco,dtca)
            rv$data.m2<-data.m2
            term<-rv$GFAG_data[,2][grep(input$gol,rv$GFAG_data[,1])]

            gg1<-suppressWarnings(eval(parse(text = paste0('ggplot(data.m2) +
                geom_path(aes(x = value, y = Probe_Gene_ID,color=G_logFC),
                arrow = arrow(length =  unit(1.5, "mm"), type = "closed")) +
                geom_text(aes(x = value, y = Probe_Gene_ID,
                label = round(value, 3)),
                hjust = ifelse(data.m2$Sign==-1,
                ifelse(data.m2$variable=="Case", 1.5, -0.4),
                ifelse(data.m2$variable=="Case", -0.4, 1.5)),
                family = "',input$Plot_font,'",
                size = ',input$Font_size1,',
                color = "gray25") +
                #---- nomes nos eixos, titulos, etc
                labs( x = "Expression",
                y = "Probe_Gene_ID",
                title = "Gene Expression Between Case and Control",
                subtitle = paste0("',input$gol,'"," ",term),
                caption = "Control= arrow back ; Case: arrow point") +
                # ----- Tema
                theme_bw() +
                theme(text = element_text(
                family = "',input$Plot_font,'", color = "gray25"),
                plot.subtitle = element_text(size = ',input$Font_size2,'),
                plot.caption = element_text(color = "gray30"),
                plot.background = element_rect(fill = "white"),
                plot.margin = unit(c(5,10,5,10), units = "mm"),
                axis.text = element_text(size = ',input$Font_size3,
                ',face = "',input$font_face1,'"),
                axis.title = element_text(size = ',input$Font_size4,
                ',face = "',input$font_face2,'"),
                title = element_text(size = ',input$Font_size2,
                ',face = "',input$font_face3,'"),
                legend.title=element_text(size=',input$Font_size3,'),
                legend.text=element_text(size=',input$Font_size3,'))  +
                # --- adicionando limites para a coordenada x
                coord_cartesian(xlim = c(min(data.m2$value)-0.7,
                max(data.m2$value)+0.7)) +
                scale_color_manual(
                name="logFC",
                values=c("-1<logFC<1"="black",
                "logFC>=1"="green",
                "logFC<=-1"="red"))'))))

            rv$ggt<-(paste0('ggplot(data.m2) +
                geom_path(aes(x = value, y = Probe_Gene_ID,color=G_logFC),
                arrow = arrow(length =  unit(1.5, "mm"), type = "closed")) +
                geom_text(aes(x = value, y = Probe_Gene_ID,
                label = round(value, 3)),
                hjust = ifelse(data.m2$Sign==-1,
                ifelse(data.m2$variable=="Case", 1.5, -0.4),
                ifelse(data.m2$variable=="Case", -0.4, 1.5)),
                family = "',input$Plot_font,'",
                size = ',input$Font_size1,',
                color = "gray25") +
                #---- nomes nos eixos, titulos, etc
                labs( x = "Expression",
                y = "Probe_Gene_ID",
                title = "Gene Expression Between Case and Control",
                subtitle = paste0("',input$gol,'"," ",term),
                caption = "Control= arrow back ; Case: arrow point") +
                # ----- Tema
                theme_bw() +
                theme(text = element_text(
                family = "',input$Plot_font,'", color = "gray25"),
                plot.subtitle = element_text(size = ',input$Font_size2,'),
                plot.caption = element_text(color = "gray30"),
                plot.background = element_rect(fill = "white"),
                plot.margin = unit(c(5,10,5,10), units = "mm"),
                axis.text = element_text(size = ',input$Font_size3,
                ',face = "',input$font_face1,'"),
                axis.title = element_text(size = ',input$Font_size4,
                ',face = "',input$font_face2,'"),
                title = element_text(size = ',input$Font_size2,
                ',face = "',input$font_face3,'"),
                legend.title=element_text(size=',input$Font_size3,'),
                legend.text=element_text(size=',input$Font_size3,'))  +
                # --- adicionando limites para a coordenada x
                coord_cartesian(xlim = c(min(data.m2$value)-0.7,
                max(data.m2$value)+0.7)) +
                scale_color_manual(
                name="logFC",
                values=c("-1<logFC<1"="black",
                "logFC>=1"="green",
                "logFC<=-1"="red"))'))

            gg2<-suppressWarnings(eval(parse(text ='gg1+
                geom_text(data = data.m2,
                aes(y = Probe_Gene_ID,
                x = gene_median,
                label = G_qvalue),
                vjust = 2,
                family = "Bookman",
                color = "gray25")')))

            rv$ggt2<-(('gg1+
                geom_text(data = data.m2,
                aes(y=Probe_Gene_ID,
                x = gene_median,
                label = G_qvalue),
                vjust = 2,
                family = "Bookman",
                color = "gray25")'))

            rv$gg<-gg2

            output[['myplot']]<-renderPlot({
                rv$gg
            })
        }
            else{
                rv$gg<-NULL
                rv$filtered<-NULL
                output[['myplot']]<-renderPlot({
                    rv$gg
                })
            }
    if (TRUE %in% ('try-error' %in% class(try(eval(parse(text=rv$ggt2)))))) {
        shinyjs::disabled('fplot')
    }
        else{
            shinyjs::enable('fplot')
        }
    }))


    observeEvent(input$plot_logfc,{
        if (rv$startfilter>0) {
            data.m <- melt(rv$filtered_dt, id.vars=c('Probe_Gene_ID','logFC',
                'GO.ID','q.value'))
            data.m2<-data.m%>%
                mutate(G_logFC=ifelse(data.m$logFC > -1 &
                    data.m$logFC < 1,'-1<logFC<1',
                    ifelse(logFC <= -1,'logFC<=-1','logFC>=1')))%>%
                        mutate(G_qvalue=ifelse(format(
                            data.m$q.value,scientific = F) <= 0.05 &
                            format(data.m$q.value,scientific = F) > 0.01,
                            '*',
                            ifelse(format(data.m$q.value,scientific = F)
                            <= 0.01 & format(data.m$q.value,scientific = F)
                            > 0.001,'**',
                            ifelse(format(data.m$q.value,scientific = F)
                            <= 0.001 ,'***','NS'))))
            sg<-rep(sign(data.m2$value[data.m2$variable=='Case']-
                data.m2$value[data.m2$variable=='Control']),2)
            data.m2$Sign<-sg
            gene_median<-(data.m2$value[data.m2$variable=='Case']+
                data.m2$value[data.m2$variable=='Control'])/2
            data.m2$Sign<-sg
            gene_median<-(data.m2$value[data.m2$variable=='Case']+
                data.m2$value[data.m2$variable=='Control'])/2
            data.m2$gene_median<-gene_median
            dtca<-subset(data.m2,data.m2$variable=='Case')
            dtco<-subset(data.m2,data.m2$variable=='Control')
            data.m2<-rbind(dtco,dtca)
            rv$data.m2<-data.m2
            data.m2<-subset(data.m2,(data.m2$logFC<=max(input$sld_logfc)) &
                data.m2$logFC>=min(input$sld_logfc))
            term<-rv$GFAG_data[,2][grep(input$gol,rv$GFAG_data[,1])]


            gg1<-eval(parse(text = paste0('ggplot(data.m2) +
                geom_path(aes(x = value, y = Probe_Gene_ID,color=G_logFC),
                arrow = arrow(length =  unit(1.5, "mm"), type = "closed")) +
                geom_text(aes(x = value, y = Probe_Gene_ID,
                label = round(value, 3)),
                hjust = ifelse(data.m2$Sign==-1,
                ifelse(data.m2$variable=="Case", 1.5, -0.4),
                ifelse(data.m2$variable=="Case", -0.4, 1.5)),
                family = "',input$Plot_font,'",
                size = ',input$Font_size1,',
                color = "gray25") +
                #---- nomes nos eixos, titulos, etc
                labs( x = "Expression",
                y = "Probe_Gene_ID",
                title = "Gene Expression Between Case and Control",
                subtitle = paste0("',input$gol,'"," ",term),
                caption = "Control= arrow back ; Case: arrow point") +
                # ----- Tema
                theme_bw() +
                theme(text = element_text(
                family = "',input$Plot_font,'", color = "gray25"),
                plot.subtitle = element_text(size = ',input$Font_size2,'),
                plot.caption = element_text(color = "gray30"),
                plot.background = element_rect(fill = "white"),
                plot.margin = unit(c(5,10,5,10), units = "mm"),
                axis.text = element_text(size = ',input$Font_size3,
                ',face = "',input$font_face1,'"),
                axis.title = element_text(size = ',input$Font_size4,
                ',face = "',input$font_face2,'"),
                title = element_text(size = ',input$Font_size2,
                ',face = "',input$font_face3,'"),
                legend.title=element_text(size=',input$Font_size3,'),
                legend.text=element_text(size=',input$Font_size3,'))  +
                # --- adicionando limites para a coordenada x
                coord_cartesian(xlim = c(min(data.m2$value)-0.7,
                max(data.m2$value)+0.7)) +
                scale_color_manual(
                name="logFC",
                values=c("-1<logFC<1"="black",
                "logFC>=1"="green",
                "logFC<=-1"="red"))')))

            rv$ggt<-(paste0('ggplot(data.m2) +
                geom_path(aes(x = value, y = Probe_Gene_ID,color=G_logFC),
                arrow = arrow(length =  unit(1.5, "mm"), type = "closed")) +
                geom_text(aes(x = value, y = Probe_Gene_ID,
                label = round(value, 3)),
                hjust = ifelse(data.m2$Sign==-1,
                ifelse(data.m2$variable=="Case", 1.5, -0.4),
                ifelse(data.m2$variable=="Case", -0.4, 1.5)),
                family = "',input$Plot_font,'",
                size = ',input$Font_size1,',
                color = "gray25") +
                #---- nomes nos eixos, titulos, etc
                labs( x = "Expression",
                y = "Probe_Gene_ID",
                title = "Gene Expression Between Case and Control",
                subtitle = paste0("',input$gol,'"," ",term),
                caption = "Control= arrow back ; Case: arrow point") +
                # ----- Tema
                theme_bw() +
                theme(text = element_text(
                family = "',input$Plot_font,'", color = "gray25"),
                plot.subtitle = element_text(size = ',input$Font_size2,'),
                plot.caption = element_text(color = "gray30"),
                plot.background = element_rect(fill = "white"),
                plot.margin = unit(c(5,10,5,10), units = "mm"),
                axis.text = element_text(size = ',input$Font_size3,
                ',face = "',input$font_face1,'"),
                axis.title = element_text(size = ',input$Font_size4,
                ',face = "',input$font_face2,'"),
                title = element_text(size = ',input$Font_size2,
                ',face = "',input$font_face3,'"),
                legend.title=element_text(size=',input$Font_size3,'),
                legend.text=element_text(size=',input$Font_size3,'))  +
                # --- adicionando limites para a coordenada x
                coord_cartesian(xlim = c(min(data.m2$value)-0.7,
                max(data.m2$value)+0.7)) +
                scale_color_manual(
                name="logFC",
                values=c("-1<logFC<1"="black",
                "logFC>=1"="green",
                "logFC<=-1"="red"))'))

            gg2<-eval(parse(text ='gg1+
                geom_text(data = data.m2,
                aes(y = Probe_Gene_ID,
                x = gene_median,
                label = G_qvalue),
                vjust = 2,
                family = "Bookman",
                color = "gray25")'))

            rv$ggt2<-(('gg1+
                geom_text(data = data.m2,
                aes(y=Probe_Gene_ID,
                x = gene_median,
                label = G_qvalue),
                vjust = 2,
                family = "Bookman",
                color = "gray25")'))

            rv$gg<-gg2

            output[['myplot']]<-renderPlot({
                rv$gg
            })
        }
            else{
                rv$gg<-NULL
                rv$filtered<-NULL
                output[['myplot']]<-renderPlot({
                rv$gg
                })}
    if (TRUE %in% ('try-error' %in% class(try(eval(parse(text=rv$ggt2)))))) {
        shinyjs::disabled('fplot')
    }
        else{
            shinyjs::enable('fplot')
        }
    })


    observeEvent(input$plot_qval,{
       if (rv$startfilter>0) {
            data.m <- melt(rv$filtered_dt,
                id.vars=c('Probe_Gene_ID','logFC','GO.ID','q.value'))
            data.m2<-data.m%>%
            mutate(G_logFC=ifelse(data.m$logFC > -1 &
                data.m$logFC < 1,'-1<logFC<1',ifelse(logFC <= -1,
                'logFC<=-1','logFC>=1')))%>%
                    mutate(G_qvalue=ifelse(format(
                        data.m$q.value,scientific = F) <= 0.05 &
                        format(data.m$q.value,scientific = F) > 0.01,
                        '*',
                        ifelse(format(data.m$q.value,scientific = F) <= 0.01 &
                        format(data.m$q.value,scientific = F) > 0.001,'**',
                        ifelse(format(data.m$q.value,scientific = F) <= 0.001 ,
                        '***','NS'))))
            sg<-rep(sign(data.m2$value[data.m2$variable=='Case']-
                data.m2$value[data.m2$variable=='Control']),2)
            data.m2$Sign<-sg
            gene_median<-(data.m2$value[data.m2$variable=='Case']+
            data.m2$value[data.m2$variable=='Control'])/2
            data.m2$gene_median<-gene_median
            dtca<-subset(data.m2,data.m2$variable=='Case')
            dtco<-subset(data.m2,data.m2$variable=='Control')
            data.m2<-rbind(dtco,dtca)
            rv$data.m2<-data.m2
            data.m2<-subset(data.m2,(data.m2$q.value<=max(input$sld_qval)) &
                data.m2$q.value>=min(input$sld_qval))
            term<-rv$GFAG_data[,2][grep(input$gol,rv$GFAG_data[,1])]


            gg1<-eval(parse(text = paste0('ggplot(data.m2) +
                geom_path(aes(x = value, y = Probe_Gene_ID,color=G_logFC),
                arrow = arrow(length =  unit(1.5, "mm"), type = "closed")) +
                geom_text(aes(x = value, y = Probe_Gene_ID,
                label = round(value, 3)),
                hjust = ifelse(data.m2$Sign==-1,
                ifelse(data.m2$variable=="Case", 1.5, -0.4),
                ifelse(data.m2$variable=="Case", -0.4, 1.5)),
                family = "',input$Plot_font,'",
                size = ',input$Font_size1,',
                color = "gray25") +
                #---- nomes nos eixos, titulos, etc
                labs( x = "Expression",
                y = "Probe_Gene_ID",
                title = "Gene Expression Between Case and Control",
                subtitle = paste0("',input$gol,'"," ",term),
                caption = "Control= arrow back ; Case: arrow point") +
                # ----- Tema
                theme_bw() +
                theme(text = element_text(
                family = "',input$Plot_font,'", color = "gray25"),
                plot.subtitle = element_text(size = ',input$Font_size2,'),
                plot.caption = element_text(color = "gray30"),
                plot.background = element_rect(fill = "white"),
                plot.margin = unit(c(5,10,5,10), units = "mm"),
                axis.text = element_text(size = ',input$Font_size3,
                ',face = "',input$font_face1,'"),
                axis.title = element_text(size = ',input$Font_size4,
                ',face = "',input$font_face2,'"),
                title = element_text(size = ',input$Font_size2,
                ',face = "',input$font_face3,'"),
                legend.title=element_text(size=',input$Font_size3,'),
                legend.text=element_text(size=',input$Font_size3,'))  +
                # --- adicionando limites para a coordenada x
                coord_cartesian(xlim = c(min(data.m2$value)-0.7,
                max(data.m2$value)+0.7)) +
                scale_color_manual(
                name="logFC",
                values=c("-1<logFC<1"="black",
                "logFC>=1"="green",
                "logFC<=-1"="red"))')))

            rv$ggt<-(paste0('ggplot(data.m2) +
                geom_path(aes(x = value, y = Probe_Gene_ID,color=G_logFC),
                arrow = arrow(length =  unit(1.5, "mm"), type = "closed")) +
                geom_text(aes(x = value, y = Probe_Gene_ID,
                label = round(value, 3)),
                hjust = ifelse(data.m2$Sign==-1,
                ifelse(data.m2$variable=="Case", 1.5, -0.4),
                ifelse(data.m2$variable=="Case", -0.4, 1.5)),
                family = "',input$Plot_font,'",
                size = ',input$Font_size1,',
                color = "gray25") +
                #---- nomes nos eixos, titulos, etc
                labs( x = "Expression",
                y = "Probe_Gene_ID",
                title = "Gene Expression Between Case and Control",
                subtitle = paste0("',input$gol,'"," ",term),
                caption = "Control= arrow back ; Case: arrow point") +
                # ----- Tema
                theme_bw() +
                theme(text = element_text(
                family = "',input$Plot_font,'", color = "gray25"),
                plot.subtitle = element_text(size = ',input$Font_size2,'),
                plot.caption = element_text(color = "gray30"),
                plot.background = element_rect(fill = "white"),
                plot.margin = unit(c(5,10,5,10), units = "mm"),
                axis.text = element_text(size = ',input$Font_size3,
                ',face = "',input$font_face1,'"),
                axis.title = element_text(size = ',input$Font_size4,
                ',face = "',input$font_face2,'"),
                title = element_text(size = ',input$Font_size2,
                ',face = "',input$font_face3,'"),
                legend.title=element_text(size=',input$Font_size3,'),
                legend.text=element_text(size=',input$Font_size3,'))  +
                # --- adicionando limites para a coordenada x
                coord_cartesian(xlim = c(min(data.m2$value)-0.7,
                max(data.m2$value)+0.7)) +
                scale_color_manual(
                name="logFC",
                values=c("-1<logFC<1"="black",
                "logFC>=1"="green",
                "logFC<=-1"="red"))'))

            gg2<-eval(parse(text ='gg1+
                geom_text(data = data.m2,
                aes(y = Probe_Gene_ID,
                x = gene_median,
                label = G_qvalue),
                vjust = 2,
                family = "Bookman",
                color = "gray25")'))

            rv$ggt2<-(('gg1+
                geom_text(data = data.m2,
                aes(y=Probe_Gene_ID,
                x = gene_median,
                label = G_qvalue),
                vjust = 2,
                family = "Bookman",
                color = "gray25")'))
            rv$gg<-gg2
            output[['myplot']]<-renderPlot({
                rv$gg
            })
        }
            else{
                rv$gg<-NULL
                rv$filtered<-NULL
                output[['myplot']]<-renderPlot({
                    rv$gg
                })
            }
    if (TRUE %in% ('try-error' %in% class(try(eval(parse(text=rv$ggt2)))))) {
        shinyjs::disabled('fplot')
    }
        else{
            shinyjs::enable('fplot')
        }
    })


    observeEvent(input$plot_qvlg,{
        if (rv$startfilter>0) {
            data.m <- melt(rv$filtered_dt, id.vars=c('Probe_Gene_ID',
                'logFC','GO.ID','q.value'))
            data.m2<-data.m%>%
            mutate(G_logFC=ifelse(data.m$logFC > -1 & data.m$logFC < 1,
                '-1<logFC<1',ifelse(logFC <= -1,'logFC<=-1','logFC>=1')))%>%
                 mutate(G_qvalue=ifelse(format(data.m$q.value,scientific = F)
                    <= 0.05 & format(data.m$q.value,scientific = F) > 0.01,'*',
                    ifelse(format(data.m$q.value,scientific = F) <= 0.01 &
                    format(data.m$q.value,scientific = F) > 0.001,'**',
                    ifelse(format(data.m$q.value,scientific = F) <= 0.001 ,
                    '***','NS'))))
            sg<-rep(sign(data.m2$value[data.m2$variable=='Case']-
            data.m2$value[data.m2$variable=='Control']),2)
            data.m2$Sign<-sg
            gene_median<-(data.m2$value[data.m2$variable=='Case']+
                data.m2$value[data.m2$variable=='Control'])/2
            data.m2$gene_median<-gene_median
            dtca<-subset(data.m2,data.m2$variable=='Case')
            dtco<-subset(data.m2,data.m2$variable=='Control')
            data.m2<-rbind(dtco,dtca)
            rv$data.m2<-data.m2
            data.m2<-subset(data.m2,(((data.m2$q.value<=max(input$sld_qval)) &
                data.m2$q.value>=min(input$sld_qval)) &
                ((data.m2$logFC<=max(input$sld_logfc)) &
                data.m2$logFC>=min(input$sld_logfc))))
            term<-rv$GFAG_data[,2][grep(input$gol,rv$GFAG_data[,1])]


            gg1<-eval(parse(text = paste0('ggplot(data.m2) +
                geom_path(aes(x = value, y = Probe_Gene_ID,color=G_logFC),
                arrow = arrow(length =  unit(1.5, "mm"), type = "closed")) +
                geom_text(aes(x = value, y = Probe_Gene_ID,
                label = round(value, 3)),
                hjust = ifelse(data.m2$Sign==-1,
                ifelse(data.m2$variable=="Case", 1.5, -0.4),
                ifelse(data.m2$variable=="Case", -0.4, 1.5)),
                family = "',input$Plot_font,'",
                size = ',input$Font_size1,',
                color = "gray25") +
                #---- nomes nos eixos, titulos, etc
                labs( x = "Expression",
                y = "Probe_Gene_ID",
                title = "Gene Expression Between Case and Control",
                subtitle = paste0("',input$gol,'"," ",term),
                caption = "Control= arrow back ; Case: arrow point") +
                # ----- Tema
                theme_bw() +
                theme(text = element_text(
                family = "',input$Plot_font,'", color = "gray25"),
                plot.subtitle = element_text(size = ',input$Font_size2,'),
                plot.caption = element_text(color = "gray30"),
                plot.background = element_rect(fill = "white"),
                plot.margin = unit(c(5,10,5,10), units = "mm"),
                axis.text = element_text(size = ',input$Font_size3,
                ',face = "',input$font_face1,'"),
                axis.title = element_text(size = ',input$Font_size4,
                ',face = "',input$font_face2,'"),
                title = element_text(size = ',input$Font_size2,
                ',face = "',input$font_face3,'"),
                legend.title=element_text(size=',input$Font_size3,'),
                legend.text=element_text(size=',input$Font_size3,'))  +
                # --- adicionando limites para a coordenada x
                coord_cartesian(xlim = c(min(data.m2$value)-0.7,
                max(data.m2$value)+0.7)) +
                scale_color_manual(
                name="logFC",
                values=c("-1<logFC<1"="black",
                "logFC>=1"="green",
                "logFC<=-1"="red"))')))

            rv$ggt<-(paste0('ggplot(data.m2) +
                geom_path(aes(x = value, y = Probe_Gene_ID,color=G_logFC),
                arrow = arrow(length =  unit(1.5, "mm"), type = "closed")) +
                geom_text(aes(x = value, y = Probe_Gene_ID,
                label = round(value, 3)),
                hjust = ifelse(data.m2$Sign==-1,
                ifelse(data.m2$variable=="Case", 1.5, -0.4),
                ifelse(data.m2$variable=="Case", -0.4, 1.5)),
                family = "',input$Plot_font,'",
                size = ',input$Font_size1,',
                color = "gray25") +
                #---- nomes nos eixos, titulos, etc
                labs( x = "Expression",
                y = "Probe_Gene_ID",
                title = "Gene Expression Between Case and Control",
                subtitle = paste0("',input$gol,'"," ",term),
                caption = "Control= arrow back ; Case: arrow point") +
                # ----- Tema
                theme_bw() +
                theme(text = element_text(
                family = "',input$Plot_font,'", color = "gray25"),
                plot.subtitle = element_text(size = ',input$Font_size2,'),
                plot.caption = element_text(color = "gray30"),
                plot.background = element_rect(fill = "white"),
                plot.margin = unit(c(5,10,5,10), units = "mm"),
                axis.text = element_text(size = ',input$Font_size3,
                ',face = "',input$font_face1,'"),
                axis.title = element_text(size = ',input$Font_size4,
                ',face = "',input$font_face2,'"),
                title = element_text(size = ',input$Font_size2,
                ',face = "',input$font_face3,'"),
                legend.title=element_text(size=',input$Font_size3,'),
                legend.text=element_text(size=',input$Font_size3,'))  +
                # --- adicionando limites para a coordenada x
                coord_cartesian(xlim = c(min(data.m2$value)-0.7,
                max(data.m2$value)+0.7)) +
                scale_color_manual(
                name="logFC",
                values=c("-1<logFC<1"="black",
                "logFC>=1"="green",
                "logFC<=-1"="red"))'))

            gg2<-eval(parse(text ='gg1+
                geom_text(data = data.m2,
                aes(y = Probe_Gene_ID,
                x = gene_median,
                label = G_qvalue),
                vjust = 2,
                family = "Bookman",
                color = "gray25")'))

            rv$ggt2<-(('gg1+
                geom_text(data = data.m2,
                aes(y=Probe_Gene_ID,
                x = gene_median,
                label = G_qvalue),
                vjust = 2,
                family = "Bookman",
                color = "gray25")'))
            rv$gg<-gg2

            rv$dwld_plot<-gg2

            output[['myplot']]<-renderPlot({
                rv$gg
            })
        }
            else{
                rv$gg<-NULL
                rv$filtered<-NULL
                output[['myplot']]<-renderPlot({
                    rv$gg
                })
            }
    if (TRUE %in% ('try-error' %in% class(try(eval(parse(text=rv$ggt2)))))) {
        shinyjs::disabled('fplot')
    }
        else{
            shinyjs::enable('fplot')
        }
    })


    final_plot<-reactive({
        rv$gg
    })


    output$crplot<-renderUI({
        plotOutput('myplot',
            height = validateCssUnit(input$height),
            width = validateCssUnit(input$width))
    })


    output$fplot<-downloadHandler(
        filename = paste0('Plot_',Sys.Date(),'.eps'),
        content = function(file) {
        postscript(file,family = input$Plot_font)
        print(final_plot())
        #ggsave(file, plot = final_plot(), device = "eps")
        dev.off()
        })


    output$test7<-renderTable({
        head(rv$data.m2,100)
    })


 }

shinyApp(ui, server)
