suppressMessages(require(shiny))
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
suppressMessages(require(DT))
suppressMessages(require(data.table))
suppressMessages(require(gridExtra))


makePlotContainers <- function(n=1, ncol=1, prefix="plot", height=100,
width="100%", ...) {
    ## Validate inputs
    validateCssUnit(width)
    validateCssUnit(height)

    ## Construct plotOutputs
    lst <- lapply(seq.int(n), function(i)
    plotOutput(sprintf('%s_%g', prefix, i), height=height, width=width))

    ## Make columns
    lst <- lapply(split(lst, (seq.int(n))), function(x) column(ncol, x))
    do.call(tagList, lst)
}



ui <- fluidPage(
    shinyjs::useShinyjs(),
    br(),
    # Sidebar layout with a input and output definitions
    sidebarLayout(
    # Inputs
        sidebarPanel(
        # added ID which will be used for a Shiny input value,
        # and it will report which tab is selected.
            tabsetPanel(id='tabs',
                tabPanel('Data Input',value = 1,
                    br(),
                    fileInput('file','Upload non-formated data'),
                    tags$hr(),
                    #helpText("Default max. file size is 5MB"),
                    helpText("Select the read.table parameters below:"),
                    checkboxInput(inputId = 'header', label = 'Header',
                        value = TRUE),
                    checkboxInput(inputId = "stringAsFactors",
                        "stringAsFactors", FALSE),
                    tags$hr(),
                    helpText('Select the separator type:'),
                    radioButtons(inputId = 'sep', label = 'Separator',
                        choices = c(Comma=',',
                    Semicolon=';',Tab='\t', Space=''), selected = '\t')),
                tabPanel('Data Selection',value = 2,
                    br(),
                    uiOutput('ID_column_select'),
                    uiOutput('Pval_column_select'),
                    helpText('Use the two dropdown menus above to select/create
                        your dataframe.'),
                    tags$hr(),
                    actionButton('dt_creation','Create Dataframe'),
                    helpText('Click in the buttom above to store the currently
                        \'Selected Dataframe\'. Note that you can add content
                        from multiple files, just select the wanted data from
                        the currently selected file, hit the
                        \'Create Dataframe\' buttom and then select another
                        filein the \'Data Input\' tab!'),
                    tags$hr(),
                    actionButton('dt_confirmation',
                        'Dataframe is ready for formating!'),
                    helpText('Click in the buttom above if the
                        \'Created Dataframe\' is ready for plotting.'),
                    tags$hr(),
                    actionButton('Reset_dt_creation','Reset Dataframe'),
                    helpText('Click in the buttom above to reset the
                        \'Created Dataframe\'.')),
                tabPanel('Data Formating',value = 3,
                    br(),
                    uiOutput('Column_order'),
                    helpText('It is recommended to always put your ID column in
                        the first position!'),
                    actionButton('Reorder_columns','Re-order Dataframe'),
                    tags$hr(),
                    uiOutput('Column_rename_selection'),
                    br(),
                    uiOutput('Column_rename_write'),
                    br(),
                    actionButton('Confirm_new_name','Rename!'),
                    helpText('Note that renaming columns is not obrigatory,
                        however since the column names will be used in the plot
                        it is recommended that the columns
                        are properly named.'),
                    tags$hr(),
                    uiOutput('ID_column_position'),
                    tags$hr(),
                    actionButton('Formating_end_confirmation',
                        'Dataframe is ready plotting!')),
                tabPanel('Plot',value=5,
                    br(),
                    textInput('Plot_Height', 'Plot Height', value="500"),
                    textInput('Plot_Width', 'Plot Width', value="100%"),
                    uiOutput('Plot_element_font'),
                    uiOutput('Plot_element_size'),
                    uiOutput('Plot_legend_size'),
                    uiOutput('Plot_axis_text_face'),
                    uiOutput('Plot_legend_face'),
                    actionButton('Plot_Graph','Plot it!'),
                    downloadButton("Download_Plot_Buttom", "Download Plot")))),
        mainPanel(
            conditionalPanel(condition = 'input.tabs==1',
                h3('File Information:'),
                textOutput('Tab_1_main_panel_progress_info'),
                tags$hr(),
                uiOutput('Tab_1_main_panel_ui')),
            conditionalPanel(condition = 'input.tabs==2',
                h3('Selected Dataframe:'),
                uiOutput('Tab_2_main_panel_disabled_warning'),
                uiOutput('Selected_dt'),
                tags$hr(),
                h3('Created Dataframe:'),
                uiOutput('Created_dt'),
                uiOutput('Tab_2_main_panel_progress_info'),
                tags$hr()),
            conditionalPanel(condition = 'input.tabs==3',
                h3('Your currently created dataframe:'),
                uiOutput('Created_dt_2'),
                tags$hr(),h3('Your currently formated dataframe:'),
                tableOutput('Formated_dt'),
                tags$hr(),
                uiOutput('Tab_3_main_panel_progress_info')),
            conditionalPanel(condition = 'input.tabs==5',
                h3('Your current plot:'),
                uiOutput('Plot_obj'))
    )
  ))

server <- function(input, output, session) {
    rv<-reactiveValues()

    shinyjs::disable('dt_creation')
    shinyjs::disable('dt_confirmation')
    shinyjs::disable('Reset_dt_creation')
    shinyjs::disable('bars')
    shinyjs::disable('ncols')
    shinyjs::disable('Confirm_new_name')
    shinyjs::disable('Reorder_columns')
    shinyjs::disable('Formating_end_confirmation')
    shinyjs::disable('upval')
    shinyjs::disable('Plot_Height')
    shinyjs::disable('Plot_Width')
    shinyjs::disable('Plot_Graph')
    shinyjs::disable('Download_Plot_Buttom')

    ###########################################################################
                              #### TAB=1 ####
    data <- reactive({

        if(is.null(input$file)) {return()
         }
            else{
                shinyjs::enable('dt_creation')
                input_file_dt<-try(read.delim(file=input$file$datapath,
                    sep=input$sep,
                    header = input$header,
                    stringsAsFactors = input$stringAsFactors),silent = TRUE)
                if (class(input_file_dt)=='data.frame') {
                    rv$GFAG_data1<-input_file_dt
                    rv$GFAG_data2<-input_file_dt
                }
                    else{
                        rv$GFAG_data1<-NULL
                        rv$GFAG_data2<-NULL
                    }
            }
    })



    # this reactive output contains the summary of the dataset and display
    # the summary in table format
    output$Tab_1_main_panel_file_info <- renderTable({
        if(is.null(data())){return ()
        }
        input$file
    })

    # this reactive output contains the summary of the dataset and display
    # the summary in table format
    output$Tab_1_main_panel_Summary <- renderTable({
        if(is.null(data())){return ()
        }
        summary(data())
    })

    # This reactive output contains the dataset and display the dataset in
    # table format
    output$Tab_1_main_panel_dt <- DT::renderDataTable({
        if(is.null(data())){
            return (NULL)
            }
                else{
                    data()
                }
    })

    # the following renderUI is used to dynamically generate the tabsets when
    # the file is loaded. Until the file is loaded,
    # app will not show the tabset.
    output$Tab_1_main_panel_ui <- renderUI({
        if(is.null(data())){
            return()
        }
        else{tabsetPanel(tabPanel("About file",
            tableOutput("Tab_1_main_panel_file_info")),
            tabPanel("Data",
            DT::dataTableOutput("Tab_1_main_panel_dt")),
            tabPanel("Summary", tableOutput("Tab_1_main_panel_Summary")))
        }
    })


    Tab_1_main_panel_progress_info<-reactive({
        if(is.null(input$file)) {
        return('No file was uploaded! We recomend a maximun of 50 rows!')
        }
            else{
                input_file_dt<-try(read.delim(file=input$file$datapath,
                    sep=input$sep,
                    header = input$header,
                    stringsAsFactors = input$stringAsFactors))
                if (class(input_file_dt)=='data.frame' &
                    length(rownames(input_file_dt))<=50) {
                    return('The file was sucessfully uploaded! Use the "Data"
                           tab to check if your file is correctly formated!')
                }
                    else{
                        if (class(input_file_dt)=='data.frame' &
                            length(rownames(input_file_dt))>50) {
                            return('The file was sucessfully uploaded! Use the
                                   "Data" tab to check if your file is
                                   correctly formated! (PS: the uploaded
                                   file contains more than 50 rows)')
                        }
                            else{
                                if (class(input_file_dt)=='try-error') {
                                    return('The uploaded file is
                                           incompatible with the selected
                                           separator type. Please select
                                           another separator.')
                                }
                            }
                    }
            }
    })



    output$Tab_1_main_panel_progress_info<-renderText({
        return(Tab_1_main_panel_progress_info())
    })


    ###########################################################################
                               ### TAB=2 ###


    output$ID_column_select<-renderUI({
        if (is.null(input$file)) {
            shinyjs::disabled(selectInput('ID_col','Select the id column',
            choices = colnames(data()),multiple = FALSE, selected = NULL))
        }
            else{
                if (!(is.null(input$file))) {
                    selectInput('ID_col','Select the id column',
                        choices = colnames(data()),
                        multiple = FALSE,selected = NULL)
                }
            }
    })


    output$Pval_column_select<-renderUI({
        if (is.null(input$file)) {
            shinyjs::disabled(selectInput('Pval_col',
            'Select the p-values columns',
            choices = colnames(data()),
            multiple = TRUE,
            selected = NULL))
        }
            else{
                if (!(is.null(input$file))) {
                    selectInput('Pval_col',
                    'Select the p-values columns',
                    choices = colnames(data()),
                    multiple = TRUE,
                    selected = NULL)
                }
            }
    })


    Selected_dt1<-reactive({
        # If missing input, return to avoid error later in function
        if(is.null(input$file)){
            return()
        }

        # Get the data set
        dat <- data()

        # Make sure columns are correct for data set (when data set changes,
        #the columns will initially be for the previous data set)

        if (is.null(c(input$ID_col,input$Pval_col)) ||
            !(c(input$ID_col,input$Pval_col) %in% names(dat))){
            return()
        }

        # Keep the selected columns
        dat <- dat[, c(input$ID_col,input$Pval_col), drop = FALSE]
    })


    output$Tab_2_main_panel_disabled_warning<-renderUI({
        if(is.null(input$file)){
            h5('You cannot select your dataframe while there is no
                file uploaded!')
        }
            else{
                return()
            }
    })


    output$Selected_dt<-renderTable({
        head(Selected_dt1(), 20)
    })


    counter<-reactiveValues(countervalue=0)


    observeEvent(input$dt_creation,{
        if(counter2$countervalue2==1){
            counter$countervalue <- counter$countervalue + 1
        }
    })


    rv<-reactiveValues()


    pvals_c<-reactive({
        # If missing input, return to avoid error later in function
        if(is.null(input$file)){
            return()
        }

        # Get the data set
        dat <- data()

        # Make sure columns are correct for data set (when data set changes,
        # the columns will initially be for the previous data set)
        if (is.null(c(input$Pval_col)) || !(c(input$Pval_col) %in%
                                            names(dat))){
            return()
        }

        # Keep the selected columns
        dat <- dat[, c(input$Pval_col), drop = FALSE]
    })


    observeEvent(input$dt_creation, {
        if(counter$countervalue==1){
            rv$data <- Selected_dt1()
        }
            else{
                if (!(is.null(pvals_c()))) {
                    rv$data<-cbind(rv$data,pvals_c())
                }
            }

        if (!(is.null(rv$data)) & length(colnames(rv$data))>=2) {
            shinyjs::enable('dt_confirmation')
            shinyjs::enable('Reset_dt_creation')
        }
    })


    counter2<-reactiveValues(countervalue2=0)


    observeEvent(input$file,{
        counter2$countervalue2 <- 1
    })


    counter3<-reactiveValues(countervalue3=0)


    observeEvent(input$dt_confirmation,{
        if(counter2$countervalue2==1 & counter$countervalue>0){
            counter3$countervalue3 <- 1
        }
    })


    output$Tab_2_main_panel_progress_info<-renderUI({
        if (counter2$countervalue2==0 & counter$countervalue==0 |
            (counter2$countervalue2==0 & counter$countervalue>0)) {
            h5('You cannot create a dataframe while
               there is no file uploaded!')
        }
            else{
                if ( (counter2$countervalue2>0 & counter$countervalue==0)) {
                    h5('A file was uploaded, now select the desired
                        dataframe and hit the \'Create Dataframe\'
                        buttom to create your dataframe. Note that
                        in order to sucessfully create a dataframe
                        you must select at least one column besides
                        the ID column.
                        When your created dataframe is finished hit
                        the \'Dataframe is ready for formating!\' buttom!')
                }
                    else{
                        if (counter2$countervalue2>0 & !(is.null(rv$data))
                            & (is.null(rv$data_to_plot)) &
                            length(colnames(rv$data))<2) {
                            h5('To create a dataframe you must select at least
                            one column besides the id column')
                        }
                            else{
                                if (counter2$countervalue2>0 &
                                    !(is.null(rv$data))
                                    & (is.null(rv$data_to_plot)) &
                                    length(colnames(rv$data))>=2) {
                                    h5('You created the dataframe, now hit the
                                    \'Dataframe is ready for formating!\'
                                    buttom if it is finished')
                                }
                                    else{
                                        if (counter2$countervalue2>0 &
                                        !(is.null(rv$data)) &
                                        !(is.null(rv$data_to_plot))) {
                                            h5('You are all done here! Go to
                                            the \'Data Formating\' tab to
                                            format your dataframe')
                                        }
                                    }

                            }
                    }
            }
    })


    output$Created_dt<-renderTable({
        if (!(is.null(rv$data)) & length(colnames(rv$data))>=2) {
            rv$data
        }
            else{
                NULL
            }
    })


    output$Created_dt_2<-renderTable({
        if(!(is.null(rv$data)) & counter3$countervalue3>0){
            rv$data
        }
    })


    observeEvent(input$Reset_dt_creation, {
        rv$data<-NULL
        rv$data_to_plot<-NULL
        rv$formated_data<-NULL
        counter$countervalue<-0
        counter3$countervalue3<-0
        counter4$countervalue4<-0
        counter5$countervalue5<-0
        counter6$countervalue6<-0
        counter7$countervalue7<-0
        counter8$countervalue8<-0
        shinyjs::reset('ID_col')
        shinyjs::reset('Pval_col')
        shinyjs::reset('bars')
        shinyjs::reset('ncols')
        shinyjs::reset('ID_col_position')
        shinyjs::disable('dt_confirmation')
        rv$final_graph<-NULL

    })


    observeEvent(input$dt_confirmation, {
        if (is.null(rv$data)) {
            rv$data_to_plot<-NULL
        }
            else{
                if (!(is.null(rv$data)) & length(colnames(rv$data))>=2) {
                    rv$data_to_plot=rv$data
                    if (!(TRUE %in% duplicated(colnames(rv$data_to_plot)))) {
                        rv$data_to_plot=rv$data
                    }
                        else{
                            if (TRUE %in%
                                duplicated(colnames(rv$data_to_plot))) {
                                colnames(rv$data_to_plot)<-make.unique(
                                    colnames(rv$data_to_plot))
                            }
                        }
                }
            }
    })


    output$Tab_3_main_panel_progress_info<-renderUI({
        if ((is.null(rv$data_to_plot)) & is.null(rv$data) & is.null(input$file)
            | (is.null(rv$data_to_plot)) & is.null(rv$data) &
            !(is.null(input$file))) {
            shinyjs::disable('bars')
            shinyjs::disable('ncols')
            shinyjs::disable('Confirm_new_name')
            shinyjs::disable('Reorder_columns')
            shinyjs::disable('Formating_end_confirmation')
            shinyjs::disable('upval')
            shinyjs::disable('Plot_Height')
            shinyjs::disable('Plot_Width')
            shinyjs::disable('Plot_Graph')
            shinyjs::disable('Download_Plot_Buttom')
            h5('First you need to upload a file in the \'Data Input\' tab,
               select the desired dataframe, create it and set it as ready
               for formating in the \'Data selection tab\'. Only then you
               can format the dataframe!')
        }
            else{
                if ((is.null(rv$data_to_plot)) & !(is.null(rv$data)) &
                    !(is.null(input$file)) & length(colnames(rv$data))<2) {
                    shinyjs::disable('bars')
                    shinyjs::disable('ncols')
                    shinyjs::disable('Confirm_new_name')
                    shinyjs::disable('Reorder_columns')
                    shinyjs::disable('Formating_end_confirmation')
                    shinyjs::disable('upval')
                    shinyjs::disable('Plot_Height')
                    shinyjs::disable('Plot_Width')
                    shinyjs::disable('Plot_Graph')
                    shinyjs::disable('Download_Plot_Buttom')
                    h5('To create a dataframe you must select at least one
                       column besides the id column.')
                }
                    else{
                        if ((is.null(rv$data_to_plot)) & !(is.null(rv$data)) &
                            !(is.null(input$file))) {
                            shinyjs::disable('bars')
                            shinyjs::disable('ncols')
                            shinyjs::disable('Confirm_new_name')
                            shinyjs::disable('Reorder_columns')
                            shinyjs::disable('Formating_end_confirmation')
                            shinyjs::disable('upval')
                            shinyjs::disable('Plot_Height')
                            shinyjs::disable('Plot_Width')
                            shinyjs::disable('Plot_Graph')
                            shinyjs::disable('Download_Plot_Buttom')
                            h5('You have selected and sucessfully created a
                            dataframe but it is not ready for formating yet.
                            Go in the \'Data Selection\' tab and hit the
                             \'Dataframe is ready formating!\' buttom them come
                            back to next tab to format the dataframe.')
                        }


                        else{
                            if (!(is.null(rv$data_to_plot)) &
                                !(is.null(rv$data))
                                & !(is.null(input$file))) {
                                shinyjs::enable('bars')
                                shinyjs::enable('ncols')
                                shinyjs::enable('upval')
                                h5('Your dataframe is ready for formating!
                                Now format your dataframe and hit the
                                \'Dataframe is ready for plotting!\' buttom.')
                            }
                        }
                    }
            }
    })


    counter4<-reactiveValues(countervalue4=1)

    observeEvent(input$upval,{
        if(counter3$countervalue3==1){
            counter4$countervalue4 <- 1
        }
    })


    bnv<-reactiveValues()
    counter5<-reactiveValues(countervalue5=1)


    observeEvent(input$upval,{
        if ((input$bars*input$ncols)==(length(colnames(rv$data_to_plot))-1)) {
            if(counter4$countervalue4==1){
                counter5$countervalue5 <- 1
                bnv$bnames<-paste0('bar_',seq(1:input$bars))
            }
        }
            else{
                if ((input$bars*input$ncols)!=
                    (length(colnames(rv$data_to_plot))-1)) {
                    warn2 <<- showNotification(
                        'The number of bars and columns per bar you choose
                         must meet the following condition:\n
                         (bars)*(cols per bar) =
                         (number of columns in your created dataframe)-1',
                         duration = 10,
                         closeButton = TRUE,
                         type = 'warning'
                    )
                    counter5$countervalue5<-counter5$countervalue5
                }
            }
    })



    output$Column_rename_selection<-renderUI({
        if (!(is.null(rv$data_to_plot)) & counter6$countervalue6==1) {
            shinyjs::enable('Confirm_new_name')
            selectInput('to_rename_col','Select the column to be renamed:',
                choices = colnames(rv$data_to_plot),multiple = FALSE)
        }
            else{
                if (!(is.null(rv$data_to_plot)) & counter6$countervalue6==0 ) {
                    shinyjs::disabled(selectInput(
                        'to_rename_col','Select the column to be renamed:',
                        choices = c(),multiple = FALSE))
                }
                    else{
                        if (((is.null(rv$data_to_plot))) ) {
                            shinyjs::disabled(selectInput('to_rename_col',
                                'Select the column to be renamed:',
                                choices = c(),multiple = FALSE))
                        }
                    }
            }
    })


    output$Column_rename_write<-renderUI({
        if (!(is.null(rv$data_to_plot)) & counter6$countervalue6==1) {
            textInput('col_new_name','Write a name for the selected column:')
        }
            else{
                if (!(is.null(rv$data_to_plot)) & counter6$countervalue6==0) {
                    shinyjs::disabled(textInput('col_new_name',
                        'Write a name for the selected column:'))
                }
                    else{
                        if ((is.null(rv$data_to_plot))) {
                            shinyjs::disabled(textInput('col_new_name',
                                'Write a name for the selected column:'))
                        }
                    }
            }
    })


    observeEvent(input$Confirm_new_name,{
        if (!(is.null(input$col_new_name)) & (input$col_new_name!="")) {
            colnames(rv$data_to_plot)[(grep(paste0('^',input$to_rename_col,'$'),
                colnames(rv$data_to_plot)))]<-input$col_new_name
        }

    })

    name<-reactive({
        input$col_new_name
    })

    output$wname<-renderText({
       name()==""
    })

    counter6<-reactiveValues(countervalue6=0)


    observeEvent(input$Reorder_columns,{
        if(counter5$countervalue5==1){
            shinyjs::enable('Formating_end_confirmation')
            counter6$countervalue6 <- 1
        }
    })


    output$Column_order<-renderUI({
        if (!(is.null(rv$data_to_plot))) {
            shinyjs::enable('Reorder_columns')
            selectInput('colorder',
                'Select the order of the columns:',
                choices = colnames(rv$data_to_plot),
                multiple = TRUE)
        }else{
            shinyjs::disable('Reorder_columns')
            shinyjs::disabled(selectInput('colorder',
                        'Select the order of the columns:',
                        choices = c(rv$data_to_plot),
                        multiple = TRUE))
        }

    })


    output$Formated_dt<-renderTable({
        if (!(is.null(rv$data)) & !(is.null(rv$data_to_plot))) {
            rv$data_to_plot
        }
    })


    output$ID_column_position<-renderUI({
        if (!(is.null(rv$data_to_plot)) & counter6$countervalue6==1) {
            selectInput('ID_col_position',
                'Select the position of the id column:',
                choices = seq(1:length(colnames(rv$data_to_plot))),
                multiple = FALSE,
                selected = 1)
        }
            else{
                if (!(is.null(rv$data_to_plot)) & counter6$countervalue6==0) {
                    shinyjs::disabled(selectInput('ID_col_position',
                    'Select the position of the id column:',
                    choices = seq(1:length(colnames(rv$data_to_plot))),
                    multiple = FALSE,
                    selected = 1))
                }
                    else{
                        if ((is.null(rv$data_to_plot))) {
                            shinyjs::disabled(selectInput('ID_col_position',
                                'Select the position of the id column:',
                                choices = seq(1:length(colnames(
                                    rv$data_to_plot))),
                                multiple = FALSE,
                                selected = 1))
                        }
                    }
            }
    })


    counter7<-reactiveValues(countervalue7=0)


    observeEvent(input$Reorder_columns,{
        if (!(is.null(rv$data_to_plot)) & !(is.null(rv$data)) &
            !(is.null(input$file)) &
            (length(input$colorder)==length(colnames(rv$data_to_plot)))
            & counter3$countervalue3==1) {
            counter7$countervalue7 <- 1
            col_indexes_vector<-'c('
            for (i in seq(length(input$colorder))) {
                if(i<length(input$colorder)){
                    col_indexes_vector=paste0(col_indexes_vector,
                            grep(paste0('^',
                                input$colorder[i],'$'),
                            colnames(rv$data_to_plot)),',')
                }
                    else{
                        if( i==length(input$colorder)){
                            col_indexes_vector=paste0(col_indexes_vector,
                                grep(paste0('^',
                                    input$colorder[i],'$'),
                                    colnames(rv$data_to_plot)),')')
                        }
                    }
            }
            rv$data_to_plot<-rv$data_to_plot[eval(parse(
                text = col_indexes_vector))]
        }
    })


    counter8<-reactiveValues(countervalue8=0)


    observeEvent(input$Formating_end_confirmation,{
        if(counter6$countervalue6==1 & !(is.null(rv$data_to_plot))){
            shinyjs::enable('Plot_Graph')
            shinyjs::enable('Plot_Height')
            shinyjs::enable('Plot_Width')
            #shinyjs::enable('Download_Plot_Buttom')
            counter8$countervalue8 <- 1
        }

        rv$formated_data<-rv$data_to_plot
        temp_dt<-rv$formated_data

        temp_dt<-melt(temp_dt,id.vars =
                colnames(rv$data_to_plot)[as.numeric(input$ID_col_position)])

        temp_dt<-temp_dt%>%
            mutate(G_pvalue=ifelse(temp_dt$value<=1.0 & temp_dt$value>0.8,
                '0.8 < p.value <= 1.0',
                ifelse(temp_dt$value<=0.8 & temp_dt$value>0.5,
                '0.5 < p.value <= 0.8',
                ifelse(temp_dt$value<=0.5 & temp_dt$value>0.1,
                '0.1 < p.value <= 0.5',
                ifelse(temp_dt$value<=0.1 & temp_dt$value>0.05,
                '0.05 < p.value <= 0.1',
                ifelse(temp_dt$value<=0.05 & temp_dt$value>0.01,
                '0.01 < p.value <= 0.05',
                ifelse(temp_dt$value<=0.01 & temp_dt$value>0.001,
                '0.001 < p.value <= 0.01',
                ifelse(temp_dt$value<=0.001 & temp_dt$value>0.0,
                '0.0 < p.value <= 0.001',
                ifelse(temp_dt$value==0.0,
                'p.value = 0.0',
                'Ns')))))))))

        rv$formated_data<-temp_dt
    })


    output$Plot_element_font<-renderUI({
        if (counter8$countervalue8==1 & !(is.null(rv$data_to_plot)) ) {
            selectInput(inputId = 'Plot_elem_font',
                'Select an .eps font type:',
                choices = names(postscriptFonts()),
                selected = "Bookman")
        }
            else{
                if (counter8$countervalue8==0 & !(is.null(rv$data_to_plot))) {
                    shinyjs::disabled(selectInput(inputId = 'Plot_elem_font',
                        'Select an .eps font type:',
                        choices = names(postscriptFonts()),
                        selected = "Bookman"))
                }
                    else{
                        if ((is.null(rv$data_to_plot))) {
                        shinyjs::disabled(selectInput(inputId =
                            'Plot_elem_font',
                            'Select an .eps font type:',
                            choices = names(postscriptFonts()),
                            selected = "Bookman"))
                        }
                    }
            }
    })


    output$Plot_element_size<-renderUI({
        if (counter8$countervalue8==1 & !(is.null(rv$data_to_plot)) ) {
            numericInput(inputId = 'Plot_elem_size',
                label = 'Enter with a font size for plot elements:',
                min = 1,
                max = 20,
                value = 10 )
        }
            else{
                if (counter8$countervalue8==0 & !(is.null(rv$data_to_plot))) {
                    shinyjs::disabled(numericInput(inputId = 'Plot_elem_size',
                        label = 'Enter with a font size for plot elements:',
                        min = 1,
                        max = 20,
                        value = 10 ))
                }
                    else{
                        if ((is.null(rv$data_to_plot))) {
                            shinyjs::disabled(numericInput(
                                inputId = 'Plot_elem_size',
                                label = 'Font size for plot elements:',
                                min = 1,
                                max = 20,
                                value = 10 ))
                        }
                    }
            }
    })


    output$Plot_title_size<-renderUI({
        if (counter8$countervalue8==1 & !(is.null(rv$data_to_plot)) ) {
            numericInput(inputId = 'Plot_title_font_size',
                label = 'Enter with a font size for the plot title:',
                min = 1,
                max = 20,
                value = 10)
        }
            else{
                if (counter8$countervalue8==0 & !(is.null(rv$data_to_plot))) {
                    shinyjs::disabled(
                        numericInput(inputId = 'Plot_title_font_size',
                        label = 'Enter with a font size for the plot title:',
                        min = 1,max = 20,value = 10))
                }
                    else{
                        if ((is.null(rv$data_to_plot))) {
                            shinyjs::disabled(
                                numericInput(inputId = 'Plot_title_font_size',
                                label = 'Font size for the plot title:',
                                min = 1,
                                max = 20,
                                value = 10))
                        }
                    }
            }
    })


    output$Plot_legend_size<-renderUI({
        if (counter8$countervalue8==1 & !(is.null(rv$data_to_plot)) ) {
            numericInput(inputId = 'Plot_legend_font_size',
            label = 'Enter with a font size for the legend:',
            min = 1,
            max = 20,
            value = 10)
        }
            else{
                if (counter8$countervalue8==0 & !(is.null(rv$data_to_plot))) {
                    shinyjs::disabled(
                        numericInput(inputId = 'Plot_legend_font_size',
                        label = 'Enter with a font size for the legend:',
                        min = 1,max = 20,value = 10))
                }
                    else{
                        if ((is.null(rv$data_to_plot))) {
                            shinyjs::disabled(
                                numericInput(inputId = 'Plot_legend_font_size',
                                label = 'Font size for the legend:',
                                min = 1,
                                max = 20,
                                value = 10))
                        }
                    }
            }

  })


    output$Plot_axis_text_face<-renderUI({
        if (counter8$countervalue8==1 & !(is.null(rv$data_to_plot)) ) {
            selectInput(inputId = 'Plot_axtxt_face',
                label = 'Font face for axis text:',
                choices = c('bold','italic','bold.italic','plain'),
                selected = 'plain')
        }
            else{
                if (counter8$countervalue8==0 & !(is.null(rv$data_to_plot))) {
                    shinyjs::disabled(
                        selectInput(inputId = 'Plot_axtxt_face',
                            label = 'Font face for axis text:',
                            choices = c('bold','italic','bold.italic','plain'),
                            selected = 'plain'))
                }
                    else{
                        if ((is.null(rv$data_to_plot))) {
                            shinyjs::disabled(
                                selectInput(inputId = 'Plot_axtxt_face',
                                    label = 'Font face for axis text:',
                                    choices = c('bold','italic','bold.italic',
                                                'plain'),
                                    selected = 'plain'))
                        }
                    }
            }
    })

    output$Plot_axis_title_face<-renderUI({
        if (counter8$countervalue8==1 & !(is.null(rv$data_to_plot)) ) {
            selectInput(inputId = 'Plot_axtitle_face',
                label = 'Select a font face for axis title:',
                choices = c('bold','italic','bold.italic','plain'),
                selected = 'plain')
        }
            else{
                if (counter8$countervalue8==0 & !(is.null(rv$data_to_plot))) {
                    shinyjs::disabled(selectInput(inputId =
                        'Plot_axtitle_face',
                        label = 'Font face for axis title:',
                        choices = c('bold','italic','bold.italic','plain'),
                        selected = 'plain'))
                }
                    else{
                        if ((is.null(rv$data_to_plot))) {
                            shinyjs::disabled(selectInput(
                                inputId = 'Plot_axtitle_face',
                                label = 'Font face for axis title:',
                                choices = c('bold','italic','bold.italic',
                                            'plain'),
                                selected = 'plain'))
                        }
                    }
            }
    })

    output$Plot_legend_face<-renderUI({
        if (counter8$countervalue8==1 & !(is.null(rv$data_to_plot)) ) {
            selectInput(inputId = 'Plot_lgd_face',
                label = 'Select a font face for the legend:',
                choices = c('bold','italic','bold.italic','plain'),
                selected = 'plain')
        }
            else{
                if (counter8$countervalue8==0 & !(is.null(rv$data_to_plot))) {
                    shinyjs::disabled(selectInput(
                        inputId = 'Plot_lgd_face',
                        label = 'Select a font face for the legend:',
                        choices = c('bold','italic','bold.italic','plain'),
                        selected = 'plain'))
                }
                    else{
                        if ((is.null(rv$data_to_plot))) {
                            shinyjs::disabled(selectInput(
                                inputId = 'Plot_lgd_face',
                                label = 'Select a font face for the legend:',
                                choices = c('bold','italic','bold.italic',
                                            'plain'),
                                selected = 'plain'))
                        }
                    }
            }
    })


    output$Plot_obj <- renderUI({
        makePlotContainers(as.numeric(input$nplots),
            ncol=as.numeric(input$ncol),
            height=input$Plot_Height,
            width=input$Plot_Width)
    })


    output$text<-renderText({
        as.character(rv$brk)
    })


    observeEvent(input$Plot_Graph,{
        an<--(grep(colnames(rv$data_to_plot)
            [as.numeric(input$ID_col_position)],colnames(rv$data_to_plot)))

        cnames<-colnames(rv$data_to_plot)[an]

        cnt<-0
        rv$brk<-list()
        rv$labs<-list()

        for (i in seq(1:as.numeric(1))) {
            aa<-'c('
            bb<-'c('
            for (j in seq(1:as.numeric(length(colnames(rv$data_to_plot))-1))) {
                cnt<-cnt+1
                if ((length(colnames(rv$data_to_plot))-1)>1) {
                    if(j==1 |j < ((length(
                        seq(1:as.numeric(length(colnames(rv$data_to_plot)
                        )))))-1)){
                        aa<-paste0(aa,'"',cnames[cnt],'"',',')
                        bb<-paste0(bb,'"',cnames[cnt],'"','=','"',
                                   cnames[cnt],'"',',')
                    }
                    else{
                        if(j == ((length(
                            seq(1:as.numeric(length(colnames(rv$data_to_plot)
                            )))))-1)){
                            aa<-paste0(aa,'"',cnames[cnt],'"',')')
                            bb<-paste0(bb,'"',cnames[cnt],'"','=','"',
                                       cnames[cnt],'"',')')
                        }
                    }
                }
                else{
                    aa<-paste0(aa,'"',cnames[cnt],'"',')')
                    bb<-paste0(bb,'"',cnames[cnt],'"','=','"',
                               cnames[cnt],'"',')')
                }

            }
            rv$brk[[i]]<-aa
            rv$labs[[i]]<-bb
        }
        gnames<-c()
        rv$glist<-list()
        rv$glist2<-list()
        for (k in seq(1)) {
            assign(paste0('brkv'),eval(parse(text = rv$brk[[k]] )))
            assign(paste0('gg',k),eval(parse(text = paste0(
                "ggplot(rv$formated_data,aes_string(x=colnames(
                rv$formated_data)[2],y=colnames(rv$formated_data)
                [as.numeric(input$ID_col_position)])) +
                geom_tile(aes(fill=G_pvalue),
                colour='white', size=0.25) +
                scale_x_discrete(expand = c(0,0),position = 'top',breaks=",
                'brkv',", limits=brkv,labels=brkv) +
                scale_y_discrete(expand = c(0,0)) + coord_equal(ratio = 1) +
                scale_fill_manual(values =
                c(cols<- c('0.8 < p.value <= 1.0'='#e7f0fa',
                #lighter than light blue
                '0.5 < p.value <= 0.8'='#c9e2f6', #light blue
                '0.1 < p.value <= 0.5'='#95cbee', #blue
                '0.05 < p.value <= 0.1'='#0099dc', #darker blue
                '0.01 < p.value <= 0.05'='#4ab04a', #mustard
                '0.001 < p.value <= 0.01'='#e29421', #dark khaki
                '0.0 < p.value <= 0.001'='#f05336', #orange red
                'p.value = 0.0'='#ce472e')),#red))
                limits=c('0.8 < p.value <= 1.0',
                '0.5 < p.value <= 0.8',
                '0.1 < p.value <= 0.5',
                '0.05 < p.value <= 0.1',
                '0.01 < p.value <= 0.05',
                '0.001 < p.value <= 0.01',
                '0.0 < p.value <= 0.001',
                'p.value = 0.0')) +
                theme(axis.title.y = element_blank(),
                    axis.line.y = element_blank(),
                    axis.text = element_text(
                    family=input$Plot_elem_font,
                    size = input$Plot_elem_size,
                    face = input$Plot_axtxt_face),
                    axis.title = element_text(
                    family=input$Plot_elem_font,
                    size = input$Plot_title_font_size,
                    face = input$Plot_axtitle_face),
                    legend.text=element_text(
                    family=input$Plot_elem_font,
                    size = input$Plot_legend_font_size,
                    face = input$Plot_lgd_face))  +
                    labs(x='') +
                    guides(fill=guide_legend(title = 'p-value'))"))))


            gnames<-append(gnames,paste0('gg',k),length(gnames))


            rv$glist[[k]]<-eval(parse(text = paste0(
                "ggplot(rv$formated_data,aes_string(
                x=colnames(rv$formated_data)[2],
                y=(colnames(rv$formated_data)
                [as.numeric(input$ID_col_position)]))) +
                geom_tile(aes(fill=G_pvalue),colour='white', size=0.25) +
                scale_x_discrete(expand = c(0,0),position = 'top',breaks=",
                'brkv',", limits=brkv,labels=brkv) +
                scale_y_discrete(expand = c(0,0)) +
                coord_equal(ratio = 1) +
                scale_fill_manual(values =
                c(cols<- c('0.8 < p.value <= 1.0'='#e7f0fa',
                #lighter than light blue
                '0.5 < p.value <= 0.8'='#c9e2f6', #light blue
                '0.1 < p.value <= 0.5'='#95cbee', #blue
                '0.05 < p.value <= 0.1'='#0099dc', #darker blue
                '0.01 < p.value <= 0.05'='#4ab04a', #mustard
                '0.001 < p.value <= 0.01'='#e29421', #dark khaki
                '0.0 < p.value <= 0.001'='#f05336', #orange red
                'p.value = 0.0'='#ce472e')),#red))
                limits=c('0.8 < p.value <= 1.0',
                '0.5 < p.value <= 0.8',
                '0.1 < p.value <= 0.5',
                '0.05 < p.value <= 0.1',
                '0.01 < p.value <= 0.05',
                '0.001 < p.value <= 0.01',
                '0.0 < p.value <= 0.001',
                'p.value = 0.0')) +
                theme(axis.title.y = element_blank(),
                    axis.line.y = element_blank(),
                    axis.text = element_text(
                        family=input$Plot_elem_font,
                        size = input$Plot_elem_size,
                        face = input$Plot_axtxt_face),
                    axis.title = element_text(
                        family=input$Plot_elem_font,
                        size = input$Plot_title_font_size,
                        face = input$Plot_axtitle_face),
                    legend.text=element_text(
                        family=input$Plot_elem_font,
                        size = input$Plot_legend_font_size,
                        face = input$Plot_lgd_face))  +
                    labs(x='') +
                    guides(fill=guide_legend(title = 'p-value'))")))

            rv$glist2[[k]]<-eval(parse(text = paste0(
                "ggplot(rv$formated_data,aes_string(
                    x=colnames(rv$formated_data)[2],
                    y=(colnames(rv$formated_data)
                    [as.numeric(input$ID_col_position)]))) +
                    geom_tile(aes(fill=G_pvalue),colour='white', size=0.25) +
                    scale_x_discrete(expand = c(0,0),position = 'top',breaks=",
                    'brkv',", limits=brkv,labels=brkv) +
                    scale_y_discrete(expand = c(0,0)) +
                    coord_equal(ratio = 1) +
                    scale_fill_manual(values = c(cols<-
                        c('0.8 < p.value <= 1.0'='#e7f0fa',
                            #lighter than light blue
                            '0.5 < p.value <= 0.8'='#c9e2f6', #light blue
                            '0.1 < p.value <= 0.5'='#95cbee', #blue
                            '0.05 < p.value <= 0.1'='#0099dc', #darker blue
                            '0.01 < p.value <= 0.05'='#4ab04a', #mustard
                            '0.001 < p.value <= 0.01'='#e29421', #dark khaki
                            '0.0 < p.value <= 0.001'='#f05336', #orange red
                            'p.value = 0.0'='#ce472e')),#red))
                            limits=c('0.8 < p.value <= 1.0',
                            '0.5 < p.value <= 0.8',
                            '0.1 < p.value <= 0.5',
                            '0.05 < p.value <= 0.1',
                            '0.01 < p.value <= 0.05',
                            '0.001 < p.value <= 0.01',
                            '0.0 < p.value <= 0.001',
                            'p.value = 0.0')) +
                            theme(axis.title.y = element_blank(),
                                axis.line.y = element_blank(),
                                axis.text = element_text(
                                family=input$Plot_elem_font,
                                size = input$Plot_elem_size,
                                face = input$Plot_axtxt_face),
                                axis.title = element_text(
                                family=input$Plot_elem_font,
                                size = input$Plot_title_font_size,
                                face = input$Plot_axtitle_face),
                                legend.text=element_text(
                                family=input$Plot_elem_font,
                                size = input$Plot_legend_font_size,
                                face = input$Plot_lgd_face)) +
                                labs(x='') +
                                guides(
                                     fill=guide_legend(title = 'p-value'))")))

        }
        rv$grnames<-gnames
        rv$legend <- get_legend(gg1)
        legend<-rv$legend
        final_graph<-'ggarrange('
        for (l in seq(1:(as.numeric(1)))) {
            if( l==as.numeric(1)){
                final_graph<-paste0(final_graph,rv$grnames[l],
                    '+theme(legend.position =
                    "bottom",legend.text=element_text(
                    family=input$Plot_elem_font,
                    size = input$Plot_legend_font_size,
                    face = input$Plot_lgd_face)),
                    legend="bottom",common.legend=TRUE)')
            }
        }
        rv$final_graph_text<-final_graph
        rv$final_graph<-eval(parse(text = final_graph))
        output[['myplot']]<-renderPlot({
        rv$final_graph
        })
    if (class(try(eval(parse(text=rv$final_graph_text))))!='try-error') {
        shinyjs::enable('Download_Plot_Buttom')
    }
        else{
            shinyjs::disable('Download_Plot_Buttom')
        }
    })


    final_plot<-reactive({
        rv$final_graph
    })


    output$Plot_obj<-renderUI({
        plotOutput('myplot',height = validateCssUnit(input$Plot_Height),
            width = validateCssUnit(input$Plot_Width))
    })


    output$Download_Plot_Buttom <- downloadHandler(
        filename = paste0('Plot_',Sys.Date(),'.eps'),
        content = function(file) {
            postscript(file,family = input$Plot_elem_font)
            print(final_plot())
            dev.off()
        })
}

shinyApp(ui, server)
