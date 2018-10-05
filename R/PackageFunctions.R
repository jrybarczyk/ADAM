########
## Function to run diversity and activity analysis
########
makeAnalysis <- function(x,ECGObject,completeTest){
    ###
    #Building the comparison dataframe from expression data
    ###
    compID <- unlist(strsplit(x,","))

    message(paste0("Analysis: ",compID[1],"_Vs_",compID[2]))

    Comparison <- data.frame(ECGObject@ExpressionData[,1],
                            ECGObject@ExpressionData[,grep(compID[1],
                            colnames(ECGObject@ExpressionData))],
                            ECGObject@ExpressionData[,grep(compID[2],
                            colnames(ECGObject@ExpressionData))])
    colnames(Comparison) <- c("gene",compID[1],compID[2])

    ###
    #Building group of functionally associated genes (GFAGs)
    ###
    message("Building GFAGs ...")

    pboptions(style = 1, char = ">")

    GroupComparisons <- pblapply(ECGObject@DBSpeciesFunctionsSample,
                                function(y,Comparison){
                                Comparison[Comparison[,1] %in% y,]},
                                Comparison=Comparison)
    GroupComparisons <- GroupComparisons[lapply(lapply(GroupComparisons,
                                        function(y){if(nrow(y)>=
                                        ECGObject@MinGene & 
                                        nrow(y)<=ECGObject@MaxGene){y}
                                        else{NULL}}),length)>0]
    GroupComparisons <- GroupComparisons[order(names(GroupComparisons))]

    ###
    #Calculating gene activity and diversity
    ###
    message("Calculating activity and diversity ...")

    GroupStatistics <- pblapply(GroupComparisons, FUN = SampleStatistics,
                                Control = compID[1], Experiment = compID[2])

    if(!completeTest){
        GroupStatistics <- as.data.frame(do.call(rbind, lapply(GroupStatistics,
                                        function(y){as.vector(unlist(y))})))
        GroupStatistics$ID <- sapply(rownames(GroupStatistics), function(y) 
                            trimws(unlist(strsplit(y,"<==>"))[1]))
        GroupStatistics$Description <- sapply(rownames(GroupStatistics),
                                        function(y) trimws(unlist(
                                        strsplit(y,"<==>")))[2])
        GroupStatistics$Raw_Number_Genes <- sapply(GroupStatistics$ID,
                                            function(y) 
                                            ECGObject@DBSpeciesFunctionsRaw[[
                                            grep(y,names(
                                            ECGObject@DBSpeciesFunctionsRaw))]])
        rownames(GroupStatistics) <- NULL
        GroupStatistics <- GroupStatistics[,c(8,9,10,1:7)]
        colnames(GroupStatistics) <- c("ID","Description","Raw_Number_Genes",
                                        names(GroupStatistics[[1]]))
        ResultAnalysis <- GroupStatistics
    }else{
        ###
        #Determining the minimum and maximum number of genes for
        #each function group
        ###
        GroupSize <- sort(unique(sapply(GroupComparisons,function(y) nrow(y))))
        
        ###
        #Filtering GFAGs by size
        ###
        message("Filtering GFAGs by size ...")
        
        ResultBootstrap <- pblapply(GroupSize, function(k,GroupStatistics){
                            GroupStatistics[lapply(lapply(GroupStatistics,
                            function(y,z){if(y[1,1]==z){y}else{NULL}}, 
                            z = k),length)>0]},
                            GroupStatistics=GroupStatistics)
        
        ###
        #Running Bootstrap
        ###
        message("Running bootstrap. This may take a few minutes ...")
        
        ResultBootstrap <- pblapply(ResultBootstrap, FUN = PrepareBootstrap,
                                    ComparisonExpression = Comparison[,c(2,3)],
                                    ECGObject = ECGObject)
        
        ResultBootstrap <- as.data.frame(do.call(rbind,
                                        lapply(ResultBootstrap,function(y){
                                        do.call(rbind,lapply(y, function(z)
                                        as.vector(unlist(z))))})))
        
        ###
        #Manipulating Bootstrap output
        ###
        colnames(ResultBootstrap) <- c("Sample_Number_Genes", 
                                    paste0("H_",compID[1]),
                                    paste0("H_",compID[2]),
                                    paste0("N_",compID[1]),
                                    paste0("N_",compID[2]),"h","n",
                                    "pValue_h","pValue_n")
        ResultBootstrap$ID <- sapply(rownames(ResultBootstrap), function(y)
                                    trimws(unlist(strsplit(y,"<==>"))[1]))
        ResultBootstrap$Description <- sapply(rownames(ResultBootstrap), 
                                            function(y)
                                            trimws(unlist(strsplit(y,
                                            "<==>")))[2])
        rownames(ResultBootstrap) <- NULL
        ResultBootstrap <- ResultBootstrap[,c(length(ResultBootstrap)-1,
                                            length(ResultBootstrap),
                                            1:(length(ResultBootstrap)-2))]
        
        ###
        #Correcting p-values
        ###
        message(paste0("Correcting bootstrap p-values by ",
                    ECGObject@PCorrectionMethod," ..."))
        
        ResultBootstrap$qValue_h <- p.adjust(ResultBootstrap$pValue_h,
                                            method=ECGObject@PCorrectionMethod)
        ResultBootstrap$qValue_n <- p.adjust(ResultBootstrap$pValue_n,
                                            method=ECGObject@PCorrectionMethod)
        ResultBootstrap$Significance_h <- sapply(ResultBootstrap$qValue_h,
                                                function(y)
                                                ifelse(y<=ECGObject@PCorrection,
                                                "significative",
                                                "not_significative"))
        ResultBootstrap$Significance_n <- sapply(ResultBootstrap$qValue_n,
                                                function(y)
                                                ifelse(y<=ECGObject@PCorrection,
                                                "significative",
                                                "not_significative"))
        
        ###
        #Adding the raw gene numbers per function to the result dataframe
        ###
        ResultBootstrap$Raw_Number_Genes <- sapply(ResultBootstrap$ID,
                                            function(y)
                                            ECGObject@DBSpeciesFunctionsRaw[[
                                            grep(y,names(
                                            ECGObject@DBSpeciesFunctionsRaw))]])
        
        ###
        #Changing the order of the dataframe result columns
        ###
        ResultBootstrap <- ResultBootstrap[,c(1,2,length(ResultBootstrap),
                                            3:(length(ResultBootstrap)-1))]
        ResultBootstrap <- arrange(ResultBootstrap,ResultBootstrap[,1])
        
        ###
        #Wicoxon Test
        ###
        if(ECGObject@WilcoxonTest){
            message("Runing Wilcox rank sum test ...")
            ResultBootstrap$Wilcox_pvalue <- pbsapply(GroupComparisons,
                                            FUN = WCAnalysis)
            ResultBootstrap$Wilcox_qvalue <- p.adjust(as.numeric(
                                            ResultBootstrap$Wilcox_pvalue),
                                            method=ECGObject@PCorrectionMethod)
            ResultBootstrap$Wilcox_significance<-ifelse(
                                                ResultBootstrap$Wilcox_qvalue <=
                                                ECGObject@PCorrection,
                                                "significative",
                                                "not significative")
        }
        
        ###
        #Fisher Test
        ###
        if(ECGObject@FisherTest){
            message("Running Fisher exact test ...")
            ResultBootstrap$Fisher_pvalue <- pbsapply(GroupComparisons,
                                                    FUN=FSHAnalysis,
                                                    Comparison = Comparison)
            ResultBootstrap$Fisher_qvalue <- p.adjust(
                                            as.numeric(
                                            ResultBootstrap$Fisher_pvalue),
                                            method=ECGObject@PCorrectionMethod)
            ResultBootstrap$Fisher_significance <- ifelse(
                                                ResultBootstrap$Fisher_qvalue <= 
                                                ECGObject@PCorrection,
                                                "significative",
                                                "not significative")
        }
        
        ResultAnalysis <- ResultBootstrap
    }

    ###
    #Analysis concluded
    ###
    message(paste0("Analysis of ",compID[1],"_Vs_",compID[2]),
            " succesfully concluded!")
    message("")

    return(ResultAnalysis)
}

########
# Function to calculate the activity of gene set.
########
Activity <- function(x, ColumnDataNumber){
    N <- sum(x[,ColumnDataNumber])
    return(N)
}

########
# Function to calculate the diversity of gene set.
########
Diversity <- function(x, ColumnDataNumber){
    ColumnData <- x[,ColumnDataNumber]
    ColumnData <- ColumnData[ColumnData>0]
    if(length(ColumnData)>0 & length(x[,ColumnDataNumber])>1){
        H <- -(1/(log(length(x[,ColumnDataNumber]),
                base = exp(1))))*sum(ColumnData/(sum(ColumnData))*
                log(ColumnData/sum(ColumnData), base = exp(1)), na.rm = FALSE)
    }else{
        H<-0
    }
    return(H)
}

########
# Function to calculate activity (N), diversity (H), relative activity (n) and  
# relative diversity (h) of a gene set.
########
SampleStatistics <- function(y,Control,Experiment){
    TempDataFrame <- y[,c(2,3)]
    H_Control <- Diversity(TempDataFrame, 1)
    H_Experiment <- Diversity(TempDataFrame, 2)
    N_Control <- Activity(TempDataFrame, 1)
    N_Experiment <- Activity(TempDataFrame, 2)
    if((H_Experiment + H_Control)>0){
        h <- H_Experiment/(H_Experiment + H_Control)
    }else{
        h <- 0
    }
    if((N_Experiment + N_Control)>0){
        n <- N_Experiment/(N_Experiment + N_Control)
    }else{
        n <-0
    }
    Result <- data.frame(nrow(y), H_Control, H_Experiment, N_Control,
                        N_Experiment, h, n)
    colnames(Result) <- c("Sample_Number_Genes", paste0("H_",Control),
                        paste0("H_",Experiment),
                        paste0("N_",Control), paste0("N_",Experiment),
                        "h","n")
    return(Result)
}

########
## Function to prepare data expression function groups for bootstrapping
########
PrepareBootstrap <- function(k,ComparisonExpression,ECGObject){

    PValues <- MakeBootstrap(BootstrapData = as.matrix(ComparisonExpression), 
                            BootstrapNumber = ECGObject@BootstrapNumber,
                            BootstrapGroupSize = as.numeric(k[[1]][1]),
                            BootstrapSeed = 
                            ECGObject@SeedNumber*as.numeric(k[[1]][1]))

    ResultPValues <- k
    ResultPValues <- lapply(ResultPValues, function(z, PValues){
                    z$pValue_h <- sum(PValues[,1] > as.numeric(z[6]))/
                                    ECGObject@BootstrapNumber
                    z$pValue_n <- sum(PValues[,2] > as.numeric(z[7]))/
                                    ECGObject@BootstrapNumber
    return(z)}, PValues = PValues)

    return(ResultPValues)
}

########
## Function to run Wilcoxon test.
########
WCAnalysis <- function(x){
    GroupData <- x[,c(2,3)]
    if(length(unique(GroupData[,1][GroupData[,1]>0]))==length(GroupData[,1]) &
        length(unique(GroupData[,2][GroupData[,2]>0]))==length(GroupData[,2]) &
        length(unique(c(GroupData[,1][GroupData[,1]>0],
            GroupData[,2][GroupData[,2]>0])))==
            (length(GroupData[,2])*2) & length(GroupData[,1]>1)){
        return(wilcox.test(GroupData[,1],GroupData[,2], paired = TRUE)$p.value)
    }else{
        return(1)
    }
}

########
## Function to run the Fisher test.
########
FSHAnalysis <- function(x,Comparison){
    INGroupGenes <- mutate(x,relative = x[,3]/(x[,2]+x[,3]))
    OUTGroupGenes <- filter(Comparison,!(as.character(Comparison[,1]) %in%
                            as.character(INGroupGenes[,1])))
    OUTGroupGenes <- mutate(OUTGroupGenes,relative = OUTGroupGenes[,3]/
                            (OUTGroupGenes[,2]+OUTGroupGenes[,3]))

    return(fisher.test(matrix(c(sum(INGroupGenes$relative>0.5,na.rm = TRUE),
                                sum(OUTGroupGenes$relative>0.5,na.rm = TRUE),
                                sum(INGroupGenes$relative<0.5,na.rm = TRUE),
                                sum(OUTGroupGenes$relative<0.5,na.rm = TRUE)),
                                2,2))$p.value)
}