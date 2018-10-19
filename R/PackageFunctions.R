########
## Function to run diversity and activity analysis
########

makeAnalysis <- function(ComparisonID,ECGObject,completeTest){
    ###
    #Building the comparison dataframe from expression data
    ###
    compID <- unlist(strsplit(ComparisonID,","))

    message(paste0("Analysis: ",compID[1],"_Vs_",compID[2]))

    Comparison <- data.frame(ECGObject@ExpressionData$gene,
                            ECGObject@ExpressionData[,grep(compID[1],
                            colnames(ECGObject@ExpressionData))],
                            ECGObject@ExpressionData[,grep(compID[2],
                            colnames(ECGObject@ExpressionData))])
    colnames(Comparison) <- c("gene","control","experiment")

    ###
    #Building group of functionally associated genes (GFAGs)
    ###
    message("Building GFAGs ...")

    pboptions(style = 1, char = ">")

    GroupComparisons <- pblapply(ECGObject@DBSpeciesFunctionsSample,
                                function(DBFSample,Comparison){
                                Comparison[Comparison$gene %in% DBFSample,]},
                                Comparison=Comparison)
    GroupComparisons <- GroupComparisons[lapply(lapply(GroupComparisons,
                            function(GComparison)
                            {if(nrow(GComparison)>=ECGObject@MinGene & 
                            nrow(GComparison)<=ECGObject@MaxGene){GComparison}
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
                        function(GroupStats){as.vector(unlist(GroupStats))})))
        GroupStatistics$ID <- vapply(rownames(GroupStatistics),
                        function(GroupStats) 
                        trimws(unlist(strsplit(GroupStats,"<==>"))[1]), 
                        FUN.VALUE = character(1))
        GroupStatistics$Description <- vapply(rownames(GroupStatistics),
                                        function(GroupStats) trimws(unlist(
                                        strsplit(GroupStats,"<==>")))[2], 
                                        FUN.VALUE = character(1))
        
        
        GroupStatistics$Raw_Number_Genes<-vapply(GroupStatistics$ID,
                                        function(GroupStats) 
                                        ECGObject@DBSpeciesFunctionsRaw[[
                                        grep(GroupStats,names(
                                        ECGObject@DBSpeciesFunctionsRaw))]], 
                                        FUN.VALUE = numeric(1))
        
        
        rownames(GroupStatistics) <- NULL
        GroupStatistics <- GroupStatistics[,c(8,9,10,1:7)]
        colnames(GroupStatistics) <- c("ID","Description","Raw_Number_Genes",
                                    "Sample_Number_Genes", 
                                    paste0("H_",compID[1]),
                                    paste0("H_",compID[2]),
                                    paste0("N_",compID[1]),
                                    paste0("N_",compID[2]),"h","n")

        ResultAnalysis <- GroupStatistics
    }else{
        ###
        #Determining the minimum and maximum number of genes for
        #each function group
        ###
        GroupSize <- sort(unique(vapply(GroupComparisons,
                    function(GComparison) nrow(GComparison),
                    FUN.VALUE = numeric(1))))
        
        ###
        #Filtering GFAGs by size
        ###
        message("Filtering GFAGs by size ...")
        
        ResultBootstrap <- pblapply(GroupSize, function(GSize,GroupStatistics){
                            GroupStatistics[lapply(lapply(GroupStatistics,
                            function(GStats,GSizeStats){
                            if(GStats[1,1]==GSizeStats){GStats}else{NULL}}, 
                            GSizeStats = GSize),length)>0]},
                            GroupStatistics=GroupStatistics)
        
        ###
        #Running Bootstrap
        ###
        message("Running bootstrap. This may take a few minutes ...")
        
        ComparisonExpression<-data.frame(Comparison$control,
                                         Comparison$experiment)
        colnames(ComparisonExpression) <- c("control","experiment")
        
        ResultBootstrap <- pblapply(ResultBootstrap, FUN = PrepareBootstrap,
                                    ComparisonExpression = ComparisonExpression,
                                    ECGObject = ECGObject)
        

        ResultBootstrap <- as.data.frame(do.call(rbind,
                                    lapply(ResultBootstrap,function(ResultBtp){
                                    do.call(rbind,lapply(ResultBtp,
                                        function(ResultBtpBind)
                                        as.vector(unlist(ResultBtpBind))))})))
        
        ###
        #Manipulating Bootstrap output
        ###
        colnames(ResultBootstrap) <- c("Sample_Number_Genes", 
                                    paste0("H_",compID[1]),
                                    paste0("H_",compID[2]),
                                    paste0("N_",compID[1]),
                                    paste0("N_",compID[2]),"h","n",
                                    "pValue_h","pValue_n")
        ResultBootstrap$ID <- vapply(rownames(ResultBootstrap),
                                function(ResultBtp)
                                trimws(unlist(strsplit(ResultBtp,"<==>"))[1]),
                                FUN.VALUE = character(1))
        ResultBootstrap$Description <- vapply(rownames(ResultBootstrap), 
                                            function(ResultBtp)
                                            trimws(unlist(strsplit(ResultBtp,
                                            "<==>")))[2],FUN.VALUE=character(1))
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
        ResultBootstrap$Significance_h <- vapply(ResultBootstrap$qValue_h,
                                    function(ResultQValue)
                                    ifelse(ResultQValue<=ECGObject@PCorrection,
                                        "significative","not_significative"),
                                    FUN.VALUE = character(1))
        ResultBootstrap$Significance_n <- vapply(ResultBootstrap$qValue_n,
                                    function(ResultQValue)
                                    ifelse(ResultQValue<=ECGObject@PCorrection,
                                        "significative","not_significative"), 
                                    FUN.VALUE = character(1))
        
        ###
        #Adding the raw gene numbers per function to the result dataframe
        ###
        ResultBootstrap$Raw_Number_Genes <- vapply(ResultBootstrap$ID,
                                        function(GroupID)
                                        ECGObject@DBSpeciesFunctionsRaw[[
                                        grep(GroupID,names(
                                        ECGObject@DBSpeciesFunctionsRaw))]],
                                        FUN.VALUE = numeric(1))
        
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
            ResultBootstrap$Wilcox_pvalue <- vapply(GroupComparisons,
                                            FUN = WCAnalysis, 
                                            FUN.VALUE = numeric(1))
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
            ResultBootstrap$Fisher_pvalue <- vapply(GroupComparisons,
                                                    FUN=FSHAnalysis,
                                                    Comparison = Comparison, 
                                                    FUN.VALUE = numeric(1))
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
Activity <- function(DataGroup, ColumnDataNumber){
    # N = gene activity of a GFAG
    N <- sum(DataGroup[,ColumnDataNumber])
    return(N)
}

########
# Function to calculate the diversity of gene set.
########
Diversity <- function(DataGroup, ColumnDataNumber){
    ColumnData <- DataGroup[,ColumnDataNumber]
    ColumnData <- ColumnData[ColumnData>0]
    # H = gene diversity of a GFAG
    if(length(ColumnData)>0 & length(DataGroup[,ColumnDataNumber])>1){
        H <- -(1/(log(length(DataGroup[,ColumnDataNumber]),
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
SampleStatistics <- function(GroupComparisons,Control,Experiment){
    TempDataFrame <- GroupComparisons[,c(2,3)]
    H_Control <- Diversity(TempDataFrame, 1)
    H_Experiment <- Diversity(TempDataFrame, 2)
    N_Control <- Activity(TempDataFrame, 1)
    N_Experiment <- Activity(TempDataFrame, 2)
    if((H_Experiment + H_Control)>0){
        # h = gene relative diversity between control and experiment samples
        h <- H_Experiment/(H_Experiment + H_Control)
    }else{
        h <- 0
    }
    if((N_Experiment + N_Control)>0){
        # n = gene relative activity between control and experiment samples
        n <- N_Experiment/(N_Experiment + N_Control)
    }else{
        n <-0
    }
    Result<-data.frame(nrow(GroupComparisons),H_Control,H_Experiment,N_Control,
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
PrepareBootstrap <- function(GroupFunction,ComparisonExpression,ECGObject){

    PValues <- MakeBootstrap(BootstrapData = as.matrix(ComparisonExpression), 
                        BootstrapNumber = ECGObject@BootstrapNumber,
                        BootstrapGroupSize = as.numeric(GroupFunction[[1]][1]),
                        BootstrapSeed = 
                        ECGObject@SeedNumber*as.numeric(GroupFunction[[1]][1]))

    ResultPValues <- GroupFunction
    ResultPValues <- lapply(ResultPValues, function(RPValues, PValues){
            RPValues$pValue_h <- sum(PValues[,1] > as.numeric(RPValues[6]))/
                                    ECGObject@BootstrapNumber
            RPValues$pValue_n <- sum(PValues[,2] > as.numeric(RPValues[7]))/
                                    ECGObject@BootstrapNumber
            return(RPValues)}, PValues = PValues)

    return(ResultPValues)
}

########
## Function to run Wilcoxon test.
########
WCAnalysis <- function(pValuesGroups){
    GroupData <- pValuesGroups[,c(2,3)]
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
FSHAnalysis <- function(pValuesGroups,Comparison){
    INGroupGenes <- mutate(pValuesGroups,relative = pValuesGroups[,3]/
(pValuesGroups[,2]+pValuesGroups[,3]))
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