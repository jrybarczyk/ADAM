############
#Arguments checking and object building
############
ECGMainData <- function(ComparisonID,
                        ExpressionData,
                        MinGene,
                        MaxGene,
                        SeedNumber,
                        BootstrapNumber,
                        PCorrection,
                        DBSpecies,
                        PCorrectionMethod,
                        WilcoxonTest,
                        FisherTest,
                        AnalysisDomain,
                        GeneIdentifier,
                        completeTest){

    #Check the expression data. ExpressionData must be a data frame or a path
    #for a text file tab separated containing at least 3
    #columns. First column mandatory corresponds to gene names,
    #according to GeneIdentifier argument. Second, third and the
    #others correspond to the gene samples expression.
    if(!is.null(ExpressionData)){
        GeneExpressionData <- .checkExpressionData(ExpressionData)
    }else{
        stop("Please inform a valid expression data!")
    }

    #Check the comparison IDs. ComparisonID argument must be a vector in wich 
    #each element corresponds to 2 sample columns from
    #the expression data. The data sample columns in each element from the
    #vector are comma separated.
    
    if(!is.null(ComparisonID)){
        ComparisonID <- .checkComparisonID(ComparisonID,GeneExpressionData)
    }else{
        stop("Please inform a valid comparison ID!")
    }

    #Check the minimum and maximum number of genes per group of a GFAG 
    #(Group of Functionally Associated Genes). Both must be integer
    #positive values and different from zero. Besides, the argument MaxGene 
    #must be allways greater than MinGene.
    if(!is.null(MinGene) & !is.null(MaxGene)){
        GeneNumbers <- .checkGeneNumbers(MinGene,MaxGene)
        InfGeneLimit <- GeneNumbers[1]
        SupGeneLimit <- GeneNumbers[2]
    }
    if(is.null(MinGene)){
        InfGeneLimit <- 3
    }
    if(is.null(MaxGene)){
        SupGeneLimit <- 2000
    }

    #Check the domain analysis, species reference database and gene 
    #identifier. The argument AnalysisDomain must be a character
    #corresponding to a domain (gobp, gocc, gomf, kegg or own). The argument
    #DBSpecies must be a character corresponding to an
    #OrgDb species package (org.Hs.eg.db, org.Dm.eg.db ...) or a character 
    #path for an own gene annotation file containing 3 columns:
    #gene name, term annotation and description of the term annotation. 
    #The GeneIdentifier argument must be a character containing
    #the nomenclature to be used (symbol or entrez).
    if(!is.null(AnalysisDomain) & !is.null(DBSpecies) & 
        !is.null(GeneIdentifier)){
        Analysis<-.checkAnalysisDomain_DBSpecies_GeneIdentifier(
                    AnalysisDomain,DBSpecies,GeneIdentifier,
                    GeneExpressionData)
        DomainGroup <- Analysis[[1]]
        DataSpeciesFunctionsSample <- Analysis[[2]]
        DataSpeciesFunctionsRaw <- Analysis[[3]]
        GeneNomenclature <- Analysis[[4]]
    }
    if(is.null(AnalysisDomain)){
        stop("Please inform a valid domain!")
    }
    if(is.null(DBSpecies)){
        stop("Please inform a valid database species!")
    }
    if(is.null(GeneIdentifier)){
        stop("Please inform a valid gene identifier!")
    }

    if(completeTest){
        #Check the seed for random generation numbers. The argument SeedNumber
        #must be a numeric value allways positive, greater than or
        #equal to zero.
        if(!is.null(SeedNumber)){
            SeedNumber <- .checkSeedNumber(SeedNumber)
        }else{
            SeedNumber <- 10049
        }
        
        #Check the number of bootstraps necessary for defining GFAG p-values.
        #The argument BootstrapNumber must be an integer number
        #greater than zero.
        if(!is.null(BootstrapNumber)){
            BootstrapNumber <- .checkBootstrapNumber(BootstrapNumber)
        }else{
            BootstrapNumber <- 1000
        }
        
        #Check the cutoff to be used for one of the p-value correction 
        #methods. The PCorrection argument must be a numeric value between
        #zero and one.
        if(!is.null(PCorrection)){
            PCorrection <- .checkPCorrection(PCorrection)
        }else{
            PCorrection <- 0.05
        }
        
        #Check the p-value correction method. The PCorrectionMethod argument
        #must be a character corresponding to one of the p.adjust function
        #correction methods (holm, hochberg, hommel, bonferroni, BH, BY, fdr).
        if(!is.null(PCorrectionMethod)){
            PCorrectionMethod <- .checkPCorrectionMethod(PCorrectionMethod)
        }else{
            PCorrectionMethod <- "fdr"
        }
        
        #Check if it will be performed the Wilcoxon Rank Sum Test. The 
        #WilcoxonTest argument should be TRUE for running the test or FALSE
        #if the test won't be performed.
        if(!is.null(WilcoxonTest)){
            WilcoxonTest <- .checkWilcoxonTest(WilcoxonTest)
        }else{
            WilcoxonTest <- FALSE
        }
        
        #Check if it will be performed the Fisher Exact Test. The FisherTest
        #argument should be TRUE for running the test or FALSE
        #if the test won't be performed.
        if(!is.null(FisherTest)){
            FisherTest <- .checkFisherTest(FisherTest)
        }else{
            FisherTest <- FALSE
        }
        
        #Object building, according to the necessary arguments for running 
        #complete analysis ADAM.
        inputObject <- new(Class = "ECGMainData",
                        ComparisonID = ComparisonID,
                        ExpressionData = GeneExpressionData,
                        MinGene = InfGeneLimit,
                        MaxGene = SupGeneLimit,
                        SeedNumber = SeedNumber,
                        BootstrapNumber = BootstrapNumber,
                        PCorrection = PCorrection,
                        DBSpeciesFunctionsSample = DataSpeciesFunctionsSample,
                        DBSpeciesFunctionsRaw = DataSpeciesFunctionsRaw,
                        PCorrectionMethod = PCorrectionMethod,
                        WilcoxonTest = WilcoxonTest,
                        FisherTest = FisherTest,
                        AnalysisDomain = DomainGroup,
                        GeneIdentifier = GeneNomenclature)
    }else{
        #Object building, according to the necessary arguments for running 
        #parcial analysis with ADAM.
        inputObject <- new(Class = "ECGMainData",
                        ComparisonID = ComparisonID,
                        ExpressionData = GeneExpressionData,
                        MinGene = InfGeneLimit,
                        MaxGene = SupGeneLimit,
                        DBSpeciesFunctionsSample = DataSpeciesFunctionsSample,
                        DBSpeciesFunctionsRaw = DataSpeciesFunctionsRaw,
                        AnalysisDomain = DomainGroup,
                        GeneIdentifier = GeneNomenclature)
    }
    return(inputObject)
}
