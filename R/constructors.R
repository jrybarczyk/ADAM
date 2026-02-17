############
#Arguments checking and object building
############

.adam_validate_required_inputs <- function(
    ExpressionData,
    ComparisonID,
    AnalysisDomain,
    DBSpecies,
    GeneIdentifier
) {
    if (is.null(ExpressionData)) {
        stop("Please provide a valid expression data input.")
    }
    if (is.null(ComparisonID)) {
        stop("Please provide a valid comparison ID.")
    }
    if (is.null(AnalysisDomain)) {
        stop("Please provide a valid analysis domain.")
    }
    if (is.null(DBSpecies)) {
        stop("Please provide a valid species database.")
    }
    if (is.null(GeneIdentifier)) {
        stop("Please provide a valid gene identifier.")
    }
}

.adam_resolve_gene_limits <- function(MinGene, MaxGene) {
    if (!is.null(MinGene) && !is.null(MaxGene)) {
        GeneNumbers <- .checkGeneNumbers(MinGene, MaxGene)
        return(list(MinGene = GeneNumbers[1], MaxGene = GeneNumbers[2]))
    }
    if (is.null(MinGene) && !is.null(MaxGene)) {
        GeneNumbers <- .checkGeneNumbers(3L, MaxGene)
        return(list(MinGene = GeneNumbers[1], MaxGene = GeneNumbers[2]))
    }
    if (!is.null(MinGene) && is.null(MaxGene)) {
        GeneNumbers <- .checkGeneNumbers(MinGene, 2000L)
        return(list(MinGene = GeneNumbers[1], MaxGene = GeneNumbers[2]))
    }
    list(
        MinGene = 3L,
        MaxGene = 2000L
    )
}

.adam_resolve_analysis_inputs <- function(
    AnalysisDomain,
    DBSpecies,
    GeneIdentifier,
    GeneExpressionData
) {
    Analysis <- .checkAnalysisDomain_DBSpecies_GeneIdentifier(
        AnalysisDomain,
        DBSpecies,
        GeneIdentifier,
        GeneExpressionData
    )
    list(
        DomainGroup = Analysis[[1]],
        DataSpeciesFunctionsSample = Analysis[[2]],
        DataSpeciesFunctionsRaw = Analysis[[3]],
        GeneNomenclature = Analysis[[4]]
    )
}

.adam_resolve_complete_options <- function(
    SeedNumber,
    BootstrapNumber,
    PCorrection,
    PCorrectionMethod,
    WilcoxonTest,
    FisherTest
) {
    list(
        SeedNumber = if (is.null(SeedNumber)) 1049 else .checkSeedNumber(SeedNumber),
        BootstrapNumber = if (is.null(BootstrapNumber)) 1000L else .checkBootstrapNumber(BootstrapNumber),
        PCorrection = if (is.null(PCorrection)) 0.05 else .checkPCorrection(PCorrection),
        PCorrectionMethod = if (is.null(PCorrectionMethod)) "fdr" else .checkPCorrectionMethod(PCorrectionMethod),
        WilcoxonTest = if (is.null(WilcoxonTest)) FALSE else .checkWilcoxonTest(WilcoxonTest),
        FisherTest = if (is.null(FisherTest)) FALSE else .checkFisherTest(FisherTest)
    )
}

.adam_build_ecg_object <- function(
    ComparisonID,
    GeneExpressionData,
    GeneLimits,
    AnalysisData,
    completeTest,
    CompleteOptions = NULL
) {
    base_slots <- list(
        ComparisonID = ComparisonID,
        ExpressionData = GeneExpressionData,
        MinGene = GeneLimits$MinGene,
        MaxGene = GeneLimits$MaxGene,
        DBSpeciesFunctionsSample = AnalysisData$DataSpeciesFunctionsSample,
        DBSpeciesFunctionsRaw = AnalysisData$DataSpeciesFunctionsRaw,
        AnalysisDomain = AnalysisData$DomainGroup,
        GeneIdentifier = AnalysisData$GeneNomenclature
    )

    if (!completeTest) {
        return(do.call(new, c(list(Class = "ECGMainData"), base_slots)))
    }

    complete_slots <- c(
        base_slots,
        list(
            SeedNumber = CompleteOptions$SeedNumber,
            BootstrapNumber = CompleteOptions$BootstrapNumber,
            PCorrection = CompleteOptions$PCorrection,
            PCorrectionMethod = CompleteOptions$PCorrectionMethod,
            WilcoxonTest = CompleteOptions$WilcoxonTest,
            FisherTest = CompleteOptions$FisherTest
        )
    )
    do.call(new, c(list(Class = "ECGMainData"), complete_slots))
}

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
    .adam_validate_required_inputs(
        ExpressionData = ExpressionData,
        ComparisonID = ComparisonID,
        AnalysisDomain = AnalysisDomain,
        DBSpecies = DBSpecies,
        GeneIdentifier = GeneIdentifier
    )

    GeneExpressionData <- .checkExpressionData(ExpressionData)
    ComparisonID <- .checkComparisonID(ComparisonID, GeneExpressionData)
    GeneLimits <- .adam_resolve_gene_limits(MinGene, MaxGene)
    AnalysisData <- .adam_resolve_analysis_inputs(
        AnalysisDomain = AnalysisDomain,
        DBSpecies = DBSpecies,
        GeneIdentifier = GeneIdentifier,
        GeneExpressionData = GeneExpressionData
    )

    if (!completeTest) {
        return(.adam_build_ecg_object(
            ComparisonID = ComparisonID,
            GeneExpressionData = GeneExpressionData,
            GeneLimits = GeneLimits,
            AnalysisData = AnalysisData,
            completeTest = FALSE
        ))
    }

    CompleteOptions <- .adam_resolve_complete_options(
        SeedNumber = SeedNumber,
        BootstrapNumber = BootstrapNumber,
        PCorrection = PCorrection,
        PCorrectionMethod = PCorrectionMethod,
        WilcoxonTest = WilcoxonTest,
        FisherTest = FisherTest
    )

    .adam_build_ecg_object(
        ComparisonID = ComparisonID,
        GeneExpressionData = GeneExpressionData,
        GeneLimits = GeneLimits,
        AnalysisData = AnalysisData,
        completeTest = TRUE,
        CompleteOptions = CompleteOptions
    )
}
