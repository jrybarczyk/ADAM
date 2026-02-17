########
## Function to run diversity and activity analysis
########

.adam_split_comparison_id <- function(ComparisonID) {
    if (!is.character(ComparisonID) || length(ComparisonID) != 1L || is.na(ComparisonID)) {
        stop("ComparisonID must be a single non-missing character string.")
    }
    compID <- trimws(unlist(strsplit(ComparisonID, ",")))
    if (length(compID) != 2L || any(compID == "")) {
        stop("ComparisonID must contain exactly two sample names: control,experiment.")
    }
    compID
}

.adam_get_sample_column <- function(ExpressionData, SampleName) {
    match_idx <- which(colnames(ExpressionData) == SampleName)
    if (length(match_idx) == 0L) {
        stop(sprintf("Sample '%s' was not found in ExpressionData.", SampleName))
    }
    if (length(match_idx) > 1L) {
        stop(sprintf("Sample '%s' appears multiple times in ExpressionData.", SampleName))
    }
    ExpressionData[, match_idx]
}

.adam_extract_raw_id <- function(NameID) {
    parts <- trimws(unlist(strsplit(NameID, "<==>", fixed = TRUE)))
    if (length(parts) >= 1L) {
        parts[[1]]
    } else {
        trimws(NameID)
    }
}

.adam_build_raw_gene_lookup <- function(DBSpeciesFunctionsRaw) {
    raw_names <- names(DBSpeciesFunctionsRaw)
    if (is.null(raw_names)) {
        stop("DBSpeciesFunctionsRaw must be a named list.")
    }
    raw_ids <- vapply(raw_names, .adam_extract_raw_id, character(1))
    if (anyDuplicated(raw_ids)) {
        dup <- unique(raw_ids[duplicated(raw_ids)])
        stop(sprintf("Duplicated raw group IDs found: %s", paste(dup, collapse = ", ")))
    }
    names(DBSpeciesFunctionsRaw) <- raw_ids
    DBSpeciesFunctionsRaw
}

.adam_get_raw_gene_count <- function(GroupID, RawLookup) {
    if (!GroupID %in% names(RawLookup)) {
        stop(sprintf("Group ID '%s' is missing in DBSpeciesFunctionsRaw.", GroupID))
    }
    as.numeric(RawLookup[[GroupID]])
}

.adam_validate_analysis_table <- function(ResultTable, compID, completeTest) {
    base_cols <- c(
        "ID", "Description", "Raw_Number_Genes", "Sample_Number_Genes",
        paste0("H_", compID[1]), paste0("H_", compID[2]),
        paste0("N_", compID[1]), paste0("N_", compID[2]), "h", "n"
    )
    complete_cols <- c(
        "pValue_h", "pValue_n", "qValue_h", "qValue_n",
        "Significance_h", "Significance_n"
    )
    required <- if (completeTest) c(base_cols, complete_cols) else base_cols
    missing <- setdiff(required, colnames(ResultTable))
    if (length(missing) > 0L) {
        stop(sprintf("Result table is missing required columns: %s", paste(missing, collapse = ", ")))
    }
    if (nrow(ResultTable) > 0L && anyDuplicated(as.character(ResultTable$ID))) {
        stop("Result table contains duplicated IDs.")
    }
    ResultTable
}

.adam_build_comparison <- function(ComparisonID, ECGObject) {
    compID <- .adam_split_comparison_id(ComparisonID)
    Comparison <- data.frame(
        ECGObject@ExpressionData$gene,
        .adam_get_sample_column(ECGObject@ExpressionData, compID[1]),
        .adam_get_sample_column(ECGObject@ExpressionData, compID[2])
    )
    colnames(Comparison) <- c("gene", "control", "experiment")
    list(compID = compID, Comparison = Comparison)
}

.adam_build_group_comparisons <- function(Comparison, ECGObject) {
    GroupComparisons <- pblapply(
        ECGObject@DBSpeciesFunctionsSample,
        function(DBFSample, ComparisonData) {
            ComparisonData[ComparisonData$gene %in% DBFSample, ]
        },
        ComparisonData = Comparison
    )

    keep_group <- vapply(GroupComparisons, function(GComparison) {
        nrow(GComparison) >= ECGObject@MinGene &&
            nrow(GComparison) <= ECGObject@MaxGene
    }, logical(1))

    GroupComparisons <- GroupComparisons[keep_group]
    GroupComparisons[order(names(GroupComparisons))]
}

.adam_group_stats_to_dataframe <- function(GroupStatistics) {
    as.data.frame(do.call(rbind, lapply(GroupStatistics, function(GroupStats) {
        as.vector(unlist(GroupStats))
    })))
}

.adam_extract_term_fields <- function(RowNames) {
    IDs <- vapply(RowNames, function(GroupName) {
        trimws(unlist(strsplit(GroupName, "<==>"))[1])
    }, character(1))
    Descriptions <- vapply(RowNames, function(GroupName) {
        trimws(unlist(strsplit(GroupName, "<==>")))[2]
    }, character(1))
    list(ID = IDs, Description = Descriptions)
}

.adam_partial_analysis_result <- function(GroupStatistics, ECGObject, compID) {
    GroupStatistics <- .adam_group_stats_to_dataframe(GroupStatistics)
    GroupFields <- .adam_extract_term_fields(rownames(GroupStatistics))
    RawLookup <- .adam_build_raw_gene_lookup(ECGObject@DBSpeciesFunctionsRaw)
    GroupStatistics$ID <- GroupFields$ID
    GroupStatistics$Description <- GroupFields$Description
    GroupStatistics$Raw_Number_Genes <- vapply(GroupStatistics$ID, function(GroupID) {
        .adam_get_raw_gene_count(GroupID, RawLookup)
    }, numeric(1))
    rownames(GroupStatistics) <- NULL
    GroupStatistics <- GroupStatistics[, c(8, 9, 10, seq_len(7))]
    colnames(GroupStatistics) <- c(
        "ID", "Description", "Raw_Number_Genes", "Sample_Number_Genes",
        paste0("H_", compID[1]), paste0("H_", compID[2]),
        paste0("N_", compID[1]), paste0("N_", compID[2]), "h", "n"
    )
    GroupStatistics
}

.adam_prepare_bootstrap_groups <- function(GroupStatistics, GroupComparisons) {
    GroupSize <- sort(unique(vapply(GroupComparisons, nrow, numeric(1))))
    pblapply(GroupSize, function(GSize, StatsList) {
        StatsList[lapply(lapply(StatsList, function(GStats, GSizeStats) {
            if (GStats[1, 1] == GSizeStats) {
                GStats
            } else {
                NULL
            }
        }, GSizeStats = GSize), length) > 0]
    }, StatsList = GroupStatistics)
}

.adam_run_bootstrap <- function(ResultBootstrap, Comparison, ECGObject) {
    ComparisonExpression <- data.frame(Comparison$control, Comparison$experiment)
    colnames(ComparisonExpression) <- c("control", "experiment")
    pblapply(
        ResultBootstrap,
        FUN = PrepareBootstrap,
        ComparisonExpression = ComparisonExpression,
        ECGObject = ECGObject
    )
}

.adam_bootstrap_to_dataframe <- function(ResultBootstrap, compID) {
    ResultBootstrap <- as.data.frame(do.call(rbind, lapply(ResultBootstrap, function(ResultBtp) {
        do.call(rbind, lapply(ResultBtp, function(ResultBtpBind) {
            as.vector(unlist(ResultBtpBind))
        }))
    })))

    colnames(ResultBootstrap) <- c(
        "Sample_Number_Genes", paste0("H_", compID[1]), paste0("H_", compID[2]),
        paste0("N_", compID[1]), paste0("N_", compID[2]), "h", "n",
        "pValue_h", "pValue_n"
    )
    GroupFields <- .adam_extract_term_fields(rownames(ResultBootstrap))
    ResultBootstrap$ID <- GroupFields$ID
    ResultBootstrap$Description <- GroupFields$Description
    rownames(ResultBootstrap) <- NULL
    ResultBootstrapColumns <- ncol(ResultBootstrap)
    ResultBootstrap[, c(length(ResultBootstrap) - 1, length(ResultBootstrap),
                        seq_len(ResultBootstrapColumns - 2))]
}

.adam_add_bootstrap_qvalues <- function(ResultBootstrap, ECGObject) {
    ResultBootstrap$qValue_h <- p.adjust(
        ResultBootstrap$pValue_h,
        method = ECGObject@PCorrectionMethod
    )
    ResultBootstrap$qValue_n <- p.adjust(
        ResultBootstrap$pValue_n,
        method = ECGObject@PCorrectionMethod
    )
    ResultBootstrap$Significance_h <- vapply(ResultBootstrap$qValue_h, function(ResultQValue) {
        ifelse(ResultQValue <= ECGObject@PCorrection, "significative", "not_significative")
    }, character(1))
    ResultBootstrap$Significance_n <- vapply(ResultBootstrap$qValue_n, function(ResultQValue) {
        ifelse(ResultQValue <= ECGObject@PCorrection, "significative", "not_significative")
    }, character(1))
    ResultBootstrap
}

.adam_reorder_bootstrap_columns <- function(ResultBootstrap) {
    ResultBootstrap <- ResultBootstrap[, c(1, 2, length(ResultBootstrap),
                                        3:(length(ResultBootstrap) - 1))]
    arrange(ResultBootstrap, ResultBootstrap[, 1])
}

.adam_apply_wilcoxon <- function(ResultBootstrap, GroupComparisons, ECGObject) {
    message("Running Wilcoxon rank-sum test ...")
    ResultBootstrap$Wilcox_pvalue <- vapply(GroupComparisons, WCAnalysis, numeric(1))
    ResultBootstrap$Wilcox_qvalue <- p.adjust(
        as.numeric(ResultBootstrap$Wilcox_pvalue),
        method = ECGObject@PCorrectionMethod
    )
    ResultBootstrap$Wilcox_significance <- ifelse(
        ResultBootstrap$Wilcox_qvalue <= ECGObject@PCorrection,
        "significative",
        "not significative"
    )
    ResultBootstrap
}

.adam_apply_fisher <- function(ResultBootstrap, GroupComparisons, Comparison, ECGObject) {
    message("Running Fisher exact test ...")
    ResultBootstrap$Fisher_pvalue <- vapply(
        GroupComparisons,
        FUN = FSHAnalysis,
        Comparison = Comparison,
        FUN.VALUE = numeric(1)
    )
    ResultBootstrap$Fisher_qvalue <- p.adjust(
        as.numeric(ResultBootstrap$Fisher_pvalue),
        method = ECGObject@PCorrectionMethod
    )
    ResultBootstrap$Fisher_significance <- ifelse(
        ResultBootstrap$Fisher_qvalue <= ECGObject@PCorrection,
        "significative",
        "not significative"
    )
    ResultBootstrap
}

.adam_finalize_bootstrap <- function(
    ResultBootstrap,
    ECGObject,
    GroupComparisons,
    Comparison,
    compID
) {
    RawLookup <- .adam_build_raw_gene_lookup(ECGObject@DBSpeciesFunctionsRaw)
    message(sprintf("Correcting bootstrap p-values by %s ...",
                    ECGObject@PCorrectionMethod))

    ResultBootstrap <- .adam_add_bootstrap_qvalues(ResultBootstrap, ECGObject)

    ResultBootstrap$Raw_Number_Genes <- vapply(ResultBootstrap$ID, function(GroupID) {
        .adam_get_raw_gene_count(GroupID, RawLookup)
    }, numeric(1))

    ResultBootstrap <- .adam_reorder_bootstrap_columns(ResultBootstrap)

    if (ECGObject@WilcoxonTest) {
        ResultBootstrap <- .adam_apply_wilcoxon(
            ResultBootstrap = ResultBootstrap,
            GroupComparisons = GroupComparisons,
            ECGObject = ECGObject
        )
    }

    if (ECGObject@FisherTest) {
        ResultBootstrap <- .adam_apply_fisher(
            ResultBootstrap = ResultBootstrap,
            GroupComparisons = GroupComparisons,
            Comparison = Comparison,
            ECGObject = ECGObject
        )
    }

    ResultBootstrap
}

makeAnalysis <- function(ComparisonID,ECGObject,completeTest){
    ParsedComparison <- .adam_build_comparison(ComparisonID, ECGObject)
    compID <- ParsedComparison$compID
    Comparison <- ParsedComparison$Comparison

    message(sprintf("Analysis: %s_Vs_%s", compID[1], compID[2]))
    message("Building GFAGs ...")

    pboptions(style = 1, char = ">")
    GroupComparisons <- .adam_build_group_comparisons(Comparison, ECGObject)

    message("Calculating activity and diversity ...")
    GroupStatistics <- pblapply(
        GroupComparisons,
        FUN = SampleStatistics,
        Control = compID[1],
        Experiment = compID[2]
    )

    if (!completeTest) {
        ResultAnalysis <- .adam_partial_analysis_result(GroupStatistics, ECGObject, compID)
    } else {
        message("Filtering GFAGs by size ...")
        ResultBootstrap <- .adam_prepare_bootstrap_groups(GroupStatistics, GroupComparisons)

        message("Running bootstrap. This may take a few minutes ...")
        ResultBootstrap <- .adam_run_bootstrap(ResultBootstrap, Comparison, ECGObject)
        ResultBootstrap <- .adam_bootstrap_to_dataframe(ResultBootstrap, compID)
        ResultAnalysis <- .adam_finalize_bootstrap(
            ResultBootstrap = ResultBootstrap,
            ECGObject = ECGObject,
            GroupComparisons = GroupComparisons,
            Comparison = Comparison,
            compID = compID
        )
    }

    message(sprintf("Analysis of %s_Vs_%s successfully concluded!",
                    compID[1], compID[2]))
    message("")

    .adam_validate_analysis_table(ResultAnalysis, compID = compID, completeTest = completeTest)
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
    Result<-data.frame(nrow(GroupComparisons),H_Control,H_Experiment,
                        N_Control, N_Experiment, h, n)
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
        length(unique(GroupData[,2][GroupData[,2]>0]))==length(GroupData[,2])
        & length(unique(c(GroupData[,1][GroupData[,1]>0],
        GroupData[,2][GroupData[,2]>0])))==(length(GroupData[,2])*2) &
        length(GroupData[,1]) > 1){
        return(wilcox.test(GroupData[,1],GroupData[,2], paired=TRUE)$p.value)
    }else{
        return(1)
    }
}

########
## Function to run the Fisher test.
########
FSHAnalysis <- function(pValuesGroups,Comparison){
    INGroupGenes <- mutate(
        pValuesGroups,
        relative = ifelse(
            (pValuesGroups[,2] + pValuesGroups[,3]) > 0,
            pValuesGroups[,3] / (pValuesGroups[,2] + pValuesGroups[,3]),
            NA_real_
        )
    )
    OUTGroupGenes <- filter(Comparison,!(as.character(Comparison[,1]) %in%
                            as.character(INGroupGenes[,1])))
    if (nrow(OUTGroupGenes) == 0) {
        return(1)
    }
    OUTGroupGenes <- mutate(
        OUTGroupGenes,
        relative = ifelse(
            (OUTGroupGenes[,2] + OUTGroupGenes[,3]) > 0,
            OUTGroupGenes[,3] / (OUTGroupGenes[,2] + OUTGroupGenes[,3]),
            NA_real_
        )
    )

    FisherCounts <- c(
        sum(INGroupGenes$relative > 0.5, na.rm = TRUE),
        sum(OUTGroupGenes$relative > 0.5, na.rm = TRUE),
        sum(INGroupGenes$relative < 0.5, na.rm = TRUE),
        sum(OUTGroupGenes$relative < 0.5, na.rm = TRUE)
    )
    if (sum(FisherCounts) == 0) {
        return(1)
    }
    FisherMatrix <- matrix(FisherCounts, 2, 2)
    if (any(!is.finite(FisherMatrix))) {
        return(1)
    }

    FisherPValue <- tryCatch(
        fisher.test(FisherMatrix)$p.value,
        error = function(...) 1
    )
    as.numeric(FisherPValue)
}
