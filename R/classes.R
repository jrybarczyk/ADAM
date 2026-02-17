setClass(
    "ECGMainData",
    slots = c(
        ComparisonID = "character",
        ExpressionData = "data.frame",
        MinGene = "integer",
        MaxGene = "integer",
        SeedNumber = "numeric",
        BootstrapNumber = "integer",
        PCorrection = "numeric",
        DBSpeciesFunctionsSample = "list",
        DBSpeciesFunctionsRaw = "list",
        PCorrectionMethod = "character",
        WilcoxonTest = "logical",
        FisherTest = "logical",
        AnalysisDomain = "character",
        GeneIdentifier = "character"
    ),
    prototype = list(
        ComparisonID = character(),
        ExpressionData = data.frame(),
        MinGene = as.integer(NA),
        MaxGene = as.integer(NA),
        SeedNumber = NA_real_,
        BootstrapNumber = as.integer(NA),
        PCorrection = NA_real_,
        DBSpeciesFunctionsSample = list(),
        DBSpeciesFunctionsRaw = list(),
        PCorrectionMethod = NA_character_,
        WilcoxonTest = NA,
        FisherTest = NA,
        AnalysisDomain = NA_character_,
        GeneIdentifier = NA_character_
    )
)

setValidity("ECGMainData", function(object) {
    if (length(object@MinGene) != 1L || is.na(object@MinGene) || object@MinGene <= 0L) {
        return("MinGene must be a positive integer scalar.")
    }
    if (length(object@MaxGene) != 1L || is.na(object@MaxGene) || object@MaxGene <= 0L) {
        return("MaxGene must be a positive integer scalar.")
    }
    if (object@MinGene > object@MaxGene) {
        return("MinGene must be less than or equal to MaxGene.")
    }
    if (!is.na(object@PCorrection) &&
        (length(object@PCorrection) != 1L || object@PCorrection < 0 || object@PCorrection > 1)) {
        return("PCorrection must be a numeric scalar between 0 and 1, or NA.")
    }
    if (!is.na(object@BootstrapNumber) &&
        (length(object@BootstrapNumber) != 1L || object@BootstrapNumber <= 0L)) {
        return("BootstrapNumber must be a positive integer scalar, or NA.")
    }
    if (!is.na(object@SeedNumber) &&
        (length(object@SeedNumber) != 1L || object@SeedNumber < 0)) {
        return("SeedNumber must be a non-negative numeric scalar, or NA.")
    }
    TRUE
})
