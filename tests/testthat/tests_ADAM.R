### Tests
test_that("data size from ADAnalysis is ok", {
    expect_silent(data(ExpressionAedes))
    expect_silent(data(KeggPathwaysAedes))
    expect_equal(length(as.data.frame(suppressMessages(
            ADAnalysis(ComparisonID = c("control1,experiment1"),
            ExpressionData = ExpressionAedes, MinGene = 3,
            DBSpecies = KeggPathwaysAedes, MaxGene = 1000,
            AnalysisDomain = "own",
            GeneIdentifier = "geneStableID"))[[2]])$ID),
            length(unique(KeggPathwaysAedes$pathwayID)))
})


test_that("data size from GFAGAnalysis is ok", {
    expect_silent(data(ExpressionAedes))
    expect_silent(data(KeggPathwaysAedes))
    expect_equal(length(as.data.frame(suppressMessages(
            GFAGAnalysis(ComparisonID  = c("control1,experiment1"),
            ExpressionData = ExpressionAedes, MinGene = 3, MaxGene = 20, 
            SeedNumber = 1049, BootstrapNumber = 1000, PCorrection = 0.05,
            DBSpecies = KeggPathwaysAedes, PCorrectionMethod = "fdr", 
            WilcoxonTest = FALSE, FisherTest = FALSE,
            AnalysisDomain = "own",GeneIdentifier="geneStableID"))[[2]])$ID),
            length(unique(KeggPathwaysAedes$pathwayID)))
})


