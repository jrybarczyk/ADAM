#' @title Group of Functionally Associated Genes (GFAG) analysis
#' @rdname GFAGAnalysis
#' @import Rcpp GO.db KEGGREST
#' @importFrom SummarizedExperiment assay
#' @import pbapply methods knitr
#' @importFrom utils installed.packages
#' @importFrom Rcpp evalCpp
#' @importFrom stats wilcox.test fisher.test
#' @importFrom DT datatable
#' @importFrom dplyr arrange mutate filter
#' @useDynLib ADAM
#' @importFrom stats p.adjust
#' @importFrom utils read.table write.table
#' @usage GFAGAnalysis(ComparisonID, ExpressionData, MinGene, MaxGene, 
#' SeedNumber, BootstrapNumber, PCorrection, DBSpecies, PCorrectionMethod, 
#' WilcoxonTest, FisherTest, AnalysisDomain, GeneIdentifier)
#' @description Analysis of functionally associated gene groups, based on
#' gene diversity and activity, for different species according to existing 
#' annotation packages or user's personal annotations.
#' @param ComparisonID Sample comparisons identification. It must be a vector
#' in which each element corresponds to 2 sample columns from the expression 
#' data. The data's sample columns in each element from the vector are comma 
#' separated. The first one refers to the control sample, while the second
#' refers to the experiment.
#' @param ExpressionData Gene expression data (microarray or RNA-seq, for
#' example). It must be a SummarizedExperiment object, a data frame or a path
#' for a text file tab separated containing at least 3 columns. First column 
#' mandatory corresponds to the gene identification, according to 
#' GeneIdentifier argument. Second, third and the other columns correspond to
#' the gene expression values realated to the genes in the first column and
#' each of these columns correspond to a different sample (control versus 
#' experiment).
#' @param MinGene Minimum number of genes per GFAG. It must be a positive
#' integer value different from zero and lower than MaxGene argument.
#' Default is 3.
#' @param MaxGene Maximum number of genes per GFAG. It must be a positive
#' integer value different from zero and greater than MinGene argument.
#' Default is 2000.
#' @param SeedNumber Seed for bootstrap. A numeric positive value used as seed
#' for random number generating by bootstrap function. Default is 1049.
#' @param BootstrapNumber Number of bootstraps. A numeric value greater than 
#' zero, used by bootstrap function generates p-values for each GFAG. Default
#' is 1000. 
#' @param PCorrection Cutoff for p-value correction. A numeric value between 
#' 0 and 1. Default is 0.05.
#' @param DBSpecies A string corresponding to the name of an OrgDb species 
#' gene annotation package: org.Ag.eg.db (Anopheles gambiae), org.At.tair.db 
#' (Arabdopsis thaliana), org.Bt.eg.db (Bos taurus), org.Ce.eg.db
#' (Caenorhabditis elegans), org.Cf.eg.db (Canis familiaris), org.Dm.eg.db 
#' (Drosophila melanogaster), org.Dr.eg.db (Danio rerio), org.EcK12.eg.db 
#' (Escherichia coli K12), org.EcSakai.eg.db (Escherichia coli Sakai),
#' org.Gg.eg.db (Gallus gallus), org.Hs.eg.db (Homo sapiens), org.Mm.eg.db
#' (Mus musculus), org.Mmu.eg.db (Macaca mulatta), org.Pf.plasmo.db
#' (Plasmodium falciparum), org.Pt.eg.db (Pan troglodytes), org.Rn.eg.db
#' (Rattus norvegicus), org.Sc.sgd.db (Saccharomyces cerevisiae),
#' org.Ss.eg.db (Sus scrofa) and org.Xl.eg.db (Xenopus laevis).
#' If there is no package, it's possible for the user to create a personal 
#' gene annotation file, tab separated, containing 3 columns: gene, term 
#' annotation code and description of the term annotation. So, istead of a 
#' string with an OrgDb name, inform a data frame or a path for the file.
#' @param PCorrectionMethod Method p-value correction: holm, hochberg, hommel,
#' bonferroni, bh, by or fdr. Default is "fdr".
#' @param WilcoxonTest A logical value indicating whether or not to perform
#' Wilcoxon test (TRUE for performing and FALSE for not performing).
#' Default is FALSE.
#' @param FisherTest  A logical value indicating whether or not to perform
#' Fisher's exact test (TRUE for performing and FALSE for not performing).
#' Default is FALSE.
#' @param AnalysisDomain Analysis domain to be considered for building GFAGs, 
#' according: gobp (Gene Ontology - Biological Processes), gocc (Gene Ontology
#'  - Celular Components), gomf (Gene Ontology - Molecular Functions), kegg
#' (KEGG Pathways) or own (if there is no annotation package - the annotations 
#' were defined in a file by the user).
#' @param GeneIdentifier Gene nomenclature to be used: symbol, entrez, tair
#' (for Arabdopsis thaliana) or orf (for Saccharomyces cerevisiae or 
#' Plasmodium falciparum). If there is no annotation package, just put the 
#' gene nomenclature present in the user's personal annotations.
#' @details The genes present in the expression data are grouped by their 
#' respective functions according to the domains described by 
#' AnalysisDomain argument. The relationship between genes and functions are 
#' made based on the species annotation package. If there is no annotation 
#' package, a three column file (gene, function and function description) must
#' be provided. For each GFAG, gene diversity and activity in each sample are 
#' calculated. As the package always compare two samples (control versus 
#' experiment), relative gene diversity and activity for each GFAG are
#' calculated. Using bootstrap method, for each GFAG, according to relative 
#' gene diversity and activity, two p-values are calculated. The p-values are 
#' then corrected, according to the correction method defined by 
#' PCorrectionMethod argument, generating a q-value. The significative GFAGs 
#' will be those whoose q-value stay under the cutoff set by PCorrection 
#' argument. Optionally, it's possible to run Wilcoxon test and/or Fisher's
#' exact test. These tests also provide a corrected p-value, 
#' and siginificative groups can be seen through them. GFAGAnalysis function 
#' allows to run a complete analysis, using all available arguments. 
#' ADAnalysis function allows to run a partial analysis, where is calculated
#' just gene diversity and activity of each GFAG with no signicance by 
#' bootrstrap, Wilcoxon or Fisher.
#' @author André Luiz Molan (andre.molan@unesp.br)
#' @references CASTRO, M. A., RYBARCZYK-FILHO, J. L., et al. Viacomplex:
#' software for landscape analysis of gene expression networks in genomic
#' context. Bioinformatics, Oxford Univ
#' Press, v. 25, n. 11, p. 1468–1469, 2009.
#' @return Return a list with two elements. The first one refers to a 
#' data frame with the GFAGs and their respective genes. The second one is a 
#' a list where each position is a data frame presenting the result of the 
#' analysis, according to ComparisonID argument.
#' @examples
#' ##
#' #Complete Analysis with Aedes aetypti through GFAGAnalysis function
#' ##
#' data(ExpressionAedes)
#' data(KeggPathwaysAedes)
#' ResultAnalysis <- GFAGAnalysis(ComparisonID = c("control1,experiment1", 
#' "control2,experiment2"), ExpressionData = ExpressionAedes, MinGene = 3L,
#' MaxGene = 20L, SeedNumber = 1049, BootstrapNumber = 1000L,
#' PCorrection = 0.05, DBSpecies = KeggPathwaysAedes,
#' PCorrectionMethod = "fdr", WilcoxonTest = TRUE, FisherTest = TRUE,
#' AnalysisDomain = "own", GeneIdentifier = "gene")
#' \dontrun{
#' head(ResultAnalysis[[1]]) #Relation between genes and functions
#' head(ResultAnalysis[[2]][1]) #Result comparison 1
#' head(ResultAnalysis[[2]][2]) #Result comparison 2
#' }
#' ##
#' #Partial Analysis with Aedes aegypti through ADAnalysis function
#' ##
#' data(ExpressionHs)
#' ResultAnalysis <- ADAnalysis(ComparisonID = c("control1,experiment1"),
#' ExpressionData = ExpressionAedes, MinGene = 3L, MaxGene = 20L, 
#' DBSpecies = KeggPathwaysAedes, AnalysisDomain = "own",
#' GeneIdentifier = "geneStableID")
#' \dontrun{
#' head(ResultAnalysis[[1]]) #Relation between genes and functions
#' head(ResultAnalysis[[2]][1]) #Result comparison 1
#' }
#' @export

GFAGAnalysis <- function(ComparisonID, ExpressionData, MinGene = 3,
                        MaxGene = 2000, SeedNumber = 1049,
                        BootstrapNumber = 1000, PCorrection = 0.05,
                        DBSpecies, PCorrectionMethod = "fdr",
                        WilcoxonTest = FALSE, FisherTest = FALSE,
                        AnalysisDomain, GeneIdentifier){
    message("Creating object ...")
    ECGObject <- ECGMainData(ComparisonID = ComparisonID, ExpressionData =
                                ExpressionData, MinGene = MinGene,
                                MaxGene = MaxGene,
                                SeedNumber = SeedNumber, BootstrapNumber = 
                                BootstrapNumber, PCorrection = PCorrection,
                                DBSpecies = DBSpecies, PCorrectionMethod = 
                                PCorrectionMethod, WilcoxonTest=WilcoxonTest,
                                FisherTest = FisherTest, AnalysisDomain =
                                AnalysisDomain,
                                GeneIdentifier = GeneIdentifier,
                                completeTest = TRUE)
    #GFAG Analysis
    ResultAnalysis <- lapply(ECGObject@ComparisonID, FUN = makeAnalysis,
                            ECGObject = ECGObject, completeTest = TRUE)
    
    #File relating genes and functions
    GenesFile <- as.data.frame(do.call(rbind,
                lapply(ECGObject@DBSpeciesFunctionsSample,
                    function(DBFunctionsSample){
                    as.data.frame(as.vector(unlist(DBFunctionsSample)))})))
    GenesFile$GroupID  <- row.names(GenesFile)
    GroupID <- data.frame(do.call('rbind',
                        strsplit(as.character(GenesFile$GroupID),
                        '<==>',fixed=TRUE)))
    rownames(GenesFile) <- NULL
    colnames(GenesFile) <- c("gene","GroupID")
    GenesFile$GroupID <- trimws(GroupID$X1)
    GenesFile <- unique(GenesFile)
    rm(GroupID)
    
    return(list(GenesFile,ResultAnalysis))
}

#' @rdname GFAGAnalysis
#' @usage ADAnalysis(ComparisonID, ExpressionData, MinGene, MaxGene, 
#' DBSpecies, AnalysisDomain, GeneIdentifier)
#' @export
ADAnalysis <- function(ComparisonID, ExpressionData, MinGene=3, MaxGene=2000, 
                    DBSpecies, AnalysisDomain, GeneIdentifier){
    message("Creating object ...")
    ECGObject <- ECGMainData(ComparisonID = ComparisonID, 
                            ExpressionData = ExpressionData, 
                            MinGene = MinGene,
                            MaxGene = MaxGene, DBSpecies = DBSpecies, 
                            AnalysisDomain = AnalysisDomain,
                            GeneIdentifier = GeneIdentifier,
                            completeTest = FALSE)
    
    #GFAG Analysis
    ResultAnalysis <- lapply(ECGObject@ComparisonID, FUN = makeAnalysis,
                            ECGObject = ECGObject, completeTest = FALSE)
    
    #File relating genes and functions
    GenesFile <- as.data.frame(do.call(rbind,
                    lapply(ECGObject@DBSpeciesFunctionsSample,
                    function(DBFunctionsSample){
                    as.data.frame(as.vector(unlist(DBFunctionsSample)))})))
    GenesFile$GroupID  <- row.names(GenesFile)
    GroupID <- data.frame(do.call('rbind', 
                                strsplit(as.character(GenesFile$GroupID),
                                '<==>',fixed=TRUE)))
    rownames(GenesFile) <- NULL
    colnames(GenesFile) <- c("gene","GroupID")
    GenesFile$GroupID <- trimws(GroupID$X1)
    GenesFile <- unique(GenesFile)
    rm(GroupID)
    
    return(list(GenesFile,ResultAnalysis))
}
