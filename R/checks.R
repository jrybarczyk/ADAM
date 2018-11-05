##########
#Argument check functions
##########

if(getRversion() >= "3.5"){
    utils::globalVariables(c("GOID","Ontology","Term","mappedkeys",
                            "org.At.tairENTREZID","org.Sc.sgdENTREZID",
                            "org.Pf.plasmoALIAS2ORF"))
}

.checkExpressionData <- function(ExpressionData){
    
    DataFile <- ExpressionData
    
    if(is(DataFile,"SummarizedExperiment") | 
        is(DataFile,"RangedSummarizedExperiment")){
        DataFile <- as.data.frame(assay(DataFile))
        colnames(DataFile)[1] <- "gene"
        GeneColumn <- as.data.frame(as.character(DataFile$gene))
        DataFile$gene <- NULL
        DataFile <- as.matrix(DataFile)
        DataFile <- apply(DataFile,2,as.numeric)
        DataFile <- as.data.frame(DataFile)
        DataFile <- cbind.data.frame(GeneColumn,DataFile)
    }else if(!is.data.frame(DataFile)){
        if(!file.exists(DataFile)){
            stop("Please, check the format of the expression data. 
                ExpressionData argument must be a valid file or a 
                datra frame.")
        }else{
            DataFile <- read.table(ExpressionData, header = TRUE, sep = "\t",
                                    quote = "", stringsAsFactors = FALSE)
        }
    }
    colnames(DataFile)[1] <- "gene"
    
    if(ncol(DataFile)<3){
        stop("Check the data format to Conditions! The Expression data must 
            have at list 3 columns.")
    }else{
        ExpressionColumns <- DataFile    
        ExpressionColumns$gene <- NULL
        ExpressionColumns <- as.list(ExpressionColumns)
        ExpressionColumns <- sapply(ExpressionColumns,
                                function(ExpColumn) is.numeric(ExpColumn))
        if(sum(ExpressionColumns)!=(ncol(DataFile)-1)){
            stop("Check the format data! Expression values must be numeric.")
        }else if(length(unique(as.character(DataFile$gene))) != 
                length(as.character(DataFile$gene))){
            stop("Gene identifiers must be unique. Check for duplicated
            gene IDs.")
        }
    }
    return(DataFile)
}

.checkComparisonID <- function(ComparisonID,ExpressionData){
    ExpressionData <- .checkExpressionData(ExpressionData)
    ListComparisons <- as.list(ComparisonID)
    ListComparisons <- lapply(ListComparisons,function(CompID) length(unlist(
                                strsplit(CompID,","))))
    ListComparisons <- sum(ifelse(ListComparisons!=2,1,0))
    ExpDataColumns <- colnames(ExpressionData)[-1]
    if(!is.vector(ComparisonID)){
        stop("Check the format data. ComparisonID argument must be a vector!")
    }else if(ListComparisons>0){
        stop("Comparison IDs must have 2 elements!")
    }else if(sum(ifelse(unique(unlist(strsplit(ComparisonID,","))) %in% 
            ExpDataColumns,0,1)) > 0){
        stop("Invalid Comparison IDs! The IDs must correspond to the
                Expression Data column names.")
    }
    return(ComparisonID)
}


.checkGeneNumbers <- function(MinGene,MaxGene){
    if(!is.integer(MinGene) | !is.integer(MaxGene)){
        if(!is.numeric(MinGene) | !is.numeric(MaxGene)){
            stop("Please, check the format of Conditions! MinGene and MaxGene 
            arguments must be integer.")
        }else{
            MinGene <- as.integer(MinGene)
            MaxGene <- as.integer(MaxGene)
        }
    }
    if(MinGene<=0 | MaxGene<=0){
        stop("Check the Data format to Conditions! The number of genes must be
            positive and different from zero.")
    }else if(MinGene>MaxGene){
        stop("Check the Data format to Conditions! The maximum number of genes
            must be greater than the minimum number.")
    }
    return(c(MinGene,MaxGene))
}

.checkSeedNumber <- function(SeedNumber){
    if(!is.numeric(SeedNumber)){
        stop("Please, check the format of Conditions! SeedNumber argument must
            be numeric.")
    }else if(SeedNumber<0){
            stop("Check the Data format Conditions! The seed must be a positive
                value.")
    }
    return(SeedNumber)
}

.checkBootstrapNumber <- function(BootstrapNumber){
    if(!is.integer(BootstrapNumber)){
        if(!is.numeric(BootstrapNumber)){
            stop("Please, check the format of Conditions! BootstrapNumber
                argument must be integer.")
        }else{
            BootstrapNumber<-as.integer(BootstrapNumber)
        }
    }
    if(BootstrapNumber<=0){
        stop("Check the Data format Conditions! The number of bootstrap steps
            genes must be positive and different from zero.")
    }
    return(BootstrapNumber)
}

.checkPCorrection <- function(PCorrection){
    if(!is.numeric(PCorrection)){
        stop("Please, check the format of Conditions! PCorrection argument
            must be numeric.")
    }else if(PCorrection<0 || PCorrection>1){
        stop("Check the Data format Conditions! Correction value must be 
            between zero and one.")
    }
    return(PCorrection)
}

.checkPCorrectionMethod <- function(PCorrectionMethod){
    options <- c("holm", "hochberg", "hommel", "bonferroni", "bh", "by","fdr")
    OptionsPAdjust <- c("holm", "hochberg", "hommel", "bonferroni", "BH",
                        "BY","fdr")
    if  (!(tolower(PCorrectionMethod) %in% options)){
        stop(paste("Correction method",PCorrectionMethod,"not found."))
    }
    return(OptionsPAdjust[grep(tolower(PCorrectionMethod),options)])
}

.checkWilcoxonTest <- function(WilcoxonTest){
    if(!is.logical(WilcoxonTest)){
        stop("Please, check the format of Conditions! WilcoxonTest argument 
            must be logical.")
    }
    return(WilcoxonTest)
}

.checkFisherTest <- function(FisherTest){
    if(!is.logical(FisherTest)){
        stop("Please, check the format of Conditions! FisherTest argument must
            be logical.")
    }
    return(FisherTest)
}

###
# Creates a data frame relating Species and gene annotation package. This 
# data frame is used by .checkAnalysisDomain_DBSpecies_GeneIdentifier funcion.
###
.CreateDataSpecies <-function(){
    SpeciesDescription <- c("Anopheles gambiae", "Arabdopsis thaliana",
                        "Bos taurus", "Caenorhabditis elegans",
                        "Canis familiaris", "Drosophila melanogaster",
                        "Danio rerio", "Escherichia coli K12",
                        "Escherichia coli Sakai", "Gallus gallus",
                        "Homo sapiens", "Mus musculus", "Macaca mulatta",
                        "Plasmodium falciparum", "Pan troglodytes",
                        "Rattus norvegicus", "Saccharomyces cerevisiae",
                        "Sus scrofa","Xenopus laevis")
    SpeciesID <- c("aga","ath","bta","cel","cfa","dme","dre","eco","ecs",
                "gga","hsa","mmu","mcc","pfa","ptr","rno","sce","ssc","xla")
    DBSpeciesList <- c("org.Ag.eg.db", "org.At.tair.db", "org.Bt.eg.db",
                    "org.Ce.eg.db", "org.Cf.eg.db","org.Dm.eg.db",
                    "org.Dr.eg.db", "org.EcK12.eg.db", "org.EcSakai.eg.db",
                    "org.Gg.eg.db", "org.Hs.eg.db", "org.Mm.eg.db",
                    "org.Mmu.eg.db", "org.Pf.plasmo.db", "org.Pt.eg.db",
                    "org.Rn.eg.db", "org.Sc.sgd.db", "org.Ss.eg.db",
                    "org.Xl.eg.db")
    
    DataFrameSpecies <- data.frame(DBSpeciesList, SpeciesDescription,
                                    SpeciesID)
    return(DataFrameSpecies)
}

.checkAnalysisDomain_DBSpecies_GeneIdentifier <- function(AnalysisDomain,
                                                        DBSpecies,
                                                        GeneIdentifier,
                                                        GeneExpressionData){
    #Creates a data frame relating Species and gene annotation packages
    DataFrameSpecies <- .CreateDataSpecies()
    
    #Available domains
    options <- c("gobp","gocc","gomf","kegg","own")
    
    #Available nomenclatures
    GeneNomenclatures <- c("symbol","entrez","tair","orf")
    
    SpeciesID <- DBSpecies
    if(!(tolower(AnalysisDomain) %in% options)){
        stop(paste("Analysis domain",AnalysisDomain,"does not exist!"))
    }else if(AnalysisDomain!=options[5]){
        if(!is.data.frame(SpeciesID)){
            if(!file.exists(SpeciesID)){
                if(!SpeciesID %in% DataFrameSpecies$DBSpeciesList){
                    if(tolower(SpeciesID)=="eck12"){
                        SpeciesID <- "EcK12"
                    }else if(tolower(SpeciesID)=="ecsakai"){
                        SpeciesID <- "EcSakai"
                    }else if(tolower(SpeciesID)=="mmu"){
                        SpeciesID <- "Mmu"
                    }else{
                        SpeciesID <- unlist(strsplit(SpeciesID,""))
                        SpeciesID <- paste0(toupper(SpeciesID[1]),
                            paste0(tolower(SpeciesID[2:length(SpeciesID)]),
                            collapse = ""))
                    }
                    SpeciesID <- ifelse(length(DataFrameSpecies$DBSpeciesList[
                        grep(SpeciesID,DataFrameSpecies$DBSpeciesList)])==0,
                        DBSpecies,
                        as.character(DataFrameSpecies$DBSpeciesList[
                        grep(SpeciesID,DataFrameSpecies$DBSpeciesList)]))
                    if(!SpeciesID %in% DataFrameSpecies$DBSpeciesList){
                        stop("Gene annotation package does not supported! 
                            Please check it.")
                    }
                }
            }else{
                stop("Please inform a valid annotation package name!")
            }
        }else{
            stop("Please inform a valid annotation package name!")
        }
    }

    ExpressionData <- GeneExpressionData
    if(!(is.data.frame(GeneExpressionData))){
        ExpressionData <- read.table(GeneExpressionData, sep = "\t",
                                    header = TRUE, quote = "",
                                    stringsAsFactors = FALSE)
    }
    colnames(ExpressionData)[1] <- "gene"
    
    if(tolower(AnalysisDomain)==options[5]){
        if(!is.data.frame(SpeciesID)){
            if(!file.exists(SpeciesID)){
                stop("Own functions file does not exist!")
            }else{
                DBFile <- read.table(SpeciesID, header = TRUE,
                                    sep = "\t", quote = "",
                                    stringsAsFactors = FALSE)
                colnames(DBFile) <- c("gene","ID","Description")
                DBFile <- .checkDB(DBFile = DBFile,
                        ExpressionData = ExpressionData)
            }
        }else{
            DBFile <- SpeciesID
            colnames(DBFile) <- c("gene","ID","Description")
            DBFile <- .checkDB(DBFile=DBFile,
                                ExpressionData = ExpressionData)
        }
    }else if(.checkPackage(NamePackage = SpeciesID)){
            if(!GeneIdentifier %in% GeneNomenclatures){
                stop("Gene identifier invalid!")
            }else if(GeneIdentifier==GeneNomenclatures[3] & 
                    SpeciesID!="org.At.tair.db"){
                stop("Incorrect GeneIdentifier!")
            }else if(GeneIdentifier==GeneNomenclatures[4] & 
                    SpeciesID!="org.Sc.sgd.db" & 
                    SpeciesID!="org.Pf.plasmo.db"){
                    stop("Incorrect GeneIdentifier!")
            }else{
                if(SpeciesID=="org.At.tair.db"){
                    SymbolEGClass <- org.At.tairENTREZID #EntrezID and tairID
                    SymbolEGList <- unlist(as.list(SymbolEGClass[mappedkeys(
                                    SymbolEGClass)]))
                    SymbolEGList <-data.frame(tair=names(SymbolEGList),
                                            entrez=as.vector(SymbolEGList))
            }else if(SpeciesID=="org.Sc.sgd.db"){
                    SymbolEGClass <- org.Sc.sgdENTREZID #EntrezID/ORF numbers
                    SymbolEGList <- unlist(as.list(SymbolEGClass[mappedkeys(
                                    SymbolEGClass)]))
                    SymbolEGList <-data.frame(orf=names(SymbolEGList),
                                    entrez=as.vector(SymbolEGList))
            }else if(SpeciesID=="org.Pf.plasmo.db"){
                    SymbolEGClass<-org.Pf.plasmoALIAS2ORF #GeneSymbol / ORF ID
                    SymbolEGList <- unlist(as.list(SymbolEGClass[mappedkeys(
                                    SymbolEGClass)]))
                    SymbolEGList <-data.frame(symbol=names(SymbolEGList),
                                    orf=as.vector(SymbolEGList))
            }else{
                SymbolEGClass <- get(paste0("org.",unlist(strsplit(SpeciesID,
                                "[.]"))[2],".egSYMBOL2EG")) #GeneSymbol/Entrez
                SymbolEGList <- unlist(as.list(SymbolEGClass[mappedkeys(
                                SymbolEGClass)]))
                SymbolEGList <-data.frame(symbol=names(SymbolEGList),
                                entrez=as.vector(SymbolEGList))
            }

            if(sum(ifelse(ExpressionData$gene %in% 
                            SymbolEGList[,c(grep(tolower(GeneIdentifier),
                            colnames(SymbolEGList)))] |
                            tolower(ExpressionData$gene) %in% 
                            SymbolEGList[,c(grep(tolower(GeneIdentifier),
                            colnames(SymbolEGList)))] |
                            toupper(ExpressionData$gene) %in% 
                            SymbolEGList[,c(grep(tolower(GeneIdentifier),
                            colnames(SymbolEGList)))],1,0)) < 1){
                stop("The genes in the Expression Data are not in the species 
                    database. Check the gene identifiers!")
            }else{
                DBFile <- .checkGenesDomain(AnalysisDomain = 
                                        tolower(AnalysisDomain), 
                                        DBSpecies = SpeciesID, 
                IDGene = tolower(GeneIdentifier), 
                ListGenes = SymbolEGList,
                keggID=as.character(DataFrameSpecies$SpeciesID[grep(SpeciesID,
                                    DataFrameSpecies$DBSpeciesList)]),
                                    ExpressionData = ExpressionData)
            }
        }
    }

    return(list(options[grep(tolower(AnalysisDomain),options)],DBFile[[1]],
            DBFile[[2]], tolower(GeneIdentifier)))
}

.checkPackage <- function(NamePackage){
    if(!suppressMessages(requireNamespace(NamePackage))){
        stop(paste("Annotation package",NamePackage,"not availlable.
                    Please install it!"))
    }else{
        suppressMessages(require(NamePackage, character.only = TRUE,
                                warn.conflicts = FALSE))
    }
    return(TRUE)
}

.checkDB <- function(DBFile,ExpressionData){
    if(ncol(DBFile)!=3){
        stop("Own functions file must have 3 columns!")
    }else if(nrow(DBFile)<1){
        stop("Own functions file is empty!")
    }else if(sum(ifelse(ExpressionData$gene %in% DBFile$gene | 
        tolower(ExpressionData$gene) %in% DBFile$gene |
        toupper(ExpressionData$gene) %in% DBFile$gene,1,0)) < 1){
        stop("The genes in the Expression Data are not in the species
            database. Check the gene identifiers!")
        }else if(length(colnames(DBFile))==length(colnames(ExpressionData))){
                if(sum(colnames(DBFile)==colnames(ExpressionData))==3){
                    stop("Check the ExpressionData and DBSpecies arguments.
                        The column names of both files are the same!!")
                }
    }

    DBFunctionsList <- lapply(as.vector(unique(as.character(DBFile$ID))), 
                        FUN = .ListConversionDB, DBFile = DBFile)
    names(DBFunctionsList) <- as.vector(unique(as.character(DBFile$ID)))
    TermOntologies <- DBFile
    TermOntologies$gene <- NULL 
    TermOntologies <- unique(TermOntologies)

    DBFunctionsListSample <- lapply(DBFunctionsList, function(ListElement,
                                    SampleGenes) {SampleGenes[SampleGenes %in% 
                                    ListElement]},
                                    SampleGenes = ExpressionData$gene)
    names(DBFunctionsListSample) <- vapply(names(DBFunctionsList),
                function(ListElement,TermOntologies){
                SubTermOntologies<-subset(TermOntologies,
                            TermOntologies$ID == ListElement)
                return(paste0(as.character(SubTermOntologies$ID),
                    " <==> ",as.character(SubTermOntologies$Description)))},
                TermOntologies=TermOntologies, FUN.VALUE = character(1))

    DBFunctionsListRaw <- lapply(DBFunctionsList,function(ListElement){
                                length(ListElement)})
    names(DBFunctionsListRaw) <- names(DBFunctionsList)

    return(list(DBFunctionsListSample,DBFunctionsListRaw))
}

.ListConversionDB <- function(ListElement,DBFile){
    ListResult <- subset(DBFile, DBFile$ID == ListElement)
    ListResult <- ListResult$gene
    return(ListResult)
}

.checkGenesDomain<-function(AnalysisDomain,DBSpecies,IDGene,ListGenes,keggID,
                            ExpressionData){

    #Each annotation package has its own set of functions that relates genes
    #to ontologies and kegg pathways. This way, we create a list linking the 
    #package to its functions. Besides that, the list also presents the 
    #default and alternative nomenclature supported according to each DB 
    #species. For example, let's consider Arabdopsis thaliana: List element
    #name: org.At.tair.db; Function to relate genes to gene ontologies terms:
    #org.At.tairGO2ALLTAIRS; Default gene nomenclature: "tair" ID;
    #Alternative gene nomenclature: "entrez" ID; Function to relate genes to
    #kegg pathways: org.At.tairARACYC.
    ListDBPackages <- list(org.At.tair.db = c("org.At.tairGO2ALLTAIRS","tair",
                                            "entrez","org.At.tairARACYC"),
                org.Sc.sgd.db = c("org.Sc.sgdGO2ALLORFS","orf",
                                "entrez","org.Sc.sgdPATH"),
                org.Pf.plasmo.db = c("org.Pf.plasmoGO2ALLORFS","orf",
                                "symbol","org.Pf.plasmoPATH"),
                otherDB=c(paste0("org.",unlist(strsplit(DBSpecies,"[.]"))[2],
                        ".egGO2ALLEGS"),"entrez","symbol",
                        paste0("org.",unlist(strsplit(DBSpecies,"[.]"))[2],
                        ".egPATH2EG")))

    if(DBSpecies %in% names(ListDBPackages)){
        OntologyFunction <- c(ListDBPackages[[grep(DBSpecies,names(
                            ListDBPackages))]][1],
                            ListDBPackages[[grep(DBSpecies,names(
                            ListDBPackages))]][2],
                            ListDBPackages[[grep(DBSpecies,names(
                            ListDBPackages))]][3],
                            ListDBPackages[[grep(DBSpecies,names(
                            ListDBPackages))]][4])
    }else{
        OntologyFunction <- c(ListDBPackages[[4]][1],ListDBPackages[[4]][2],
        ListDBPackages[[4]][3],ListDBPackages[[4]][4])
    }

    if(AnalysisDomain == "kegg"){
        keggPathways <- as.data.frame(keggList("pathway", keggID[2]))
        colnames(keggPathways)[1] <- "ID"
        keggPathwaysID <- unlist(lapply(rownames(keggPathways),
                        function(KPathway)
                        unlist(strsplit(KPathway,paste0(":",keggID[2])))[2]))
        
        keggPathwaysDescription<-unlist(lapply(as.character(keggPathways$ID),
                                function(KPathway) unlist(strsplit(KPathway,
                                paste0("- ",keggID[1])))[1]))
        
        TermOntologies <- data.frame(ID = keggPathwaysID, 
                                    Description = keggPathwaysDescription)

        DBKeggClass <- get(OntologyFunction[4])
        DBFunctionsList <- as.list(DBKeggClass[mappedkeys(DBKeggClass)])

        DBFunctionsList <- lapply(DBFunctionsList, function(DBFunction) 
                                unique(as.vector(DBFunction)))
        DBFunctionsList <- mapply(.CheckTerm, ListElement = DBFunctionsList, 
                                ListElementName = names(DBFunctionsList),
                                MoreArgs =  list(VectorTerms = 
                                as.character(TermOntologies$ID)))

        TermOntologies <- TermOntologies[TermOntologies$ID %in% 
                            names(DBFunctionsList),]
        TermOntologies$ID <- paste0(keggID[2],TermOntologies$ID)

        DBFunctionsList<-DBFunctionsList[lapply(DBFunctionsList, length) > 0]
        names(DBFunctionsList) <- paste0(keggID[2],names(DBFunctionsList))
        
        
    }else{
        infoGO <- GOTERM
        TermOntologies <- data.frame(t(vapply(infoGO, function(GOntology) 
                    c(trimws(GOID(GOntology)), trimws(Ontology(GOntology)), 
                    trimws(Term(GOntology))), FUN.VALUE = character(3))))

        TermOntologies <- subset(TermOntologies, TermOntologies$X2 == 
                            toupper(unlist(strsplit(AnalysisDomain,"go"))[2]))
        colnames(TermOntologies) <- c("ID","Domain","Description")
        TermOntologies$Domain <- NULL
        rownames(TermOntologies) <- NULL

        DBFunctions <- get(OntologyFunction[1])
        DBFunctionsList <- as.list(DBFunctions[mappedkeys(DBFunctions)])
        DBFunctionsList <- lapply(DBFunctionsList, function(DBFunction)
                            unique(as.vector(DBFunction)))
        DBFunctionsList <- mapply(.CheckTerm, ListElement = DBFunctionsList,
                            ListElementName = names(DBFunctionsList),
                            MoreArgs =  list(VectorTerms = 
                            as.character(TermOntologies$ID)))
        DBFunctionsList<-DBFunctionsList[lapply(DBFunctionsList, length) > 0]
    }

    if(IDGene != OntologyFunction[2]){
        DBFunctionsList <- lapply(DBFunctionsList, FUN = .GeneIDConversion, 
                            ListGenes = ListGenes,
                            StandardGeneID = OntologyFunction[2], 
                            AlternativeGeneID = OntologyFunction[3])
    }

    TermOntologies <- subset(TermOntologies, TermOntologies$ID %in%
                        names(DBFunctionsList))

    DBFunctionsListSample <- lapply(DBFunctionsList, 
                            function(ListElement,SampleGenes){
                            SampleGenes[SampleGenes %in% 
                            ListElement]},
                            SampleGenes = ExpressionData$gene)
    
    names(DBFunctionsListSample) <- vapply(names(DBFunctionsList),
                function(ListElement,TermOntologies){
                SubTermOntologies<-subset(TermOntologies, 
                            TermOntologies$ID==ListElement)
                return(paste0(as.character(SubTermOntologies$ID),
                    " <==> ",as.character(SubTermOntologies$Description)))},
                TermOntologies=TermOntologies, FUN.VALUE = character(1))

    DBFunctionsListRaw <- lapply(DBFunctionsList,function(ListElement){
                                length(ListElement)})
    names(DBFunctionsListRaw) <- names(DBFunctionsList)

    DBFunctionsListSample<-DBFunctionsListSample[lapply(DBFunctionsListSample,
                                                length) > 0]

    return(list(DBFunctionsListSample,DBFunctionsListRaw))
}

.CheckTerm <- function(ListElement,ListElementName,VectorTerms){
    if(as.character(ListElementName) %in% as.character(VectorTerms)){
        result <- ListElement
    }else{
        result <- NULL
    }
    return(result)
}

.GeneIDConversion <- function(ListElement,ListGenes,StandardGeneID,
                            AlternativeGeneID){
    SubListGenes <- subset(ListGenes, ListGenes[,grep(StandardGeneID,
                    colnames(ListGenes))] %in% ListElement)
    return(as.character(SubListGenes[,grep(AlternativeGeneID,
        colnames(SubListGenes))]))
}
