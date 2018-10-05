##########
#Argument check functions
##########
#x = ComparisonID
#y = ExpressionData
checkComparisonID <- function(x,y){
    ListComparisons <- as.list(x)
    ListComparisons <- lapply(ListComparisons,function(k) length(unlist(
                                strsplit(k,","))))
    ListComparisons <- sum(ifelse(ListComparisons!=2,1,0))
    if(!is.vector(x)){
        stop("Check the format data. ComparisonID argument must be a vector!")
    }else if(ListComparisons>0){
        stop("Comparison IDs must have 2 elements!")
    }else if(sum(ifelse(unique(unlist(strsplit(x,","))) %in% 
            colnames(y)[2:ncol(y)],0,1)) > 0){
        stop("Invalid Comparison IDs! The IDs must correspond to the Expression
                Data column names.")
    }
    return(x)
}

#x = ExpressionData
checkExpressionData <- function(x){
    DataFile <- x
    if(!is.data.frame(DataFile)){
        if(!file.exists(DataFile)){
            stop("Please, check the format of the expression data. 
                ExpressionData argument must be a valid file or a 
                datra.frame.")
        }else{
            DataFile <- read.table(x, header = TRUE, sep = "\t", quote = "", 
                        stringsAsFactors = FALSE)
        }
    }
    if(ncol(DataFile)<3){
        stop("Check the data format to Conditions! The Expression data must 
            have at list 3 columns.")
    }else{
    ExpressionColumns <- as.list(DataFile[c(2:ncol(DataFile))])
    ExpressionColumns <- lapply(ExpressionColumns,function(y) is.numeric(y))
    if(sum(ifelse(ExpressionColumns==TRUE,1,0))!=(ncol(DataFile)-1)){
        stop("Check the format data! Expression values must be numeric.")
    }else if(length(unique(as.character(DataFile[,1]))) != 
               length(as.character(DataFile[,1]))){
        stop("Gene identifiers must be unique. Check for duplicated
            gene IDs.")
        }
    }
    return(DataFile)
}

#x = MinGene
#y = MaxGene
checkGeneNumbers <- function(x,y){
    if(!is.integer(x) | !is.integer(y)){
        stop("Please, check the format of Conditions! MinGene and MaxGene 
            arguments must be integer.")
    }else if(x<=0 | y<=0){
        stop("Check the Data format to Conditions! The number of genes must be
            positive and different from zero.")
    }else if(x>y){
        stop("Check the Data format to Conditions! The maximum number of genes
            must be greater than the minimum number.")
    }
    return(c(x,y))
}

#x = SeedNumber
checkSeedNumber <- function(x){
    if(!is.numeric(x)){
        stop("Please, check the format of Conditions! SeedNumber argument must
            be numeric.")
    }else if(x<0){
         stop("Check the Data format Conditions! The seed must be a positive
            value.")
    }
    return(x)
}

#x = BootstrapNumber
checkBootstrapNumber <- function(x){
    if(!is.integer(x)){
        stop("Please, check the format of Conditions! BootstrapNumber argument
            must be integer.")
    }else if(x<=0){
        stop("Check the Data format Conditions! The number of bootstrap steps
            genes must be positive and different from zero.")
    }
    return(x)
}

#x = PCorrection
checkPCorrection <- function(x){
    if(!is.numeric(x)){
        stop("Please, check the format of Conditions! PCorrection argument
             must be numeric.")
    }else if(x<0 || x>1){
        stop("Check the Data format Conditions! Correction value must be 
            between zero and one.")
    }
    return(x)
}

#x = PCorrectionMethod
checkPCorrectionMethod <- function(x){
    options <- c("holm", "hochberg", "hommel", "bonferroni", "bh", "by","fdr")
    OptionsPAdjust <- c("holm", "hochberg", "hommel", "bonferroni", "BH",
                        "BY","fdr")
    if  (!(tolower(x) %in% options)){
        stop(paste("Correction method",x,"not found."))
    }
    return(OptionsPAdjust[grep(tolower(x),options)])
}

#x = WilcoxonTest
checkWilcoxonTest <- function(x){
    if(!is.logical(x)){
        stop("Please, check the format of Conditions! WilcoxonTest argument 
            must be logical.")
    }
    return(x)
}

#x = FisherTest
checkFisherTest <- function(x){
    if(!is.logical(x)){
        stop("Please, check the format of Conditions! FisherTest argument must
             be logical.")
    }
    return(x)
}

#x = domain (gobp, gocc, gomf, kegg or own) => AnalysisDomain
#y = DB (Own DB or OrgDB) => DBSpecies
#z = Gene identification (symbol or entrez) => GeneIdentifier
#k = database of gene expression => ExpressionData
checkAnalysisDomain_DBSpecies_GeneIdentifier <- function(x,y,z,k){
    options <- c("gobp","gocc","gomf","kegg","own")

    GeneNomenclatures <- c("symbol","entrez","tair","orf")

    SpeciesList <- list(c("Anopheles gambiae","aga"),
                        c("Arabdopsis thaliana","ath"),
                        c("Bos taurus","bta"),
                        c("Caenorhabditis elegans","cel"),
                        c("Canis familiaris","cfa"),
                        c("Drosophila melanogaster","dme"),
                        c("Danio rerio","dre"),
                        c("Escherichia coli K12","eco"),
                        c("Escherichia coli Sakai","ecs"),
                        c("Gallus gallus","gga"),
                        c("Homo sapiens","hsa"),
                        c("Mus musculus","mmu"),
                        c("Macaca mulatta","mcc"),
                        c("Plasmodium falciparum","pfa"),
                        c("Pan troglodytes","ptr"),
                        c("Rattus norvegicus","rno"),
                        c("Saccharomyces cerevisiae","sce"),
                        c("Sus scrofa","ssc"),
                        c("Xenopus laevis","xla"))

    DBSpeciesList <- c("org.Ag.eg.db", #Anopheles gambiae
                    "org.At.tair.db", #Arabdopsis thaliana
                    "org.Bt.eg.db", #Bos taurus
                    "org.Ce.eg.db", #Caenorhabditis elegans
                    "org.Cf.eg.db", #Canis familiaris
                    "org.Dm.eg.db", #Drosophila melanogaster
                    "org.Dr.eg.db", #Danio rerio (zebrafish)
                    "org.EcK12.eg.db", #Escherichia coli (strain K12)
                    "org.EcSakai.eg.db", #Escherichia coli (strain Sakai)
                    "org.Gg.eg.db", #Gallus gallus
                    "org.Hs.eg.db", #Homo sapiens
                    "org.Mm.eg.db", #Mus musculus (mouse)
                    "org.Mmu.eg.db", #Macaca mulatta (rhesus monkey)
                    "org.Pf.plasmo.db", #Plasmodium falciparum (malaria)
                    "org.Pt.eg.db", #Pan troglodytes (chimp)
                    "org.Rn.eg.db", #Rattus norvegicus (rattus)
                    "org.Sc.sgd.db", #Saccharomyces cerevisiae
                    "org.Ss.eg.db", #Sus scrofa (pig)
                    "org.Xl.eg.db") #Xenopus laevis (African Clawed Frog)

    SpeciesID <- y
    if(!is.data.frame(SpeciesID)){
        if(!file.exists(SpeciesID)){
            if(x!=options[5]){
                if(tolower(SpeciesID)=="eck12"){
                    SpeciesID <- "EcK12"
                }else if(tolower(SpeciesID)=="ecsakai"){
                    SpeciesID <- "EcSakai"
                }else{
                    SpeciesID <- unlist(strsplit(SpeciesID,""))
                    SpeciesID <- paste0(toupper(SpeciesID[1]),
                                paste0(tolower(SpeciesID[2:length(SpeciesID)]),
                                collapse = ""))
                }
                SpeciesID <- ifelse(length(DBSpeciesList[grep(SpeciesID,
                DBSpeciesList)])==0,y,
                DBSpeciesList[grep(SpeciesID,DBSpeciesList)])
            }
        }
    }

    ExpressionData <- k
    if(!(is.data.frame(k))){
        ExpressionData <- read.table(k, sep = "\t", header = TRUE, quote = "",
                                    stringsAsFactors = FALSE)
    }

    if(!(tolower(x) %in% options)){
        stop(paste("Analysis domain",x,"does not exist!"))
    }else if(tolower(x)==options[5]){
            if(!is.data.frame(SpeciesID)){
                if(file.exists(SpeciesID)){
                    if(!file.exists(SpeciesID)){
                        stop("Own functions file does not exist!")
                    }else{
                        DBFile <- read.table(SpeciesID, header = TRUE,
                                            sep = "\t", quote = "",
                                            stringsAsFactors = FALSE)
                        DBFile <- checkDB(DBFile = DBFile,
                                        ExpressionData = ExpressionData)
                    }
                }else{
                    stop("Own functions file does not exist!")
                }
        }else{
            DBFile <- checkDB(DBFile=SpeciesID, ExpressionData = ExpressionData)
        }
    }else if(is.data.frame(SpeciesID)){
        stop("DBSpecies does not exist!")
    }else if(!(SpeciesID %in% DBSpeciesList)){
        stop("DBSpecies does not exist!")
    }else if(checkPackage(NamePackage = SpeciesID)){
            if(!z %in% GeneNomenclatures){
                stop("Gene identifier invalid!")
            }else if(z==GeneNomenclatures[3] & SpeciesID!="org.At.tair.db"){
                stop("Incorrect GeneIdentifier!")
            }else if(z==GeneNomenclatures[4] & SpeciesID!="org.Sc.sgd.db" & 
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
                    SymbolEGClass <- org.Sc.sgdENTREZID #EntrezID / ORF numbers
                    SymbolEGList <- unlist(as.list(SymbolEGClass[mappedkeys(
                                    SymbolEGClass)]))
                    SymbolEGList <-data.frame(orf=names(SymbolEGList),
                                    entrez=as.vector(SymbolEGList))
            }else if(SpeciesID=="org.Pf.plasmo.db"){
                    SymbolEGClass <- org.Pf.plasmoALIAS2ORF #GeneSymbol / ORF ID
                    SymbolEGList <- unlist(as.list(SymbolEGClass[mappedkeys(
                                    SymbolEGClass)]))
                    SymbolEGList <-data.frame(symbol=names(SymbolEGList),
                                    orf=as.vector(SymbolEGList))
            }else{
                SymbolEGClass <- get(paste0("org.",unlist(strsplit(SpeciesID,
                                "[.]"))[2],".egSYMBOL2EG")) #GeneSymbol / Entrez
                SymbolEGList <- unlist(as.list(SymbolEGClass[mappedkeys(
                                SymbolEGClass)]))
                SymbolEGList <-data.frame(symbol=names(SymbolEGList),
                                entrez=as.vector(SymbolEGList))
            }

            rm(SymbolEGClass)

            if(sum(ifelse(ExpressionData[,c(1)] %in% 
                            SymbolEGList[,c(grep(tolower(z),
                            colnames(SymbolEGList)))] |
                            tolower(ExpressionData[,c(1)]) %in% 
                            SymbolEGList[,c(grep(tolower(z),
                            colnames(SymbolEGList)))] |
                            toupper(ExpressionData[,c(1)]) %in% 
                            SymbolEGList[,c(grep(tolower(z),
                            colnames(SymbolEGList)))],1,0)) < 1){
                stop("The genes in the Expression Data are not in the species 
                    database. Check the gene identifiers!")
            }else{
                DBFile <- checkGenesDomain(x = tolower(x), y = SpeciesID, 
                IDGene = tolower(z), 
                ListGenes = SymbolEGList,
                keggID = SpeciesList[[grep(SpeciesID,
                DBSpeciesList)]],
                ExpressionData = ExpressionData)
            }
        }
        rm(SymbolEGList)
    }

    rm(DBSpeciesList)
    rm(ExpressionData)

    return(list(options[grep(tolower(x),options)],DBFile[[1]],DBFile[[2]],
            tolower(z)))
}

checkPackage <- function(NamePackage){
    if(!is.element(NamePackage,installed.packages())){
        warning(paste("Installing DBSpecies",NamePackage))
        BiocManager::install(NamePackage)
    }
    if(!suppressWarnings(require(NamePackage, character.only = TRUE, 
        quietly = TRUE))){
        stop(paste("Package",NamePackage,"not installed. Check it and try
        again!"))
    }
    return(TRUE)
}

checkDB <- function(DBFile,ExpressionData){
    if(ncol(DBFile)!=3){
        stop("Own functions file must have 3 columns!")
    }else if(nrow(DBFile)<1){
        stop("Own functions file is empty!")
    }else if(sum(ifelse(ExpressionData[,c(1)] %in% DBFile[,c(1)] | 
        tolower(ExpressionData[,c(1)]) %in% DBFile[,c(1)] |
        toupper(ExpressionData[,c(1)]) %in% DBFile[,c(1)],1,0)) < 1){
        stop("The genes in the Expression Data are not in the species database.
            Check the gene identifiers!")
        }else if(length(colnames(DBFile))==length(colnames(ExpressionData))){
                if(sum(colnames(DBFile)==colnames(ExpressionData))==3){
                    stop("Check the ExpressionData and DBSpecies arguments.
                        The column names of both files are the same!!")
                }
    }

    DBFunctionsList <- lapply(as.vector(unique(as.character(DBFile[,2]))), 
                        FUN = ListConversionDB, DBFile = DBFile)
    names(DBFunctionsList) <- as.vector(unique(as.character(DBFile[,2])))
    TermOntologies <- unique(DBFile[,c(2,3)])
    colnames(TermOntologies) <- c("ID","Description")

    DBFunctionsListSample <- lapply(DBFunctionsList, function(ListElement,
                                    SampleGenes) {SampleGenes[SampleGenes %in% 
                                    ListElement]},
                                    SampleGenes = ExpressionData[,1])
    names(DBFunctionsListSample) <- sapply(names(DBFunctionsList),
                                    function(ListElement,TermOntologies){
                                    a<-subset(TermOntologies,
                                            TermOntologies[,1] == ListElement)
                                    return(paste0(as.character(a[,1])," <==> ",
                                            as.character(a[,2])))},
                                    TermOntologies=TermOntologies)

    DBFunctionsListRaw <- lapply(DBFunctionsList,function(ListElement){
                                length(ListElement)})
    names(DBFunctionsListRaw) <- names(DBFunctionsList)

    rm(DBFunctionsList)
    rm(TermOntologies)

    return(list(DBFunctionsListSample,DBFunctionsListRaw))
}

ListConversionDB <- function(ListElement,DBFile){
    return(as.character(subset(DBFile, DBFile[,2] == ListElement)[,1]))
}

checkGenesDomain <- function(x,y,IDGene,ListGenes,keggID,ExpressionData){
    ListDBPackages <- list(org.At.tair.db = c("org.At.tairGO2ALLTAIRS","tair",
                                            "entrez","org.At.tairARACYC "),
                        org.Sc.sgd.db = c("org.Sc.sgdGO2ALLORFS","orf",
                                        "entrez","org.Sc.sgdPATH"),
                        org.Pf.plasmo.db = c("org.Pf.plasmoGO2ALLORFS","orf",
                                        "symbol","org.Pf.plasmoPATH"),
                        otherDB = c(paste0("org.",unlist(strsplit(y,"[.]"))[2],
                                    ".egGO2ALLEGS"),"entrez","symbol",
                                    paste0("org.", unlist(strsplit(y,"[.]"))[2],
                                    ".egPATH2EG")))

    if(y %in% names(ListDBPackages)){
        OntologyFunction <- c(ListDBPackages[[grep(y,names(
                            ListDBPackages))]][1],
                            ListDBPackages[[grep(y,names(
                            ListDBPackages))]][2],
                            ListDBPackages[[grep(y,names(
                            ListDBPackages))]][3],
                            ListDBPackages[[grep(y,names(
                            ListDBPackages))]][4])
    }else{
        OntologyFunction <- c(ListDBPackages[[4]][1],ListDBPackages[[4]][2],
        ListDBPackages[[4]][3],ListDBPackages[[4]][4])
    }

    if(x == "kegg"){
        keggPathways <- as.data.frame(keggList("pathway", keggID[2]))
        keggPathwaysID <- unlist(lapply(rownames(keggPathways), function(b)
        unlist(strsplit(b,paste0(":",keggID[2])))[2]))
        
        keggPathwaysDescription <- unlist(lapply(as.character(keggPathways[,1]),
        function(b) unlist(strsplit(b,paste0("- ",keggID[1])))[1]))
        
        TermOntologies <- data.frame(ID = keggPathwaysID, 
                                    Description = keggPathwaysDescription)
        rm(keggPathways)
        rm(keggPathwaysID)
        rm(keggPathwaysDescription)

        DBKeggClass <- get(OntologyFunction[4])
        DBFunctionsList <- as.list(DBKeggClass[mappedkeys(DBKeggClass)])

        DBFunctionsList <- lapply(DBFunctionsList, function(k) 
                                unique(as.vector(k)))
        DBFunctionsList <- mapply(CheckTerm, ListElement = DBFunctionsList, 
                                ListElementName = names(DBFunctionsList),
                                MoreArgs =  list(VectorTerms = as.character(
                                TermOntologies$ID)))

        TermOntologies <- TermOntologies[TermOntologies$ID %in% 
                            names(DBFunctionsList),]
        TermOntologies$ID <- paste0(keggID[2],TermOntologies$ID)

        DBFunctionsList <- DBFunctionsList[lapply(DBFunctionsList, length) > 0]
        names(DBFunctionsList) <- paste0(keggID[2],names(DBFunctionsList))
        
        
    }else{
        infoGO <- GOTERM
        TermOntologies <- data.frame(t(sapply(infoGO, function(k) 
                    c(trimws(GOID(k)), trimws(Ontology(k)), trimws(Term(k))))))
        TermOntologies <- subset(TermOntologies, TermOntologies$X2 == 
                            toupper(unlist(strsplit(x,"go"))[2]))
        TermOntologies <- TermOntologies[,c(1,3)]
                            colnames(TermOntologies) <- c("ID","Description")
        rownames(TermOntologies) <- NULL

        DBFunctions <- get(OntologyFunction[1])
        DBFunctionsList <- as.list(DBFunctions[mappedkeys(DBFunctions)])
        DBFunctionsList <- lapply(DBFunctionsList, function(k)
                            unique(as.vector(k)))
        DBFunctionsList <- mapply(CheckTerm, ListElement = DBFunctionsList,
                            ListElementName = names(DBFunctionsList),
                            MoreArgs =  list(VectorTerms = 
                            as.character(TermOntologies$ID)))
        DBFunctionsList <- DBFunctionsList[lapply(DBFunctionsList, length) > 0]
        rm(infoGO)
        rm(DBFunctions)
    }

    if(IDGene != OntologyFunction[2]){
        DBFunctionsList <- lapply(DBFunctionsList, FUN = GeneIDConversion, 
                            ListGenes = ListGenes,
                            StandardGeneID = OntologyFunction[2], 
                            AlternativeGeneID = OntologyFunction[3])
    }

    rm(ListDBPackages)
    rm(OntologyFunction)

    TermOntologies <- subset(TermOntologies, TermOntologies$ID %in%
                        names(DBFunctionsList))

    DBFunctionsListSample <- lapply(DBFunctionsList, 
                            function(ListElement,SampleGenes){
                            SampleGenes[SampleGenes %in% 
                            ListElement]},
                            SampleGenes = ExpressionData[,1])
    
    names(DBFunctionsListSample) <- sapply(names(DBFunctionsList),
                                    function(ListElement,TermOntologies){
                                    a<-subset(TermOntologies, 
                                              TermOntologies[,1] == ListElement)
                                    return(paste0(as.character(a[,1])," <==> ",
                                                  as.character(a[,2])))},
                                    TermOntologies=TermOntologies)

    DBFunctionsListRaw <- lapply(DBFunctionsList,function(ListElement){
                                length(ListElement)})
    names(DBFunctionsListRaw) <- names(DBFunctionsList)

    DBFunctionsListSample <- DBFunctionsListSample[lapply(DBFunctionsListSample,
                                                          length) > 0]
    
    rm(DBFunctionsList)
    rm(TermOntologies)

    return(list(DBFunctionsListSample,DBFunctionsListRaw))
}

CheckTerm <- function(ListElement,ListElementName,VectorTerms){
    if(as.character(ListElementName) %in% as.character(VectorTerms)){
        result <- ListElement
    }else{
        result <- NULL
    }
    return(result)
}

GeneIDConversion <- function(ListElement,ListGenes,StandardGeneID,
                            AlternativeGeneID){
    z <- subset(ListGenes, ListGenes[,grep(StandardGeneID,colnames(ListGenes))]
                %in% ListElement)
    return(as.character(z[,grep(AlternativeGeneID,colnames(z))]))
}

