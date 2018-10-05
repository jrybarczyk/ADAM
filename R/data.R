#' @title Available species annotation packages
#' @description Contains a data.frame describing all available species
#' annotation packages from Bioconductor.
#' @format A data.frame with 19 rows and 3 variables
#' @name AvailableDBSpecies
#' @references GENTLEMAN, Robert C. et al. Bioconductor: open software 
#' development for computational biology and bioinformatics. Genome biology,
#' v. 5, n. 10, p. R80, 2004.
#' @examples 
#' data(AvailableDBSpecies)
NULL

#' @title Aedes aegypti RNA-seq differential expression data example
#' @description  A sample fragment of differential gene expression from an
#' RNA-seq experiment of Aedes aegypti mosquito.
#' @format A data.frame with 1963 rows and 3 variables
#' @name DiffAedes
#' @references AKBARI, O. S. et al. The developmental transcriptome of the 
#' mosquito aedes aegypti, an invasive species and major arbovirus vector. 
#' G3: Genes— Genomes— Genetics, Genetics Society of America, v. 3, n. 9, 
#' p. 1493–1509, 2013.
#' @examples 
#' data(DiffAedes)
NULL

#' @title Homo sapiens microarray differential expression data example
#' @description  A sample fragment of differential gene expression from a
#' microarray experiment of Homo sapiens.
#' @format A data.frame with 1964 rows and 3 variables
#' @name DiffHs
#' @references UEHARA, T. et al. The japanese toxicogenomics project: 
#' application of toxicogenomics. Molecular nutrition & food research, 
#' Wiley Online Library, v. 54, n. 2, p. 218–227, 2010.
#' @examples 
#' data(DiffHs)
NULL

#' @title Aedes aegypti RNA-seq data expression
#' @description  A sample fragment of gene expression from an RNA-seq experiment
#' of Aedes aegypti mosquito.
#' @format A data.frame with 2000 rows and 5 variables
#' @name ExpressionAedes
#' @references AKBARI, O. S. et al. The developmental transcriptome of the 
#' mosquito aedes aegypti, an invasive species and major arbovirus vector. G3:
#' Genes— Genomes— Genetics, Genetics Society of America, v. 3, n. 9, 
#' p. 1493–1509, 2013.
#' #' @examples 
#' data(ExpressionAedes)
NULL

#' @title Homo sapiens microarray data expression
#' @description  A sample fragment of gene expression from an RNA-seq experiment
#' of Homo sapiens.
#' @format A data.frame with 2000 rows and 5 variables
#' @name ExpressionHs
#' @references UEHARA, T. et al. The japanese toxicogenomics project: 
#' application of toxicogenomics. Molecular nutrition & food research, 
#' Wiley Online Library, v. 54, n. 2, p. 218–227, 2010.
#' @examples 
#' data(ExpressionHs)
NULL

#' @title Relation between Aedes aegypti genes and KEGG pathways
#' @description A relation between the GFAGs present in ResultAnalysisAedes and
#' their respective genes.
#' @format A data.frame with 782 rows and 2 variables
#' @name GeneFunctionAedes
#' @references Molan, A. L. 2018. “Construction of a Tool for Multispecies Genic
#' Functional Enrichment Analysis Among Comparative Samples.” Master’s thesis, 
#' Institute of Biosciences of Botucatu – Univ. Estadual Paulista.
#' http://hdl.handle.net/11449/157105.
#' @examples 
#' data(GeneFunctionAedes)
NULL

#' @title Relation between Homo sapiens genes and Cellular Components
#' @description A relation between the GFAGs present in ResultAnalysisHs and
#' their respective genes.
#' @format A data.frame with 43292 rows and 2 variables
#' @name GeneFunctionHs
#' @references Seco, G. B. S. 2018. "Comparative analysis of methods of 
#' determination of altered biological pathways in toxicogenomic data of in 
#' vitro and in vivo studies." Master’s thesis, Institute of Biosciences of
#' Botucatu – Univ. Estadual Paulista. http://hdl.handle.net/11449/157078
#' @examples 
#' data(GeneFunctionHs)
NULL

#' @title Relation between Aedes aegypti genes and KEGG pathways as 
#' ADAM input
#' @description A relation between the genes in the ExpressionAedes data and
#' their respective KEGG pathways (GFAGs).
#' @format A data.frame with 782 rows and 2 variables
#' @name KeggPathwaysAedes
#' @references Molan, A. L. 2018. “Construction of a Tool for Multispecies Genic
#' Functional Enrichment Analysis Among Comparative Samples.” Master’s thesis, 
#' Institute of Biosciences of Botucatu – Univ. Estadual Paulista.
#' http://hdl.handle.net/11449/157105.
#' @examples
#' data(KeggPathwaysAedes)
NULL

#' @title Result from an example of ADAM Aedes aegypti analysis
#' @description Result from ADAM Aedes aegypti analysis according to 
#' ExpressionAedes data.
#' @format A data.frame with 87 rows and 22 variables
#' @name ResultAnalysisAedes
#' @references Molan, A. L. 2018. “Construction of a Tool for Multispecies Genic
#' Functional Enrichment Analysis Among Comparative Samples.” Master’s thesis, 
#' Institute of Biosciences of Botucatu – Univ. Estadual Paulista.
#' http://hdl.handle.net/11449/157105.
#' @examples
#' data(ResultAnalysisAedes)
NULL

#' @title Result from an example of ADAM Homo sapiens analysis
#' @description Result from ADAM Homo sapiens analysis according to 
#' ExpressionHs data.
#' @format A data.frame with 554 rows and 22 variables
#' @name ResultAnalysisHs
#' @references Molan, A. L. 2018. “Construction of a Tool for Multispecies Genic
#' Functional Enrichment Analysis Among Comparative Samples.” Master’s thesis, 
#' Institute of Biosciences of Botucatu – Univ. Estadual Paulista.
#' http://hdl.handle.net/11449/157105.
#' @examples 
#' data(ResultAnalysisHs)
NULL
