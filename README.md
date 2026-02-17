![Bioconductor availability](https://bioconductor.org/shields/availability/release/ADAM.svg)
![Bioconductor downloads](https://bioconductor.org/shields/downloads/release/ADAM.svg)
![Years in Bioconductor](https://bioconductor.org/shields/years-in-bioc/ADAM.svg)
![Bioconductor build status](https://bioconductor.org/shields/build/release/bioc/ADAM.svg)
![Bioconductor dependencies](https://bioconductor.org/shields/dependencies/release/ADAM.svg)


# ADAM: Activity and Diversity Analysis Module

**Authors:**  André L. Molan, Giordano B. S. Seco, Agnes A. S. Takeda, Jose L. Rybarczyk-Filho  
**Maintainer:** José L. Rybarczyk-Filho

## Installation
`ADAM` is available on Bioconductor.

To install this package, start R and enter:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ADAM")
```

### Alternative installation method using devtools
If you prefer, you can also install `ADAM` using `devtools`. This method might be useful if you want to install the development version directly from a **GitHub repository** or another source. To install using `devtools`, you'll first need to ensure that `devtools` is installed. If it's not, you can install it using the following command:

```r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("jrybarczyk/ADAM")

```

## Quick start

```r
library(ADAM)
data(ExpressionAedes)
data(KeggPathwaysAedes)

result <- ADAnalysis(
  ComparisonID = c("control1,experiment1"),
  ExpressionData = ExpressionAedes,
  MinGene = 3L,
  MaxGene = 20L,
  DBSpecies = KeggPathwaysAedes,
  AnalysisDomain = "own",
  GeneIdentifier = "geneStableID"
)
```
