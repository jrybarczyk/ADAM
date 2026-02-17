\name{run_adam}
\alias{run_adam}
\title{Run ADAM with Bioconductor Containers}
\usage{
run_adam(
  x,
  mode = c("complete", "partial"),
  assay_name = NULL,
  group_col = NULL,
  contrast = NULL,
  reference = NULL,
  comparison_id = NULL,
  gene_id_source = NULL,
  aggregate_by = NULL,
  aggregate_fun = c("sum", "mean"),
  MinGene = 3,
  MaxGene = 2000,
  SeedNumber = 1049,
  BootstrapNumber = 1000,
  PCorrection = 0.05,
  DBSpecies = NULL,
  PCorrectionMethod = "fdr",
  WilcoxonTest = FALSE,
  FisherTest = FALSE,
  AnalysisDomain = NULL,
  GeneIdentifier = NULL
)
}
\arguments{
\item{x}{Expression container: \code{SummarizedExperiment},
\code{SingleCellExperiment}, \code{ExpressionSet}, \code{data.frame}, or a
tab-delimited file path.}

\item{mode}{\code{"complete"} for \code{GFAGAnalysis} or \code{"partial"} for
\code{ADAnalysis}.}

\item{assay_name}{Assay name used for
\code{SummarizedExperiment}/\code{SingleCellExperiment}. If \code{NULL},
the first assay is used.}

\item{group_col}{Sample metadata column used to create comparisons when
\code{comparison_id} is not supplied.}

\item{contrast}{Group labels as \code{c(control, experiment)}.}

\item{reference}{Reference (control) group used to infer
\code{contrast} when there are two groups.}

\item{comparison_id}{Explicit ADAM comparison vector in format
\code{"sampleA,sampleB"}.}

\item{gene_id_source}{Optional rowData column used as gene identifier.}

\item{aggregate_by}{Optional metadata columns used for pseudo-bulk aggregation
before running ADAM (for example \code{c("sample_id", "cluster")}).}

\item{aggregate_fun}{Aggregation method: \code{"sum"} or \code{"mean"}.}

\item{MinGene, MaxGene, SeedNumber, BootstrapNumber, PCorrection, DBSpecies, PCorrectionMethod, WilcoxonTest, FisherTest, AnalysisDomain, GeneIdentifier}{Standard ADAM arguments forwarded to \code{GFAGAnalysis}/\code{ADAnalysis}.}
}
\value{
A \code{SummarizedExperiment} with one assay (\code{metrics}) and row-level
annotations for GFAG results. \code{metadata()} contains:
\code{legacy_result}, \code{gene_function_map}, \code{comparison_id},
and \code{mode}.
}
\description{
High-level Bioconductor-oriented wrapper for ADAM. It accepts common
Bioconductor expression containers, derives comparisons from sample metadata,
and returns results in \code{SummarizedExperiment} for easier downstream
integration.
}
\examples{
data(ExpressionAedes)
data(KeggPathwaysAedes)

se <- run_adam(
  x = ExpressionAedes,
  mode = "partial",
  comparison_id = "control1,experiment1",
  DBSpecies = KeggPathwaysAedes,
  AnalysisDomain = "own",
  GeneIdentifier = "geneStableID"
)
}
