#' Simulated Cancer Genome Dataset
#'
#' An artificial cancer dataset describing mutations found in 9 different tumour samples.
#' Rows represent mutations.
#'
#'
#' @format ## `oncosim`
#' A data frame with 143 rows and 3 columns:
#' \describe{
#'   \item{Samples}{Sample containing mutation in the specified gene}
#'   \item{Genes}{Mutated  gene}
#'   \item{VariantType}{Type of mutation in gene}
#' }
#'
#' @source not applicable, simulated
"oncosim"


#' Simulated Cancer Dataset Metadata
#'
#' A sample‐level metadata table for the \code{oncosim} simulated cancer dataset.
#' Contains assorted numeric, categorical, clinical, and logical features for each sample.
#'
#' @format ## `oncosim_metadata`
#' A data frame with 11 rows and 6 columns:
#' \describe{
#'   \item{Samples}{Unique sample identifiers}
#'   \item{numeric_feature}{Numeric variable including zeros, positive and negative values, \code{NA}, and \code{Inf}/\code{-Inf}}
#'   \item{categorical_feature4levels}{Categorical variable with four levels (\code{"cat"}, \code{"dog"}, \code{"magpie"}, \code{"giraffe"}), may contain empty strings or \code{NA}}
#'   \item{clinical_feature2levels}{Clinical categorical variable indicating biological sex with two levels (\code{"male"}, \code{"female"}), may contain \code{NA}}
#'   \item{logical_feature}{Logical variable with \code{TRUE}, \code{FALSE}, or \code{NA}}
#'   \item{numeric_that_could_be_logical}{Integer variable coded as \code{0} or \code{1} (and \code{NA}) that could be interpreted as logical}
#' }
#'
#' @source not applicable, simulated
"oncosim_metadata"
