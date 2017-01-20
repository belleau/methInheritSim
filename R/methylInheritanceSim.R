#' methylInheritanceSim: Simulation TODO
#'
#' This package does a simulation of multigeneration of bisulfite data
#'
#' @docType package
#'
#' @name methylInheritanceSim-package
#'
#' @aliases methylInheritanceSim-package methylPed
#'
#' @author Pascal Belleau,
#' Astrid Deschênes and 
#' Arnaud Droit
#'
#' Maintainer:
#' Pascal Belleau <pascal_belleau@hotmail.com>
#'
#' @seealso
#' \itemize{
#' \item \code{\link{runSim}} { TODO }
#' }
#'
#' @keywords package
NULL


#' All samples information, formated by \code{methylKit}, in a
#' \code{methylBase} format (for demo purpose).
#'
#' The object is a \code{methylBase}.
#' There is 12 samples (6 controls and 6 cases). Each
#' sample information is stored in a \code{methylRaw} object.
#' 
#' This dataset can be
#' used to test the \code{runSim} function.
#' 
#' @name samplesForChrSynthetic
#'
#' @docType data
#'
#' @aliases samplesForChrSynthetic
#'
#' @format A \code{methylBase} object contains the information for one 
#' generation. Each sample information is
#' stored in a \code{methylRaw} object. There is \code{methylRaw} objects
#'  (6 controls and 6 cases).
#'
#' @return A \code{methylBase} contains the information for one generation.
#' Each sample information is stored in a \code{methylRaw} object. 
#' There is \code{methylRaw} objects
#' (6 controls and 6 cases).
#'
#' @seealso
#' \itemize{
#' \item \code{\link{runSim}} {for running a
#' simulation analysis using methylKit info entry}
#' }
#'
#' @usage data(samplesForChrSynthetic)
#'
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading dataset
#' data(samplesForChrSynthetic)
#'
#' ## Run a permutation analysis
#' \dontrun{runSim(TODO)}
#'
NULL

