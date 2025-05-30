#' Euro Area Macroeconomic Data from Banbura and Modugno 2014
#'
#' A data extract from BM 2014 replication files. Some proprietary series (mostly PMI's) are excluded. The dataset \code{BM14_Models} provides information about all series
#' and their inclusion in the 'small', 'medium' and 'large' sized dynamic factor models estimated by BM 2014. The actual data is contained in \emph{xts} format in \code{BM14_M} for monthly data and \code{BM14_Q} for quarterly data.
#'
#' @aliases BM14_Models BM14_M BM14_Q
#' @usage
#' BM14_Models
#' BM14_M
#' BM14_Q
#'
#' @format \code{BM14_Models} is a data frame with 101 obs. (series) and 8 columns:
#' \describe{
#'   \item{series}{BM14 series code (converted to snake case for R)}
#'   \item{label}{BM14 series label}
#'   \item{code}{original series code from data source}
#'   \item{freq}{series frequency}
#'   \item{log_trans}{logical indicating whether the series was transformed by the natural log before differencing. Note that all data are provided in untransformed levels, and all data was (log-)differenced by BM14 before estimation.}
#'   \item{small}{logical indicating series included in the 'small' model of BM14. Proprietary series are excluded.}
#'   \item{medium}{logical indicating series included in the 'medium' model of BM14. Proprietary series are excluded.}
#'   \item{large}{logical indicating series included in the 'large' model of BM14. This comprises all series, thus the variable is redundant but included for completeness. Proprietary series are excluded.}
#' }
#' @source
#' Banbura, M., & Modugno, M. (2014). Maximum likelihood estimation of factor models on datasets with arbitrary pattern of missing data. \emph{Journal of Applied Econometrics, 29}(1), 133-160.
#'
#' @seealso \link{dfms-package}
#'
#' @examples
#' library(magrittr)
#' library(xts)
#'
#' # Constructing the database for the large model
#' BM14 = merge(BM14_M, BM14_Q)
#' BM14[, BM14_Models$log_trans] %<>% log()
#' BM14[, BM14_Models$freq == "M"] %<>% diff()
#' BM14[, BM14_Models$freq == "Q"] %<>% diff(3)
#'
#' # Small Model Database
#' head(BM14[, BM14_Models$small])
#'
#' # Medium-Sized Model Database
#' head(BM14[, BM14_Models$medium])
#'
"BM14_Models"
