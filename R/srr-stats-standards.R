#' srr_stats
#'
#' All of the following standards initially have `@srrstatsTODO` tags.
#' These may be moved at any time to any other locations in your code.
#' Once addressed, please modify the tag from `@srrstatsTODO` to `@srrstats`,
#' or `@srrstatsNA`, ensuring that references to every one of the following
#' standards remain somewhere within your code.
#' (These comments may be deleted at any time.)
#'
#' @srrstatsVerbose TRUE
#'
#' @srrstats {G1.0} *Statistical Software should list at least one primary reference from published academic literature.*
#' @srrstats {G1.1} *Statistical Software should document whether the algorithm(s) it implements are:* - *The first implementation of a novel algorithm*; or - *The first implementation within **R** of an algorithm which has previously been implemented in other languages or contexts*; or - *An improvement on other implementations of similar algorithms in **R***.
#' -> See README, DFM() implements simple baseline versions of algorithms that have been around for a while in Matlab, and in other langaues (R, Python, Julia), but inside more elaborate nowcasting codes - thus not directly accessible, and less efficient.
#' @srrstats {G1.2} *Statistical Software should include a* Life Cycle Statement *describing current and anticipated future states of development.*
#' -> See README: "The package is stable, but functionality may expand in the future. In particular, mixed-frequency estimation with autoregressive errors is planned for the near future, and generation of the 'news' may be added in the further future."
#' @srrstats {G1.3} *All statistical terminology should be clarified and unambiguously defined.*
#' @srrstats {G1.4} *Software should use [`roxygen2`](https://roxygen2.r-lib.org/) to document all functions.*
#' @srrstats {G1.4a} *All internal (non-exported) functions should also be documented in standard [`roxygen2`](https://roxygen2.r-lib.org/) format, along with a final `@noRd` tag to suppress automatic generation of `.Rd` files.*
#' @srrstats {G1.6} *Software should include code necessary to compare performance claims with alternative implementations in other R packages.*
#' -> Packages *nowcasting* and *nowcastDFM* are not on CRAN, are not C++ based, and provide more specialized functionality for nowcasting applications, thus althougth performance is a key objective to *dfms*, I don't provide benchmarks.
#' @noRd
NULL

#' NA_standards
#'
#' Any non-applicable standards can have their tags changed from `@srrstatsTODO`
#' to `@srrstatsNA`, and placed together in this block, along with explanations
#' for why each of these standards have been deemed not applicable.
#' (These comments may also be deleted at any time.)
#'
#' @srrstatsNA {G1.5} *Software should include all code necessary to reproduce results which form the basis of performance claims made in associated publications.*
#' -> No publications.
#' @srrstatsNA {G4.0} *Statistical Software which enables outputs to be written to local files should parse parameters specifying file names to ensure appropriate file suffices are automatically generated where not provided.*
#' @srrstatsNA {UL3.3} *Where applicable, Unsupervised Learning Software should implement routines to predict the properties (such as numerical ordinates, or cluster memberships) of additional new data without re-running the entire algorithm.*
#' -> DFMs are always used to reduce or forecast the data used to estimate them.
#' @srrstatsNA {UL3.4} *Objects returned from Unsupervised Learning Software which labels, categorise, or partitions data into discrete groups should include, or provide immediate access to, quantitative information on intra-group variances or equivalent, as well as on inter-group relationships where applicable.*
#' @srrstatsNA {UL4.1} *Unsupervised Learning Software may enable an ability to generate a model object without actually fitting values. This may be useful for controlling batch processing of computationally intensive fitting algorithms.*
#' @srrstatsNA {UL6.1} *Where the default `plot` method is **NOT** a generic `plot` method dispatched on the class of return objects (that is, through an S3-type `plot.<myclass>` function or equivalent), that method dispatch (or equivalent) should nevertheless exist in order to explicitly direct users to the appropriate function.*
#' @srrstatsNA {UL7.5} *Batch processing routines should be explicitly tested, commonly via extended tests (see **G4.10**--**G4.12**).*
#' @srrstatsNA {UL7.5a} *Tests of batch processing routines should demonstrate that equivalent results are obtained from direct (non-batch) processing.*
#' @srrstatsNA {TS5.4} *For frequency visualization, abscissa spanning $[-\pi, \pi]$ should be avoided in favour of positive units of $[0, 2\pi]$ or $[0, 0.5]$, in all cases with appropriate additional explanation of units.*
#' @noRd
NULL
