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
#' @srrstats {G1.3} *All statistical terminology should be clarified and unambiguously defined.*
#' @srrstats {G1.4} *Software should use [`roxygen2`](https://roxygen2.r-lib.org/) to document all functions.*
#' @srrstats {G1.4a} *All internal (non-exported) functions should also be documented in standard [`roxygen2`](https://roxygen2.r-lib.org/) format, along with a final `@noRd` tag to suppress automatic generation of `.Rd` files.*
#' @srrstats {G1.6} *Software should include code necessary to compare performance claims with alternative implementations in other R packages.*
#' -> Packages *nowcasting* and *nowcastDFM* are not on CRAN, are not C++ based, and provide more specialized functionality for nowcasting applications, thus althougth performance is a key objective to *dfms*, I don't provide benchmarks.
#' @srrstats {G2.0} *Implement assertions on lengths of inputs, particularly through asserting that inputs expected to be single- or multi-valued are indeed so.*
#' @srrstats {G2.0a} Provide explicit secondary documentation of any expectations on lengths of inputs
#' @srrstats {G2.1} *Implement assertions on types of inputs (see the initial point on nomenclature above).*
#' @srrstats {G2.1a} *Provide explicit secondary documentation of expectations on data types of all vector inputs.*
#' @srrstats {G2.2} *Appropriately prohibit or restrict submission of multivariate input to parameters expected to be univariate.*
#' @srrstats {G2.3} *For univariate character input:*
#' @srrstats {G2.3a} *Use `match.arg()` or equivalent where applicable to only permit expected values.*
#' @srrstats {G2.3b} *Either: use `tolower()` or equivalent to ensure input of character parameters is not case dependent; or explicitly document that parameters are strictly case-sensitive.*
#' @srrstats {G2.4} *Provide appropriate mechanisms to convert between different data types, potentially including:*
#' @srrstats {G2.4a} *explicit conversion to `integer` via `as.integer()`*
#' @srrstats {G2.4b} *explicit conversion to continuous via `as.numeric()`*
#' @srrstats {G2.4c} *explicit conversion to character via `as.character()` (and not `paste` or `paste0`)*
#' @srrstats {G2.4d} *explicit conversion to factor via `as.factor()`*
#' @srrstats {G2.4e} *explicit conversion from factor via `as...()` functions*
#' @srrstats {G2.5} *Where inputs are expected to be of `factor` type, secondary documentation should explicitly state whether these should be `ordered` or not, and those inputs should provide appropriate error or other routines to ensure inputs follow these expectations.*
#' @srrstats {G2.6} *Software which accepts one-dimensional input should ensure values are appropriately pre-processed regardless of class structures.*
#' @srrstats {G2.7} *Software should accept as input as many of the above standard tabular forms as possible, including extension to domain-specific forms.*
#' @srrstats {G2.8} *Software should provide appropriate conversion or dispatch routines as part of initial pre-processing to ensure that all other sub-functions of a package receive inputs of a single defined class or type.*
#' @srrstats {G2.9} *Software should issue diagnostic messages for type conversion in which information is lost (such as conversion of variables from factor to character; standardisation of variable names; or removal of meta-data such as those associated with [`sf`-format](https://r-spatial.github.io/sf/) data) or added (such as insertion of variable or column names where none were provided).*
#' @srrstats {G2.10} *Software should ensure that extraction or filtering of single columns from tabular inputs should not presume any particular default behaviour, and should ensure all column-extraction operations behave consistently regardless of the class of tabular data used as input.*
#' @srrstats {G2.11} *Software should ensure that `data.frame`-like tabular objects which have columns which do not themselves have standard class attributes (typically, `vector`) are appropriately processed, and do not error without reason. This behaviour should be tested. Again, columns created by the [`units` package](https://github.com/r-quantities/units/) provide a good test case.*
#' @srrstats {G2.12} *Software should ensure that `data.frame`-like tabular objects which have list columns should ensure that those columns are appropriately pre-processed either through being removed, converted to equivalent vector columns where appropriate, or some other appropriate treatment such as an informative error. This behaviour should be tested.*
#' -> In general, all input to DFM() is converted to matrix by collapse::qM(), which also removes classes from objects that are already matrix based. All internal code is based on plain matrices.
#' @srrstats {G2.13} *Statistical Software should implement appropriate checks for missing data as part of initial pre-processing prior to passing data to analytic algorithms.*
#' @srrstats {G2.14} *Where possible, all functions should provide options for users to specify how to handle missing (`NA`) data, with options minimally including:*
#' @srrstats {G2.14a} *error on missing data*
#' @srrstats {G2.14b} *ignore missing data with default warnings or messages issued*
#' @srrstats {G2.14c} *replace missing data with appropriately imputed values*
#' @srrstats {G2.15} *Functions should never assume non-missingness, and should never pass data with potential missing values to any base routines with default `na.rm = FALSE`-type parameters (such as [`mean()`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/mean.html), [`sd()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/sd.html) or [`cor()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html)).*
#' -> DFM() estimation assumes missing data, all routines are designed for missing data.
#' @srrstats {G2.16} *All functions should also provide options to handle undefined values (e.g., `NaN`, `Inf` and `-Inf`), including potentially ignoring or removing such values.*
#' # -> NaN will be treated like NA, and Inf and -Inf will be handled as any other numeric values (in line with nearly all other software, including base R). I initially call anyNA() to efficiently check for missingness, if FALSE, the software will assume complete cases.
#' @srrstats {G3.0} *Statistical software should never compare floating point numbers for equality. All numeric equality comparisons should either ensure that they are made between integers, or use appropriate tolerances for approximate equality.*
#' @srrstats {G3.1} *Statistical software which relies on covariance calculations should enable users to choose between different algorithms for calculating covariances, and should not rely solely on covariances from the `stats::cov` function.*
#' @srrstats {G3.1a} *The ability to use arbitrarily specified covariance methods should be documented (typically in examples or vignettes).*
#' -> All PCA routines I know of, including the authors original Matlab code, use Pearsons Covariance for the eigen decomposition, and I currently see no reason to provide alternatives in *dfms*.
#' @srrstats {G5.0} *Where applicable or practicable, tests should use standard data sets with known properties (for example, the [NIST Standard Reference Datasets](https://www.itl.nist.gov/div898/strd/), or data sets provided by other widely-used R packages).*
#' @srrstats {G5.1} *Data sets created within, and used to test, a package should be exported (or otherwise made generally available) so that users can confirm tests and run examples.*
#' @srrstats {G5.2} *Appropriate error and warning behaviour of all functions should be explicitly demonstrated through tests. In particular,*
#' @srrstats {G5.2a} *Every message produced within R code by `stop()`, `warning()`, `message()`, or equivalent should be unique*
#' @srrstats {G5.2b} *Explicit tests should demonstrate conditions which trigger every one of those messages, and should compare the result with expected values.*
#' @srrstats {G5.3} *For functions which are expected to return objects containing no missing (`NA`) or undefined (`NaN`, `Inf`) values, the absence of any such values in return objects should be explicitly tested.*
#' @srrstats {G5.4} **Correctness tests** *to test that statistical algorithms produce expected results to some fixed test data sets (potentially through comparisons using binding frameworks such as [RStata](https://github.com/lbraglia/RStata)).*
#' -> I have translated the authors original matlab code into R and run tests with that code (see misc/ directory in the repo). These comparisons yielded that my implementation is equivalent to the original Matlab code.
#'  testing can be improved upon in general, an idea would be to run tests against the r-translations of those Matlab codes in the testthat framework, but I have not come round to do that yet and I would first like to gather some more substantive feedback on the software.
#' @srrstats {G5.4a} *For new methods, it can be difficult to separate out correctness of the method from the correctness of the implementation, as there may not be reference for comparison. In this case, testing may be implemented against simple, trivial cases or against multiple implementations such as an initial R implementation compared with results from a C/C++ implementation.*
#' @srrstats {G5.4b} *For new implementations of existing methods, correctness tests should include tests against previous implementations. Such testing may explicitly call those implementations in testing, preferably from fixed-versions of other software, or use stored outputs from those where that is not possible.*
#' @srrstats {G5.4c} *Where applicable, stored values may be drawn from published paper outputs when applicable and where code from original implementations is not available*
#' @srrstats {G5.5} *Correctness tests should be run with a fixed random seed*
#' @srrstats {G5.6} **Parameter recovery tests** *to test that the implementation produce expected results given data with known properties. For instance, a linear regression algorithm should return expected coefficient values for a simulated data set generated from a linear model.*
#' @srrstats {G5.6a} *Parameter recovery tests should generally be expected to succeed within a defined tolerance rather than recovering exact values.*
#' @srrstats {G5.6b} *Parameter recovery tests should be run with multiple random seeds when either data simulation or the algorithm contains a random component. (When long-running, such tests may be part of an extended, rather than regular, test suite; see G4.10-4.12, below).*
#' @srrstats {G5.7} **Algorithm performance tests** *to test that implementation performs as expected as properties of data change. For instance, a test may show that parameters approach correct estimates within tolerance as data size increases, or that convergence times decrease for higher convergence thresholds.*
#' @srrstats {G5.8} **Edge condition tests** *to test that these conditions produce expected behaviour such as clear warnings or errors when confronted with data with extreme properties including but not limited to:*
#' @srrstats {G5.8a} *Zero-length data*
#' @srrstats {G5.8b} *Data of unsupported types (e.g., character or complex numbers in for functions designed only for numeric data)*
#' @srrstats {G5.8c} *Data with all-`NA` fields or columns or all identical fields or columns*
#' @srrstats {G5.8d} *Data outside the scope of the algorithm (for example, data with more fields (columns) than observations (rows) for some regression algorithms)*
#' @srrstats {G5.9} **Noise susceptibility tests** *Packages should test for expected stochastic behaviour, such as through the following conditions:*
#' @srrstats {G5.9a} *Adding trivial noise (for example, at the scale of `.Machine$double.eps`) to data does not meaningfully change results*
#' @srrstats {G5.9b} *Running under different random seeds or initial conditions does not meaningfully change results*
#' @srrstats {G5.10} *Extended tests should included and run under a common framework with other tests but be switched on by flags such as as a `<MYPKG>_EXTENDED_TESTS="true"` environment variable.* - The extended tests can be then run automatically by GitHub Actions for example by adding the following to the `env` section of the workflow:
#' @srrstats {G5.11} *Where extended tests require large data sets or other assets, these should be provided for downloading and fetched as part of the testing workflow.*
#' @srrstats {G5.11a} *When any downloads of additional data necessary for extended tests fail, the tests themselves should not fail, rather be skipped and implicitly succeed with an appropriate diagnostic message.*
#' @srrstats {G5.12} *Any conditions necessary to run extended tests such as platform requirements, memory, expected runtime, and artefacts produced that may need manual inspection, should be described in developer documentation such as a `CONTRIBUTING.md` or `tests/README.md` file.*
#' @srrstats {UL1.0} *Unsupervised Learning Software should explicitly document expected format (types or classes) for input data, including descriptions of types or classes which are not accepted; for example, specification that software accepts only numeric inputs in `vector` or `matrix` form, or that all inputs must be in `data.frame` form with both column and row names.*
#' @srrstats {UL1.1} *Unsupervised Learning Software should provide distinct sub-routines to assert that all input data is of the expected form, and issue informative error messages when incompatible data are submitted.*
#' @srrstats {UL1.2} *Unsupervised learning which uses row or column names to label output objects should assert that input data have non-default row or column names, and issue an informative message when these are not provided.*
#' @srrstats {UL1.3} *Unsupervised Learning Software should transfer all relevant aspects of input data, notably including row and column names, and potentially information from other `attributes()`, to corresponding aspects of return objects.*
#' @srrstats {UL1.3a} *Where otherwise relevant information is not transferred, this should be explicitly documented.*
#' @srrstats {UL1.4} *Unsupervised Learning Software should document any assumptions made with regard to input data; for example assumptions about distributional forms or locations (such as that data are centred or on approximately equivalent distributional scales). Implications of violations of these assumptions should be both documented and tested, in particular:*
#' @srrstats {UL1.4a} *Software which responds qualitatively differently to input data which has components on markedly different scales should explicitly document such differences, and implications of submitting such data.*
#' @srrstats {UL1.4b} *Examples or other documentation should not use `scale()` or equivalent transformations without explaining why scale is applied, and explicitly illustrating and contrasting the consequences of not applying such transformations.*
#' -> Data is always scaled and centered before estimation to remove the need for intercept terms in the Kalman Filter and Smoother. This is stated in the documentatio e.g. ?DFM
#' @srrstats {UL2.0} *Routines likely to give unreliable or irreproducible results in response to violations of assumptions regarding input data (see UL1.6) should implement pre-processing steps to diagnose potential violations, and issue appropriately informative messages, and/or include parameters to enable suitable transformations to be applied.*
#' @srrstats {UL2.1} *Unsupervised Learning Software should document any transformations applied to input data, for example conversion of label-values to `factor`, and should provide ways to explicitly avoid any default transformations (with error or warning conditions where appropriate).*
#' @srrstats {UL2.2} *Unsupervised Learning Software which accepts missing values in input data should implement explicit parameters controlling the processing of missing values, ideally distinguishing `NA` or `NaN` values from `Inf` values.*
#' @srrstats {UL2.3} *Unsupervised Learning Software should implement pre-processing routines to identify whether aspects of input data are perfectly collinear.*
#' -> It does not matter for DFMs if series are collinear. Users may deliberately choose to duplicate quarterly series to increase their weight in mixed-frequency estimations.
#' @srrstats {UL3.0} *Algorithms which apply sequential labels to input data (such as clustering or partitioning algorithms) should ensure that the sequence follows decreasing group sizes (so labels of "1", "a", or "A" describe the largest group, "2", "b", or "B" the second largest, and so on.)*
#' @srrstats {UL3.1} *Dimensionality reduction or equivalent algorithms which label dimensions should ensure that that sequences of labels follows decreasing "importance" (for example, eigenvalues or variance contributions).*
#' @srrstats {UL3.2} *Unsupervised Learning Software for which input data does not generally include labels (such as `array`-like data with no row names) should provide an additional parameter to enable cases to be labelled.*
#' @srrstats {UL4.0} *Unsupervised Learning Software should return some form of "model" object, generally through using or modifying existing class structures for model objects, or creating a new class of model objects.*
#' @srrstats {UL4.2} *The return object from Unsupervised Learning Software should include, or otherwise enable immediate extraction of, all parameters used to control the algorithm used.*
#' @srrstats {UL4.3} *Model objects returned by Unsupervised Learning Software should implement or appropriately extend a default `print` method which provides an on-screen summary of model (input) parameters and methods used to generate results. The `print` method may also summarise statistical aspects of the output data or results.*
#' @srrstats {UL4.3a} *The default `print` method should always ensure only a restricted number of rows of any result matrices or equivalent are printed to the screen.*
#' @srrstats {UL4.4} *Unsupervised Learning Software should also implement `summary` methods for model objects which should summarise the primary statistics used in generating the model (such as numbers of observations, parameters of methods applied). The `summary` method may also provide summary statistics from the resultant model.*
#' @srrstats {UL6.0} *Objects returned by Unsupervised Learning Software should have default `plot` methods, either through explicit implementation, extension of methods for existing model objects, through ensuring default methods work appropriately, or through explicit reference to helper packages such as [`factoextra`](https://github.com/kassambara/factoextra) and associated functions.*
#' @srrstats {UL6.2} *Where default plot methods include labelling components of return objects (such as cluster labels), routines should ensure that labels are automatically placed to ensure readability, and/or that appropriate diagnostic messages are issued where readability is likely to be compromised (for example, through attempting to place too many labels).*
#' @srrstats {UL7.0} *Inappropriate types of input data are rejected with expected error messages.*
#' @srrstats {UL7.1} *Tests should demonstrate that violations of assumed input properties yield unreliable or invalid outputs, and should clarify how such unreliability or invalidity is manifest through the properties of returned objects.*
#' @srrstats {UL7.2} *Demonstrate that labels placed on output data follow decreasing group sizes (**UL3.0**)*
#' @srrstats {UL7.3} *Demonstrate that labels on input data are propagated to, or may be recovered from, output data.
#' @srrstats {UL7.4} *Demonstrate that submission of new data to a previously fitted model can generate results more efficiently than initial model fitting.*
#' @srrstats {TS1.0} *Time Series Software should use and rely on explicit class systems developed for representing time series data, and should not permit generic, non-time-series input*
#' @srrstats {TS1.1} *Time Series Software should explicitly document the types and classes of input data able to be passed to each function.*
#' @srrstats {TS1.2} *Time Series Software should implement validation routines to confirm that inputs are of acceptable classes (or represented in otherwise appropriate ways for software which does not use class systems).*
#' @srrstats {TS1.3} *Time Series Software should implement a single pre-processing routine to validate input data, and to appropriately transform it to a single uniform type to be passed to all subsequent data-processing functions (the [`tsbox` package](https://www.tsbox.help/) provides one convenient approach for this).*
#' @srrstats {TS1.4} *The pre-processing function described above should maintain all time- or date-based components or attributes of input data.*
#' @srrstats {TS1.5} *The software should ensure strict ordering of the time, frequency, or equivalent ordering index variable.*
#' @srrstats {TS1.6} *Any violations of ordering should be caught in the pre-processing stages of all functions.*
#' -> This software is at the intersection of dimensionality reduction and time series, and does not require input objects to have a certain class.
#' For all practical purposes I believe it is convenient to not require certain object types, although certain input object attributes such as the
#' time scale of 'ts' or 'xts' objects could be used in model plots. This is an area where I still want to improve the software e.g. keeping it class
#' agnostic but using some of the information contained in time series objects passed as inputs.
#' @srrstats {TS1.7} *Accept inputs defined via the [`units` package](https://github.com/r-quantities/units/) for attributing SI units to R vectors.*
#' @srrstats {TS1.8} *Where time intervals or periods may be days or months, be explicit about the system used to represent such, particularly regarding whether a calendar system is used, or whether a year is presumed to have 365 days, 365.2422 days, or some other value.*
#' @srrstats {TS2.0} *Time Series Software which presumes or requires regular data should only allow **explicit** missing values, and should issue appropriate diagnostic messages, potentially including errors, in response to any **implicit** missing values.*
#' @srrstats {TS2.1} *Where possible, all functions should provide options for users to specify how to handle missing data, with options minimally including:*
#' @srrstats {TS2.1a} *error on missing data; or.
#' @srrstats {TS2.1b} *warn or ignore missing data, and proceed to analyse irregular data, ensuring that results from function calls with regular yet missing data return identical values to submitting equivalent irregular data with no missing values; or*
#' @srrstats {TS2.1c} *replace missing data with appropriately imputed values.*
#' @srrstats {TS2.2} *Consider stationarity of all relevant moments - typically first (mean) and second (variance) order, or otherwise document why such consideration may be restricted to lower orders only.*
#' @srrstats {TS2.3} *Explicitly document all assumptions and/or requirements of stationarity*
#' @srrstats {TS2.4} *Implement appropriate checks for all relevant forms of stationarity, and either:*
#' @srrstats {TS2.4a} *issue diagnostic messages or warnings; or*
#' @srrstats {TS2.4b} *enable or advise on appropriate transformations to ensure stationarity.*
#' -> it is possible to run DFM() on non-stationary data and obtain convergent results, although this is strongly discouraged.
#' @srrstats {TS2.5} *Incorporate a system to ensure that both row and column orders follow the same ordering as the underlying time series data. This may, for example, be done by including the `index` attribute of the time series data as an attribute of the auto-covariance matrix.*
#' @srrstats {TS2.6} *Where applicable, auto-covariance matrices should also include specification of appropriate units.*
#' @srrstats {TS3.0} *Provide tests to demonstrate at least one case in which errors widen appropriately with forecast horizon.*
#' @srrstats {TS3.1} *If possible, provide at least one test which violates TS3.0*
#' -> currently I don't forecast the covariance matrices. This could be implemented in the future.
#' @srrstats {TS3.2} *Document the general drivers of forecast errors or horizons, as demonstrated via the particular cases of TS3.0 and TS3.1*
#' @srrstats {TS3.3} *Either:*
#' @srrstats {TS3.3a} *Document, preferable via an example, how to trim forecast values based on a specified error margin or equivalent; or*
#' @srrstats {TS3.3b} *Provide an explicit mechanism to trim forecast values to a specified error margin, either via an explicit post-processing function, or via an input parameter to a primary analytic function.*
#' @srrstats {TS4.0} *Return values should either:*
#' @srrstats {TS4.0a} *Be in same class as input data, for example by using the [`tsbox` package](https://www.tsbox.help/) to re-convert from standard internal format (see 1.4, above); or*
#' @srrstats {TS4.0b} *Be in a unique, preferably class-defined, format.*
#' @srrstats {TS4.1} *Any units included as attributes of input data should also be included within return values.*
#' @srrstats {TS4.2} *The type and class of all return values should be explicitly documented.*
#' @srrstats {TS4.3} *Return values should explicitly include all appropriate units and/or time scales*
#' @srrstats {TS4.4} *Document the effect of any such transformations on forecast data, including potential effects on both first- and second-order estimates.*
#' @srrstats {TS4.5} *In decreasing order of preference, either:*
#' @srrstats {TS4.5a} *Provide explicit routines or options to back-transform data commensurate with original, non-stationary input data*
#' @srrstats {TS4.5b} *Demonstrate how data may be back-transformed to a form commensurate with original, non-stationary input data.*
#' @srrstats {TS4.5c} *Document associated limitations on forecast values*
#' @srrstats {TS4.6} *Time Series Software which implements or otherwise enables forecasting should return either:*
#' @srrstats {TS4.6a} *A distribution object, for example via one of the many packages described in the CRAN Task View on [Probability Distributions](https://cran.r-project.org/web/views/Distributions.html) (or the new [`distributional` package](https://pkg.mitchelloharawild.com/distributional/) as used in the [`fable` package](https://fable.tidyverts.org) for time-series forecasting).*
#' @srrstats {TS4.6b} *For each variable to be forecast, predicted values equivalent to first- and second-order moments (for example, mean and standard error values).*
#' @srrstats {TS4.6c} *Some more general indication of error associated with forecast estimates.*
#' @srrstats {TS4.7} *Ensure that forecast (modelled) values are clearly distinguished from observed (model or input) values, either (in this case in no order of preference) by*
#' @srrstats {TS4.7a} *Returning forecast values alone*
#' @srrstats {TS4.7b} *Returning distinct list items for model and forecast values*
#' @srrstats {TS4.7c} *Combining model and forecast values into a single return object with an appropriate additional column clearly distinguishing the two kinds of data.*
#' @srrstats {TS5.0} *Implement default `plot` methods for any implemented class system.*
#' @srrstats {TS5.1} *When representing results in temporal domain(s), ensure that one axis is clearly labelled "time" (or equivalent), with continuous units.*
#' @srrstats {TS5.2} *Default to placing the "time" (or equivalent) variable on the horizontal axis.*
#' @srrstats {TS5.3} *Ensure that units of the time, frequency, or index variable are printed by default on the axis.*
#' @srrstats {TS5.5} *Provide options to determine whether plots of data with missing values should generate continuous or broken lines.*
#' @srrstats {TS5.6} *By default indicate distributional limits of forecast on plot*
#' @srrstats {TS5.7} *By default include model (input) values in plot, as well as forecast (output) values*
#' @srrstats {TS5.8} *By default provide clear visual distinction between model (input) values and forecast (output) values.*
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
