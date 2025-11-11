# Remove and Impute Missing Values in a Multivariate Time Series

This function imputes missing values in a stationary multivariate time
series using various methods, and removes cases with too many missing
values.

## Usage

``` r
tsnarmimp(
  X,
  max.missing = 0.8,
  na.rm.method = c("LE", "all"),
  na.impute = c("median.ma.spline", "median.ma", "median", "rnorm"),
  ma.terms = 3L
)
```

## Arguments

- X:

  a `T x n` numeric data matrix (incl. ts or xts objects) or data frame
  of stationary time series.

- max.missing:

  numeric. Proportion of series missing for a case to be considered
  missing.

- na.rm.method:

  character. Method to apply concerning missing cases selected through
  `max.missing`: `"LE"` only removes cases at the beginning or end of
  the sample, whereas `"all"` always removes missing cases.

- na.impute:

  character. Method to impute missing values for the PCA estimates used
  to initialize the EM algorithm. Note that data are standardized
  (scaled and centered) beforehand. Available options are:

  |                      |     |                                                                                                                                                                                                                                       |     |
  |----------------------|-----|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----|
  | `"median"`           |     | simple series-wise median imputation.                                                                                                                                                                                                 |     |
  | `"rnorm"`            |     | imputation with random numbers drawn from a standard normal distribution.                                                                                                                                                             |     |
  | `"median.ma"`        |     | values are initially imputed with the median, but then a moving average is applied to smooth the estimates.                                                                                                                           |     |
  | `"median.ma.spline"` |     | "internal" missing values (not at the beginning or end of the sample) are imputed using a cubic spline, whereas missing values at the beginning and end are imputed with the median of the series and smoothed with a moving average. |     |

- ma.terms:

  the order of the (2-sided) moving average applied in `na.impute`
  methods `"median.ma"` and `"median.ma.spline"`.

## Value

The imputed matrix `X_imp`, with attributes:

- `"missing"`:

  a missingness matrix `W` matching the dimensions of `X_imp`.

- `"rm.rows"`:

  and a vector of indices of rows (cases) with too many missing values
  that were removed.

## See also

[dfms-package](https://sebkrantz.github.io/dfms/reference/dfms-package.md)

## Examples

``` r
library(xts)
str(tsnarmimp(BM14_M))
#> An xts object on 1980-01-31 / 2009-09-30 containing: 
#>   Data:    double [357, 92]
#>   Columns: ip_total, ip_tot_cstr, ip_tot_cstr_en, ip_constr, ip_im_goods ... with 87 more columns
#>   Index:   Date [357] (TZ: "UTC")
#>   xts Attributes:
#>     $ missing: logi [1:357, 1:92] TRUE TRUE TRUE TRUE TRUE TRUE ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:92] "ip_total" "ip_tot_cstr" "ip_tot_cstr_en" "ip_constr" ...
```
