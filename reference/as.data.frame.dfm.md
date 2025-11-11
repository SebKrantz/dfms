# Extract Factor Estimates in a Data Frame

Extract Factor Estimates in a Data Frame

## Usage

``` r
# S3 method for class 'dfm'
as.data.frame(
  x,
  ...,
  method = "all",
  pivot = c("long", "wide.factor", "wide.method", "wide", "t.wide"),
  time = seq_row(x$F_pca),
  stringsAsFactors = TRUE
)
```

## Arguments

- x:

  an object class 'dfm'.

- ...:

  not used.

- method:

  character. The factor estimates to use: any of `"qml"`, `"2s"`,
  `"pca"` (multiple can be supplied) or `"all"` for all estimates.

- pivot:

  character. The orientation of the frame: `"long"`, `"wide.factor"` or
  `"wide.method"`, `"wide"` or `"t.wide"`.

- time:

  a vector identifying the time dimension, or `NULL` to omit a time
  variable.

- stringsAsFactors:

  make factors from method and factor identifiers. Same as option to
  [`as.data.frame.table`](https://rdrr.io/r/base/table.html).

## Value

A data frame of factor estimates.

## See also

[dfms-package](https://sebkrantz.github.io/dfms/reference/dfms-package.md)

## Examples

``` r
# \donttest{
library(xts)
# Fit DFM with 3 factors and 3 lags in the transition equation
mod <- DFM(diff(BM14_M), r = 3, p = 3)
#> Converged after 26 iterations.

# Taking a single estimate:
print(head(as.data.frame(mod, method = "qml")))
#>   Method Factor Time       Value
#> 1    QML     f1    1   3.0487579
#> 2    QML     f1    2  -0.4236117
#> 3    QML     f1    3  -9.4779810
#> 4    QML     f1    4 -10.8844746
#> 5    QML     f1    5  -6.0722144
#> 6    QML     f1    6  -1.0811737
print(head(as.data.frame(mod, method = "qml", pivot = "wide")))
#>   Time          f1        f2         f3
#> 1    1   3.0487579 -4.371719 -0.1134302
#> 2    2  -0.4236117 -2.190223 -6.3464355
#> 3    3  -9.4779810  3.824410 -2.1678243
#> 4    4 -10.8844746  6.942244  0.4569135
#> 5    5  -6.0722144  1.303051  0.7505150
#> 6    6  -1.0811737 -2.430000  1.5735378

# Adding a proper time variable
time <- index(BM14_M)[-1L]
print(head(as.data.frame(mod, method = "qml", time = time)))
#>   Method Factor       Time       Value
#> 1    QML     f1 1980-02-29   3.0487579
#> 2    QML     f1 1980-03-31  -0.4236117
#> 3    QML     f1 1980-04-30  -9.4779810
#> 4    QML     f1 1980-05-31 -10.8844746
#> 5    QML     f1 1980-06-30  -6.0722144
#> 6    QML     f1 1980-07-31  -1.0811737

# All estimates: different pivoting methods
for (pv in c("long", "wide.factor", "wide.method", "wide", "t.wide")) {
   cat("\npivot = ", pv, "\n")
   print(head(as.data.frame(mod, pivot = pv, time = time), 3))
}
#> 
#> pivot =  long 
#>   Method Factor       Time      Value
#> 1    PCA     f1 1980-02-29  0.8445713
#> 2    PCA     f1 1980-03-31  0.5259228
#> 3    PCA     f1 1980-04-30 -1.2107116
#> 
#> pivot =  wide.factor 
#>   Method       Time         f1         f2         f3
#> 1    PCA 1980-02-29  0.8445713 -0.7908231 -1.0289352
#> 2    PCA 1980-03-31  0.5259228 -0.6706157 -3.2251023
#> 3    PCA 1980-04-30 -1.2107116  0.0519631  0.9270935
#> 
#> pivot =  wide.method 
#>   Factor       Time        PCA    TwoStep        QML
#> 1     f1 1980-02-29  0.8445713  0.2903274  3.0487579
#> 2     f1 1980-03-31  0.5259228 -1.1656341 -0.4236117
#> 3     f1 1980-04-30 -1.2107116 -4.9014535 -9.4779810
#> 
#> pivot =  wide 
#>         Time     f1_PCA     f2_PCA     f3_PCA f1_TwoStep  f2_TwoStep f3_TwoStep
#> 1 1980-02-29  0.8445713 -0.7908231 -1.0289352  0.2903274 -1.26492938 -1.6769534
#> 2 1980-03-31  0.5259228 -0.6706157 -3.2251023 -1.1656341 -0.70548999 -5.5720725
#> 3 1980-04-30 -1.2107116  0.0519631  0.9270935 -4.9014535  0.06938226 -0.3158611
#>       f1_QML    f2_QML     f3_QML
#> 1  3.0487579 -4.371719 -0.1134302
#> 2 -0.4236117 -2.190223 -6.3464355
#> 3 -9.4779810  3.824410 -2.1678243
#> 
#> pivot =  t.wide 
#>         Time     f1_PCA f1_TwoStep     f1_QML     f2_PCA  f2_TwoStep    f2_QML
#> 1 1980-02-29  0.8445713  0.2903274  3.0487579 -0.7908231 -1.26492938 -4.371719
#> 2 1980-03-31  0.5259228 -1.1656341 -0.4236117 -0.6706157 -0.70548999 -2.190223
#> 3 1980-04-30 -1.2107116 -4.9014535 -9.4779810  0.0519631  0.06938226  3.824410
#>       f3_PCA f3_TwoStep     f3_QML
#> 1 -1.0289352 -1.6769534 -0.1134302
#> 2 -3.2251023 -5.5720725 -6.3464355
#> 3  0.9270935 -0.3158611 -2.1678243
# }
```
