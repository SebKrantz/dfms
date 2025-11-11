# (Fast) Barebones Vector-Autoregression

Quickly estimate a VAR(p) model using Armadillo's inverse function.

## Usage

``` r
.VAR(x, p = 1L)
```

## Arguments

- x:

  data numeric matrix with time series in columns - without missing
  values.

- p:

  positive integer. The lag order of the VAR.

## Value

A list containing matrices `Y = x[-(1:p), ]`, `X` which contains lags
1 - p of `x` combined column-wise, `A` which is the \\np \times n\\
transition matrix, where n is the number of series in `x`, and the VAR
residual matrix `res = Y - X %*% A`.

A list with the following elements:

- `Y`:

  `x[-(1:p), ]`.

- `X`:

  lags 1 - p of `x` combined column-wise.

- `A`:

  \\np \times n\\ transition matrix, where n is the number of series in
  `x`.

- `res`:

  VAR residual matrix: `Y - X %*% A`.

## See also

[dfms-package](https://sebkrantz.github.io/dfms/reference/dfms-package.md)

## Examples

``` r
var = .VAR(diff(EuStockMarkets), 3)
str(var)
#> List of 4
#>  $ Y  : num [1:1856, 1:4] -2.88 -7.55 20.14 9.42 -4.7 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:4] "DAX" "SMI" "CAC" "FTSE"
#>  $ X  : num [1:1856, 1:12] 14.53 -2.88 -7.55 20.14 9.42 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:12] "DAX" "SMI" "CAC" "FTSE" ...
#>  $ A  : num [1:12, 1:4] 0.00946 -0.09901 0.04544 0.09475 0.02456 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:12] "DAX" "SMI" "CAC" "FTSE" ...
#>   .. ..$ : chr [1:4] "DAX" "SMI" "CAC" "FTSE"
#>  $ res: num [1:1856, 1:4] -3.39 -5.54 21.51 5.04 -4.91 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:4] "DAX" "SMI" "CAC" "FTSE"
var$A
#>               DAX          SMI         CAC         FTSE
#> DAX   0.009458998  0.027283428 -0.01093944  0.025490277
#> SMI  -0.099010841 -0.034164357 -0.06651577 -0.094882142
#> CAC   0.045444935  0.037822109  0.03972168 -0.018984419
#> FTSE  0.094745867  0.143551499  0.09179713  0.196154434
#> DAX   0.024559851  0.001705899  0.01593365  0.023913757
#> SMI  -0.021282229  0.018452931 -0.02701831 -0.029297417
#> CAC   0.063422809  0.109479561  0.08450745  0.026931128
#> FTSE -0.076819232 -0.097909530 -0.06820239 -0.022626747
#> DAX  -0.053821302 -0.124261943 -0.01894164 -0.040102476
#> SMI   0.001254346 -0.018010725  0.01927024 -0.001231464
#> CAC   0.049037401  0.100660797 -0.06168014  0.066430851
#> FTSE  0.037908210  0.122676727  0.01831834 -0.004359654
rm(var)
```
