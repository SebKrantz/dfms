# DFM Summary Methods

Summary and print methods for class 'dfm'. `print.dfm` just prints basic
model information and the factor transition matrix \\\textbf{A}\\,
`coef.dfm` returns \\\textbf{A}\\ and \\\textbf{C}\\ in a plain list,
whereas `summary.dfm` returns all system matrices and additional
residual and goodness of fit statistics—with a print method allowing
full or compact printout.

## Usage

``` r
# S3 method for class 'dfm'
print(x, digits = 4L, ...)

# S3 method for class 'dfm'
coef(object, ...)

# S3 method for class 'dfm'
logLik(object, ...)

# S3 method for class 'dfm'
summary(object, method = switch(object$em.method, none = "2s", "qml"), ...)

# S3 method for class 'dfm_summary'
print(x, digits = 4L, compact = sum(x$info["n"] > 15, x$info["n"] > 40), ...)
```

## Arguments

- x, object:

  an object class 'dfm'.

- digits:

  integer. The number of digits to print out.

- ...:

  not used.

- method:

  character. The factor estimates to use: one of `"qml"`, `"2s"` or
  `"pca"`.

- compact:

  integer. Display a more compact printout: `0` prints everything, `1`
  omits the observation matrix \\\textbf{C}\\ and residual covariance
  matrix `cov(resid(model))`, and `2` omits all disaggregated
  information on the input data. Sensible default are chosen for
  different sizes of the input dataset so as to limit large printouts.

## Value

Summary information following a dynamic factor model estimation.
[`coef()`](https://rdrr.io/r/stats/coef.html) returns \\\textbf{A}\\ and
\\\textbf{C}\\.

## See also

[dfms-package](https://sebkrantz.github.io/dfms/reference/dfms-package.md)

## Examples

``` r
mod <- DFM(diff(BM14_Q), 2, 3)
#> Converged after 26 iterations.
print(mod)
#> Dynamic Factor Model: n = 9, T = 117, r = 2, p = 3, %NA = 7.5973
#> 
#> Factor Transition Matrix [A]
#>     L1.f1  L1.f2  L2.f1   L2.f2   L3.f1   L3.f2
#> f1 0.6789 0.2413 -0.034 -0.4640 -0.0012 -0.1988
#> f2 0.0353 0.2270 -0.026  0.0645 -0.0744  0.1802
summary(mod)
#> Dynamic Factor Model: n = 9, T = 117, r = 2, p = 3, %NA = 7.5973
#> 
#> Call:  DFM(X = diff(BM14_Q), r = 2, p = 3)
#> 
#> Summary Statistics of Factors [F]
#>       N     Mean  Median      SD      Min     Max
#> f1  117  -0.0084  0.3469  2.2931  -14.408  3.7167
#> f2  117    0.003  0.0867  0.8146  -2.4636  2.1071
#> 
#> Factor Transition Matrix [A]
#>      L1.f1  L1.f2    L2.f1    L2.f2     L3.f1   L3.f2
#> f1 0.67890 0.2413 -0.03401 -0.46403 -0.001235 -0.1988
#> f2 0.03533 0.2270 -0.02598  0.06451 -0.074449  0.1802
#> 
#> Factor Covariance Matrix [cov(F)]
#>          f1       f2
#> f1  5.2584   0.1622 
#> f2  0.1622   0.6636 
#> 
#> Factor Transition Error Covariance Matrix [Q]
#>        u1     u2
#> u1 2.7065 0.2039
#> u2 0.2039 0.6618
#> 
#> Observation Matrix [C]
#>                 f1      f2
#> gdp         0.4094 -0.1237
#> priv_cons   0.2755 -0.4353
#> invest      0.3810 -0.3022
#> export      0.3842  0.4215
#> import      0.3911  0.2106
#> empl        0.3072 -0.3443
#> prductivity 0.2894  0.0222
#> capacity    0.2933  0.0157
#> gdp_us      0.2511  0.1259
#> 
#> Observation Error Covariance Matrix [diag(R) - Restricted]
#>         gdp   priv_cons      invest      export      import        empl 
#>      0.0953      0.4616      0.1685      0.0301      0.1232      0.4160 
#> prductivity    capacity      gdp_us 
#>      0.2317      0.4737      0.6360 
#> 
#> Observation Residual Covariance Matrix [cov(resid(DFM))]
#>                   gdp priv_cons    invest    export    import      empl
#> gdp           0.0670   -0.0046   -0.0096   -0.0038*  -0.0357*  -0.0679*
#> priv_cons    -0.0046    0.4173   -0.0735*   0.0059    0.0060   -0.0777*
#> invest       -0.0096   -0.0735*   0.1263   -0.0010   -0.0210   -0.0620*
#> export       -0.0038*   0.0059   -0.0010    0.0061   -0.0119*   0.0058 
#> import       -0.0357*   0.0060   -0.0210   -0.0119*   0.1090    0.0367 
#> empl         -0.0679*  -0.0777*  -0.0620*   0.0058    0.0367    0.3816 
#> prductivity   0.0683*   0.0273    0.0173   -0.0024   -0.0748*  -0.2119*
#> capacity     -0.0470*  -0.0978*  -0.0320   -0.0201*   0.0624*   0.0628 
#> gdp_us       -0.0209   -0.0023    0.0045   -0.0095   -0.0252   -0.0236 
#>             prductivity  capacity    gdp_us
#> gdp             0.0683*  -0.0470*  -0.0209 
#> priv_cons       0.0273   -0.0978*  -0.0023 
#> invest          0.0173   -0.0320    0.0045 
#> export         -0.0024   -0.0201*  -0.0095 
#> import         -0.0748*   0.0624*  -0.0252 
#> empl           -0.2119*   0.0628   -0.0236 
#> prductivity     0.2215   -0.1119*  -0.0496 
#> capacity       -0.1119*   0.4666   -0.0059 
#> gdp_us         -0.0496   -0.0059    0.6353 
#> 
#> Residual AR(1) Serial Correlation
#>         gdp   priv_cons      invest      export      import        empl 
#>   -0.149924   -0.272110   -0.206251   -0.215038   -0.007949    0.434361 
#> prductivity    capacity      gdp_us 
#>    0.053696    0.091852    0.179658 
#> 
#> Summary of Residual AR(1) Serial Correlations
#>   N     Mean   Median      SD      Min     Max
#>   9  -0.0102  -0.0079  0.2282  -0.2721  0.4344
#> 
#> Goodness of Fit: R-Squared
#>         gdp   priv_cons      invest      export      import        empl 
#>      0.9330      0.5827      0.8737      0.9939      0.8910      0.6184 
#> prductivity    capacity      gdp_us 
#>      0.7785      0.5334      0.3647 
#> 
#> Summary of Individual R-Squared's
#>   N    Mean  Median      SD     Min     Max
#>   9  0.7299  0.7785  0.2139  0.3647  0.9939
```
