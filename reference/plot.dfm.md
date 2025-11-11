# Plot DFM

Plot DFM

## Usage

``` r
# S3 method for class 'dfm'
plot(
  x,
  method = switch(x$em.method, none = "2s", "qml"),
  type = c("joint", "individual", "residual"),
  scale.factors = TRUE,
  ...
)

# S3 method for class 'dfm'
screeplot(x, ...)
```

## Arguments

- x:

  an object class 'dfm'.

- method:

  character. The factor estimates to use: one of `"qml"`, `"2s"`,
  `"pca"` or `"all"` to plot all estimates.

- type:

  character. The type of plot: `"joint"`, `"individual"` or
  `"residual"`.

- scale.factors:

  logical. Standardize factor estimates, this usually improves the plot
  since the factor estimates corresponding to the greatest PCA
  eigenvalues tend to have a greater variance than the data.

- ...:

  for `plot.dfm`: further arguments to
  [`plot`](https://rdrr.io/r/base/plot.html),
  [`ts.plot`](https://rdrr.io/r/stats/ts.plot.html), or
  [`boxplot`](https://rdrr.io/r/graphics/boxplot.html), depending on the
  `type` of plot. For `screeplot.dfm`: further arguments to
  [`screeplot.ICr`](https://sebkrantz.github.io/dfms/reference/ICr.md).

## Value

Nothing.

## See also

[dfms-package](https://sebkrantz.github.io/dfms/reference/dfms-package.md)

## Examples

``` r
# \donttest{
# Fit DFM with 3 factors and 3 lags in the transition equation
mod <- DFM(diff(BM14_M), r = 3, p = 3)
#> Converged after 26 iterations.
plot(mod)

plot(mod, type = "individual", method = "all")

plot(mod, type = "residual")

# }
```
