# (Fast) Fixed-Interval Smoother (Kalman Smoother)

(Fast) Fixed-Interval Smoother (Kalman Smoother)

## Usage

``` r
FIS(A, F, F_pred, P, P_pred, F_0 = NULL, P_0 = NULL)
```

## Arguments

- A:

  transition matrix (\\rp \times rp\\).

- F:

  state estimates (\\T \times rp\\).

- F_pred:

  state predicted estimates (\\T \times rp\\).

- P:

  variance estimates (\\rp \times rp \times T\\).

- P_pred:

  predicted variance estimates (\\rp \times rp \times T\\).

- F_0:

  initial state vector (\\rp \times 1\\) or empty (`NULL`).

- P_0:

  initial state covariance (\\rp \times rp\\) or empty (`NULL`).

## Value

Smoothed state and covariance estimates, including initial (t = 0)
values.

- `F_smooth`:

  \\T \times rp\\ smoothed state vectors, equal to the filtered state in
  period \\T\\.

- `P_smooth`:

  \\rp \times rp \times T\\ smoothed state covariance, equal to the
  filtered covariance in period \\T\\.

- `F_smooth_0`:

  \\1 \times rp\\ initial smoothed state vectors, based on `F_0`.

- `P_smooth_0`:

  \\rp \times rp\\ initial smoothed state covariance, based on `P_0`.

## Details

The Kalman Smoother is given by:

\$\$\textbf{J}\_t = \textbf{P}\_t \textbf{A} +
inv(\textbf{P}^{pred}\_{t+1})\$\$ \$\$\textbf{F}^{smooth}\_t =
\textbf{F}\_t + \textbf{J}\_t (\textbf{F}^{smooth}\_{t+1} -
\textbf{F}^{pred}\_{t+1})\$\$ \$\$\textbf{P}^{smooth}\_t =
\textbf{P}\_t + \textbf{J}\_t (\textbf{P}^{smooth}\_{t+1} -
\textbf{P}^{pred}\_{t+1}) \textbf{J}\_t'\$\$

The initial smoothed values for period t = T are set equal to the
filtered values. If `F_0` and `P_0` are supplied, the smoothed initial
conditions (t = 0 values) are also calculated and returned. For further
details see any textbook on time series such as Shumway & Stoffer
(2017), which provide an analogous R implementation in
`astsa::Ksmooth0`.

## References

Shumway, R. H., & Stoffer, D. S. (2017). Time Series Analysis and Its
Applications: With R Examples. Springer.

Harvey, A. C. (1990). Forecasting, structural time series models and the
Kalman filter.

## See also

[`SKF`](https://sebkrantz.github.io/dfms/reference/SKF.md)
[`SKFS`](https://sebkrantz.github.io/dfms/reference/SKFS.md)
[dfms-package](https://sebkrantz.github.io/dfms/reference/dfms-package.md)

## Examples

``` r
# See ?SKFS
```
