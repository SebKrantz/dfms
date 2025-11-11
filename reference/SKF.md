# (Fast) Stationary Kalman Filter

A simple and fast C++ implementation of the Kalman Filter for stationary
data (or random walks - data should be mean zero and without a trend)
with time-invariant system matrices and missing data.

## Usage

``` r
SKF(X, A, C, Q, R, F_0, P_0, loglik = FALSE)
```

## Arguments

- X:

  numeric data matrix (\\T \times n\\).

- A:

  transition matrix (\\rp \times rp\\).

- C:

  observation matrix (\\n \times rp\\).

- Q:

  state covariance (\\rp \times rp\\).

- R:

  observation covariance (\\n \times n\\).

- F_0:

  initial state vector (\\rp \times 1\\).

- P_0:

  initial state covariance (\\rp \times rp\\).

- loglik:

  logical. Compute log-likelihood?

## Value

Predicted and filtered state vectors and covariances.

- `F`:

  \\T \times rp\\ filtered state vectors.

- `P`:

  \\rp \times rp \times T\\ filtered state covariances.

- `F_pred`:

  \\T \times rp\\ predicted state vectors.

- `P_pred`:

  \\rp \times rp \times T\\ predicted state covariances.

- `loglik`:

  value of the log likelihood.

## Details

The underlying state space model is:

\$\$\textbf{x}\_t = \textbf{C} \textbf{F}\_t + \textbf{e}\_t \sim
N(\textbf{0}, \textbf{R})\$\$ \$\$\textbf{F}\_t = \textbf{A F}\_{t-1} +
\textbf{u}\_t \sim N(\textbf{0}, \textbf{Q})\$\$

where \\x_t\\ is `X[t, ]`. The filter then first performs a time update
(prediction)

\$\$\textbf{F}\_t = \textbf{A F}\_{t-1}\$\$ \$\$\textbf{P}\_t =
\textbf{A P}\_{t-1}\textbf{A}' + \textbf{Q}\$\$

where \\P_t = Cov(F_t)\\. This is followed by the measurement update
(filtering)

\$\$\textbf{K}\_t = \textbf{P}\_t \textbf{C}' (\textbf{C P}\_t
\textbf{C}' + \textbf{R})^{-1}\$\$ \$\$\textbf{F}\_t = \textbf{F}\_t +
\textbf{K}\_t (\textbf{x}\_t - \textbf{C F}\_t)\$\$ \$\$\textbf{P}\_t =
\textbf{P}\_t - \textbf{K}\_t\textbf{C P}\_t\$\$

If a row of the data is all missing the measurement update is skipped
i.e. the prediction becomes the filtered value. The log-likelihood is
computed as \$\$1/2 \sum_t \log(\|St\|)-e_t'S_te_t-n\log(2\pi)\$\$ where
\\S_t = (C P_t C' + R)^{-1}\\ and \\e_t = x_t - C F_t\\ is the
prediction error.

For further details see any textbook on time series such as Shumway &
Stoffer (2017), which provide an analogous R implementation in
`astsa::Kfilter0`. For another fast (C-based) implementation that also
allows time-varying system matrices and non-stationary data see
`FKF::fkf`.

## References

Shumway, R. H., & Stoffer, D. S. (2017). Time Series Analysis and Its
Applications: With R Examples. Springer.

Harvey, A. C. (1990). Forecasting, structural time series models and the
Kalman filter.

Hamilton, J. D. (1994). Time Series Analysis. Princeton university
press.

## See also

[`FIS`](https://sebkrantz.github.io/dfms/reference/FIS.md)
[`SKFS`](https://sebkrantz.github.io/dfms/reference/SKFS.md)
[dfms-package](https://sebkrantz.github.io/dfms/reference/dfms-package.md)

## Examples

``` r
# See ?SKFS
```
