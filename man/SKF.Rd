% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/my_RcppExports.R
\name{SKF}
\alias{SKF}
\title{(Fast) Stationary Kalman Filter}
\usage{
SKF(X, A, C, Q, R, F_0, P_0, loglik = FALSE)
}
\arguments{
\item{X}{numeric data matrix (\eqn{T \times n}{T x n}).}

\item{A}{transition matrix (\eqn{rp \times rp}{rp x rp}).}

\item{C}{observation matrix (\eqn{n \times rp}{n x rp}).}

\item{Q}{state covariance (\eqn{rp \times rp}{rp x rp}).}

\item{R}{observation covariance (\eqn{n \times n}{n x n}).}

\item{F_0}{initial state vector (\eqn{rp \times 1}{rp x 1}).}

\item{P_0}{initial state covariance (\eqn{rp \times rp}{rp x rp}).}

\item{loglik}{logical. Compute log-likelihood?}
}
\value{
Predicted and filtered state vectors and covariances.
\item{\code{F}}{\eqn{T \times rp}{T x rp} filtered state vectors. }
\item{\code{P}}{\eqn{rp \times rp \times T}{rp x rp x T} filtered state covariances. }
\item{\code{F_pred}}{\eqn{T \times rp}{T x rp} predicted state vectors. }
\item{\code{P_pred}}{\eqn{rp \times rp \times T}{rp x rp x T} predicted state covariances. }
\item{\code{loglik}}{value of the log likelihood. }
}
\description{
A simple and fast C++ implementation of the Kalman Filter for stationary data (or random walks - data should be mean zero and without a trend) with time-invariant system matrices and missing data.
}
\details{
The underlying state space model is:

\deqn{\textbf{x}_t = \textbf{C} \textbf{F}_t + \textbf{e}_t \sim N(\textbf{0}, \textbf{R})}{x(t) = C F(t) + e(t) ~ N(0, R)}
\deqn{\textbf{F}_t = \textbf{A F}_{t-1} + \textbf{u}_t \sim N(\textbf{0}, \textbf{Q})}{F(t) = A F(t-1) + u(t) ~ N(0, Q)}

where \eqn{x_t}{x(t)} is \code{X[t, ]}. The filter then first performs a time update (prediction)

\deqn{\textbf{F}_t = \textbf{A F}_{t-1}}{F(t) = A F(t-1)}
\deqn{\textbf{P}_t = \textbf{A P}_{t-1}\textbf{A}' + \textbf{Q}}{P(t) = A P(t-1) A' + Q}

where \eqn{P_t = Cov(F_t)}{P(t) = Cov(F(t))}. This is followed by the measurement update (filtering)

\deqn{\textbf{K}_t = \textbf{P}_t \textbf{C}' (\textbf{C P}_t \textbf{C}' + \textbf{R})^{-1}}{K(t) = P(t) C' inv(C P(t) C' + R)}
\deqn{\textbf{F}_t = \textbf{F}_t + \textbf{K}_t (\textbf{x}_t - \textbf{C F}_t)}{F(t) = F(t) + K(t) (x(t) - C F(t))}
\deqn{\textbf{P}_t = \textbf{P}_t - \textbf{K}_t\textbf{C P}_t}{P(t) = P(t) - K(t) C P(t)}

If a row of the data is all missing the measurement update is skipped i.e. the prediction becomes the filtered value. The log-likelihood is
computed as
\deqn{1/2 \sum_t \log(|St|)-e_t'S_te_t-n\log(2\pi)}{1/2 sum_t[log(det(S(t))) - e(t)' S(t) e(t) - n log(2 pi)]}
where \eqn{S_t = (C P_t C' + R)^{-1}}{S(t) = inv(C P(t) C' + R)} and \eqn{e_t = x_t - C F_t}{e(t) = x(t) - C F(t)} is the prediction error.

For further details see any textbook on time series such as Shumway & Stoffer (2017), which provide an analogous R implementation in \code{astsa::Kfilter0}.
For another fast (C-based) implementation that also allows time-varying system matrices and non-stationary data see \code{FKF::fkf}.
}
\examples{
# See ?SKFS

}
\references{
Shumway, R. H., & Stoffer, D. S. (2017). Time Series Analysis and Its Applications: With R Examples. Springer.

Harvey, A. C. (1990). Forecasting, structural time series models and the Kalman filter.

Hamilton, J. D. (1994). Time Series Analysis. Princeton university press.
}
\seealso{
\code{\link{FIS}} \code{\link{SKFS}} \link{dfms-package}
}
