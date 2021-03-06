% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/event2car.R
\name{event2car}
\alias{event2car}
\title{get_car_events: Calculates \code{CAR} (cumulative abnormal returns) for firm(s) and event(s)}
\usage{
event2car(returns = returns, regressor = regressor,
  event_dates = event_dates, estimation_period = 250, car_lag = 0,
  car_lead = 5, market_model = c("mean_adj", "mrkt_adj",
  "mrkt_adj_within"))
}
\arguments{
\item{returns}{an object of \code{data.frame} or \code{zoo} containing rates of returns of securities.}

\item{regressor}{an object of the same class as \code{returns} containing regressors.
The argument can be omitted, if market model is \code{mean_adj}.}

\item{event_dates}{an object of class \code{Date} containing one event date or multiple event dates}

\item{estimation_period}{an object of class \code{intenger} stating the number of days
prior to the event over which the market model parameters are estimated. Default is 250 days.
Note that the event period itself is not included in the event period to prevent the event from influencing the normal performance model parameter estimates.}

\item{car_lag}{an object of class \code{intenger} measuring the start of the event window. The default is 1 day prior to the event date.}

\item{car_lead}{an object of class \code{intenger} measuring the end of the event window. The default is 5 days after the event date.}

\item{market_model}{market_model a character indicating the market model among
\code{mean_adj}, \code{mrkt_adj}, and \code{sim}.}
}
\value{
an object of class \code{data.frame} which contains cumulative abnormal returns, the average cumulative abnormal return
(controls for varying event period durations if non-trading days are in the period), the number of tradingdays,
significance levels of the market-return-coefficient, significance level of the event-dummy-coefficient, and model fit (rquared) per securities per event date(s).
Note that significance levels and rsqared are NA if market model is \code{mean_adj}.
}
\description{
The function applies the following steps to securities' rates of returns: A)
calculate abnormal returns for \code{estimation_period}(s), B) predict
abnormal returns for event period(s), C) aggregate the predicted abnormal
returns to cumulative abnormal return(s) for \code{event_dates}.
The package covers three models for the calculation of the cumulative abnomral returns:
- Mean-Adjusted Model (\code{mean_adj})
- Market-Adjusted Model (\emph{mrkt_adj})
- Within-Sample Marekt-Adjusted Model (\emph{mrkt_adj_within})
This is the logic suggested by multiple scholars. See references below.
}
\details{
The generic function is dispatched for such classes as
 \code{zoo}. (future versions of the package allow for classes of \code{data.frame}.)
If \code{market_model} is \emph{mrkt_adj} or \code{sim}
and \code{regressor} has the length greater than one, the first element of
\code{regressor} will be applied for each security in \code{returns}.
}
\references{
MacKinlay, A.C. \emph{Event Studies in Economics and Finance}.
Journal of Economic Literature, 35(1):13-39, 1997.

Brown S.J., Warner J.B. \emph{Using Daily Stock Returns, The Case
of Event Studies}. Journal of Financial Economics, 14:3-31, 1985.

Davies, R., Studnicka, Z. \emph{The heterogeneous impact of Brexit: Early indications from the FTSE}.
European Economic Review, 110:1-17, 2018.
}
