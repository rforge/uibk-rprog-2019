#' The Trump election and stock returns of 18 US tech firms.
#'
#' Datset containing stock returns of 18 important US tech firms (See \url{https://finance.yahoo.com/u/yahoo-finance/watchlists/tech-stocks-that-move-the-market/}),
#' the NASDAQ return.
#' The dataset covers trading days between 2015-11-09 and 2017-11-08. This results in 503 trading days.
#'
#' @docType data
#'
#' @usage data("tech_returns", package="estudy2car")
#'
#' @keywords datasets
#'
#' @format An objects of class \code{zoo} containing 503 observations and 19 variables.
#' \describe{
#'   \item{^NDX}{NASDAQ return  from 2015-11-09 to 2017-11-08.}
#'   \item{MSFT}{Microsoft's from 2015-11-09 to 2017-11-08.}
#'   \item{AMZ}{Amazon's from 2015-11-09 to 2017-11-08.}
#'   \item{AAPL}{Apple's from 2015-11-09 to 2017-11-08.}
#'   \item{GOOG}{Google's from 2015-11-09 to 2017-11-08.}
#'   \item{FB}{Facebook's from 2015-11-09 to 2017-11-08.}
#'   \item{BABA}{Alibaba's from 2015-11-09 to 2017-11-08.}
#'   \item{INTC}{Intel's from 2015-11-09 to 2017-11-08.}
#'   \item{PYPL}{PayPal's from 2015-11-09 to 2017-11-08.}
#'   \item{NVDA}{NVIDIA's from 2015-11-09 to 2017-11-08.}
#'   \item{TSLA}{Tesla's from 2015-11-09 to 2017-11-08.}
#'   \item{ATVI}{Activision Blizzard's from 2015-11-09 to 2017-11-08.}
#'   \item{AMD}{Advanced Micro's from 2015-11-09 to 2017-11-08.}
#'   \item{EA}{Electronic Arts's from 2015-11-09 to 2017-11-08.}
#'   \item{MTCH}{Match Group's from 2015-11-09 to 2017-11-08.}
#'   \item{TTD}{The Trade Desk's from 2015-11-09 to 2017-11-08.}
#'   \item{ZG}{Zillow Group's from 2015-11-09 to 2017-11-08.}
#'   \item{YELP}{Yelp's from 2015-11-09 to 2017-11-08.}
#'   \item{TIVO}{TiVo's from 2015-11-09 to 2017-11-08.}
#'   ...
#' }
#' @source \url{https://finance.yahoo.com/}
#'
#' #' @examples
#' # load data
#' data('tech_returns')
#' # mean adjusted model
#' trumpelection <- as.Date("2016-11-08")
#' event2car(returns=returns_firms,regressor=return_indx,event_dates =trumpelection,market_model="mean_adj" )
#' # market adjusted model (out-of sample estimation)
#' event2car(returns=returns_firms,regressor=return_indx,event_dates =trumpelection,market_model="mrkt_adj" )
#' # market adjusted model (within sample estimation)
#' event2car(returns=returns_firms,regressor=return_indx,event_dates =trumpelection,market_model="mrkt_adj_within" )

"tech_returns"

