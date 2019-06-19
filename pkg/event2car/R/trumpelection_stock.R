#' The Trump election and stock returns of 18 US tech firms.
#'
#' Datset containing stock returns of 18 important US tech firms (See https://finance.yahoo.com/u/yahoo-finance/watchlists/tech-stocks-that-move-the-market/),
#' the NASDAQ return and the election date of president Trump.
#' The dataset covers 356 days before the election (e.g. 2015-11-09)
#' and 365 days after the election (e.g. 2017-11-08). This results in 503 trading days.
#'
#' @docType data
#'
#' @usage data(trumpelection_stock)
#'
#' @keywords datasets
#'
#' @format Two objects of class \code{zoo} and one object of class \code{Date}.
#' \describe{
#'   \item{returns_firms}{stock return of 18 US tech firms from 2015-11-09 to 2017-11-08.}
#'   \item{return_indx}{NASDAQ return  from 2015-11-09 to 2017-11-08.}
#'   \item{trumpelection}{Date of Donald Trump Election.}
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

