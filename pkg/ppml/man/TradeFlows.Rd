\name{TradeFlows}

\alias{TradeFlows}

\title{TradeFlows}

\description{
Data on international trade in business services retrieved from the \emph{World Input-Output Database} from 2014. (WIOD, Timmer et al. (2015)). The definition of business services corresponds to the \emph{NACE Rev. 2} classifications of J62, J63, M and N. Trade flows are aggregated at the country level: The value of a trade flow from an industry \emph{k} in country \emph{i} to country \emph{j} is computed as the sum of the intermediate trade flows from \emph{k} to all importing industries in \emph{j} and the final consumption of the exported commodity in country \emph{j}. The subset of business services, where industries J62, J63, M and N are totalled, is then merged with the following datasets:


\emph{GeoDist (CEPII)} (Mayer and Zignago, 2011): GeoDist delivers a comprehensive set of gravity variables for pairs of countries. In addition to bilateral distance, the following dummy variables are added: \emph{contiguity, common official language, colony} and \emph{same country}

\emph{STRI - Service Trade Restrictiveness Index (OECD)} (Grosso et al., 2015): The STRI is a measure for service trade restrictiveness that assigns a score between zero and one to each importing country. The index reflects the regulatory non-tariff barriers to market access for foreign service providers. Zero represents an open market, one a market completely closed to foreign service providers.


\emph{Mario Larch's Regional Trade Agreements Database} (Egger and Larch, 2008): This comprehensive database delivers dummy variables to indicate if a country was part in a bilateral or multilateral RTA between 1950 and 2017. Based on all RTAs notified to the WTO, it distinguishes between free trade agreements, custom unions, economic integration agreements and partial scope agreements. The dummy variable used in this dataset is equal to one if one or more of those agreements were in force in 2014.
}


\usage{data("TradeFlows", package = "ppml")}

\format{
A data frame containing three variables with 213 observations.
\describe{
  \item{s}{Trade flows from country i to country j in relation to world production/expenditure Y}
  \item{log_dist}{The log of distance}
  \item{contig}{Dummy for contiguity}
  \item{comlang_off}{Dummy for common official language}
  \item{colony}{Dummy for colony}
  \item{smctry}{Dummy for same country (1 if country i and j have been the same country once)}
  \item{i_j}{Dummy variable equal to 1 for i != j}
  \item{log_stri_i_j}{The log of (1+STRI_j)*i_j}
  \item{rta}{Dummy variable equal to 1 if both countries had a trade agreement in 2014}
  \item{i}{Exporting country}
  \item{j}{Importing country}
}
}

\source{
Data provided by \url{http://www.wiod.org/database/wiots16}, \url{https://www.ewf.uni-bayreuth.de/en/research/RTA-data/index.html}, \url{http://www.cepii.fr/CEPII/en/bdd_modele/presentation.asp?id=6} and \url{https://stats.oecd.org/Index.aspx?DataSetCode=STRI}.
}

\references{

  Egger, P. and Larch, M. (2008). \dQuote{Interdependent preferential trade agreement memberships: An empirical analysis}, \emph{Journal of International Economics}, \bold{76}(2), 384--399

  Grosso, M. G., Gonzales, F., Miroudot, S., Nordas, H. K., Rouzet, D. and Ueno, A. (2015). \dQuote{Services trade restrictiveness index (STRI): Scoring and weighting methodology}, (177)
  
  Kleinheinz, L. (2019). \dQuote{General equilibrium effects of trade liberalization in business services}
  
  Mayer, T. and Zignago, S. (2011), \dQuote{Notes on cepiis distances measures: The geodist database}, \emph{Working papers}, CEPII
  
  Timmer, M. P., Dietzenbacher, E., Los, B., Stehrer, R. and de Vries, G. J. (2015). \dQuote{An illustrated user guide to the world input output database: the case of global automotive production},
  \emph{Review of International Economics}, \bold{23}(3), 575--605. \doi{10.1111/roie.12178}

}

\examples{
data("TradeFlows", package = "ppml")

}

\keyword{datasets}
