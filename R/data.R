# Data documentation for Roxygen

#' 1router data from Cao et al. (2000)
#'
#' Data from 4-node network with star topology collected from Bell Labs; used in
#' Cao et al. (2000).
#'
#' @section Objects:
#' This imports several objects:
#' \itemize{
#'  \item \code{onerouter}, a data.frame with all data
#'  \item \code{X}, a matrix of origin-destination flows formatted for analysis
#'  \item \code{Y}, a matrix of link loads formatted for analysis
#'  \item \code{tvec}, a vector of times
#'  \item \code{A}, the routing matrix for this network (truncated for full row
#'      rank)
#' }
#' In this data, we have \code{A \%*\% t(X) == t(Y)}.
#'
#' @section Variables:
#' The data.frame \code{onerouter} contains the following:
#' \itemize{
#'  \item value, level of traffic recorded
#'  \item nme, name of flow or load
#'  \item method, whether flow was directly observered or inferred
#'      (all observed)
#'  \item time, time of observation
#'  \item od, flag for origin-destination vs. link loads
#'  \item orig, origin of flow or load
#'  \item dest, destination of flow or load
#'  \item node, node involved in flow or load
#' }
#'
#' @docType data
#' @name onerouter
#' @usage onerouter
#' @references J. Cao, D. Davis, S. Van Der Viel, and B. Yu.
#' Time-varying network tomography: router link data.
#' Journal of the American Statistical Association, 95:1063-75, 2000.
#' @keywords datasets
#' @family onerouter
NULL

#' 1router OD data from Cao et al. (2000)
#'
#' Data from 4-node network with star topology collected from Bell Labs; used in
#' Cao et al. (2000). The X matrix contains only the OD flows, formatted for
#' analysis.
#'
#' @section Structure:
#' The columns of the X matrix correspond to individual OD flows, and the rows
#' correspond to observations.
#' In this data, we have \code{A \%*\% t(X) == t(Y)}.
#'
#' @docType data
#' @name X
#' @usage X
#' @references J. Cao, D. Davis, S. Van Der Viel, and B. Yu.
#' Time-varying network tomography: router link data.
#' Journal of the American Statistical Association, 95:1063-75, 2000.
#' @keywords datasets
#' @family onerouter
NULL

#' 1router link data from Cao et al. (2000)
#'
#' Data from 4-node network with star topology collected from Bell Labs; used in
#' Cao et al. (2000). The Y matrix contains only the link loads, formatted for
#' analysis. The last (8th) link load has been removed to ensure that the
#' resulting routing matrix is of full row rank.
#'
#' @section Structure:
#' The columns of the Y matrix correspond to individual link loads, and the rows
#' correspond to observations.
#' In this data, we have \code{A \%*\% t(X) == t(Y)}.
#'
#' @docType data
#' @name Y
#' @usage Y
#' @references J. Cao, D. Davis, S. Van Der Viel, and B. Yu.
#' Time-varying network tomography: router link data.
#' Journal of the American Statistical Association, 95:1063-75, 2000.
#' @keywords datasets
#' @family onerouter
NULL

#' 1router time data from Cao et al. (2000)
#'
#' Data from 4-node network with star topology collected from Bell Labs; used in
#' Cao et al. (2000). The numeric vector \code{tvec} contains the time in
#' decimal hours since midnight for each observation. 
#'
#' @docType data
#' @name tvec
#' @usage tvec
#' @references J. Cao, D. Davis, S. Van Der Viel, and B. Yu.
#' Time-varying network tomography: router link data.
#' Journal of the American Statistical Association, 95:1063-75, 2000.
#' @keywords datasets
#' @family onerouter
NULL

#' 1router routing matrix from Cao et al. (2000)
#'
#' Data from 4-node network with star topology collected from Bell Labs; used in
#' Cao et al. (2000). The A matrix contains the routing structure of the given
#' network. It is 7 x 16 due to the full row rank constraint.
#' In this data, we have \code{A \%*\% t(X) == t(Y)}.
#'
#' @section Structure:
#' The columns of this matrix correspond to individual OD flows (the columns of
#' X), and its rows correspond to individual link loads (the columns of Y).
#'
#' @docType data
#' @name A
#' @usage A
#' @references J. Cao, D. Davis, S. Van Der Viel, and B. Yu.
#' Time-varying network tomography: router link data.
#' Journal of the American Statistical Association, 95:1063-75, 2000.
#' @keywords datasets
#' @family onerouter
NULL
