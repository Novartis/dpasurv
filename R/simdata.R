#' This is a simulated data set in the (start,stop] format. The subject column
#' refers to subject ID and note that there are multiple rows per subject.
#' The columns x and dose correspond to a treatment and dose variable respectively, while
#' M refers to the longitudinal mediator values. The triplet (start, stop, event)
#' corresponds to the time-to-event data in the required (start,stop] format.
#'
#' @format A data frame with the following variables: \code{subject}, \code{x},
#' \code{dose} (factors), and \code{M}, \code{start}, \code{stop}, \code{event} (numeric).
#'
#' @name simdata
#' @docType data
#' @keywords data
NULL
