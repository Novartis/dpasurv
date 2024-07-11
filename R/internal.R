#' Find variables of a given formula
#'
#' Do not call this function on its own
#'
#' @param formula
#'
#' @return character vector of the formula's variable names
#' @export
#'
#' @keywords internal
find.variables <- function(formula) {

  `%>%` <- dplyr::`%>%`

  variables <- base::deparse(formula) %>%
    base::paste(collapse = "") %>%
    (function(x) {base::gsub("[[:space:]]", "", x)}) %>%
    base::strsplit("~") %>%
    base::unlist()

  return(variables)

}

#' Internal check of whether the formulas passed to the function dpa follow a directed acyclic graph (DAG)
#'
#' Do not call this function on its own
#'
#' @param meta meta data obtained internally inside the dpasurv::dpa function
#' @param data input data to the dpasurv::dpa function
#'
#' @return this function doesn't return anything
#' @export
#'
#' @keywords internal
check.dag <- function(meta, data) {

  `%>%` <- dplyr::`%>%`

  # Mediator k should not be a covariate of another mediator regression j > k
  if (base::length(meta$mediator$y) == 1) {
    is.dag <- TRUE
  } else {
    for (ii in 1:(base::length(meta$mediator$y)-1)) {
      upstream.covariates <- base::strsplit(base::paste(base::gsub("~","", meta$mediator$xreg[(ii+1):base::length(meta$mediator$y)]), collapse="+"), "+", fixed=TRUE) %>%
        base::unlist() %>%
        base::unique() %>%
        base::setdiff(base::c("+","~"))
      if (meta$mediator$y[ii] %in% upstream.covariates) {
        is.dag <- FALSE
        which.equation <- base::lapply(base::strsplit(base::gsub("~","", meta$mediator$xreg[(ii+1):base::length(meta$mediator$y)]), "+", fixed=TRUE), function(x) meta$mediator$y[ii] %in% x) %>% base::unlist() %>% base::which()
        stop(base::paste0("Mediator ", ii, " (", meta$mediator$y[ii], ") appears as covariate in mediator equation ", ii + which.equation ," > ", ii, " violating the DAG structure. Please refine DAG or order mediator equations differently."))
      }
    }
  }

  # Make sure all mediators are numeric:
  if(sum(base::apply(data[meta$mediator$y], 2, base::class)!="numeric") > 0)
    stop(base::paste0("mediators (", base::paste(meta$mediator$y, collapse=", "), ") should be numeric"))

}

#' Internal function that gathers meta data on the input to function dpa
#'
#' Do not call this function on its own
#'
#' @param out.formula obtained directly from corresponding input to dpasurv::dpa
#' @param mediator.formulas obtained directly from corresponding input to dpasurv::dpa
#' @param data obtained directly from corresponding input to dpasurv::dpa
#'
#' @return this function returns meta data associated with the call to dpasurv::dpa
#' @export
#'
#' @keywords internal
get.meta <- function(out.formula, mediator.formulas, data) {

  `%>%` <- dplyr::`%>%`

  # list to store all the individual components:
  comp <- base::list()

  # Parse the out.formula:
  out.xy = base::deparse(out.formula) %>% base::paste(collapse = "") %>%
    (function(x) {base::gsub("[[:space:]]", "", x)}) %>% base::strsplit("~") %>% base::unlist()

  comp[["outcome"]] <- base::list(startt = base::substr(out.xy[1], 5 + base::unlist(base::gregexpr("Surv\\(", out.xy[1])), base::unlist(base::gregexpr(",", out.xy[1]))[1]-1),
                            stopt = base::substr(out.xy[1], 1 + base::unlist(base::gregexpr(",", out.xy[1]))[1], base::unlist(base::gregexpr(",", out.xy[1]))[2]-1),
                            event = base::substr(out.xy[1], 1 + base::unlist(base::gregexpr(",", out.xy[1]))[2], base::unlist(base::gregexpr("\\)", out.xy[1]))-1),
                            xreg = base::paste0("~", out.xy[2]))

  comp[["mediator"]] <- base::list(y = NULL, xreg = NULL)

  for (ii in 1:base::length(mediator.formulas)) {

    # Deparse each formula in medformula to derive input the mediator regression models:
    med.xy = base::deparse(mediator.formulas[[ii]]) %>% base::paste(collapse = "") %>%
      (function(x) {base::gsub("[[:space:]]", "", x)}) %>% base::strsplit("~") %>% base::unlist()

    comp[["mediator"]]$y <- base::c(comp[["mediator"]]$y, med.xy[1])
    comp[["mediator"]]$xreg <- base::c(comp[["mediator"]]$xreg, base::paste0("~", med.xy[2]))

  }

  # Define class (and levels) of all variables:
  all.vars <- base::c(comp$mediator$y, base::strsplit(base::paste(base::gsub("~","", comp$mediator$xreg), collapse="+"), "+", fixed=TRUE) %>%
                  base::unlist() %>%
                  base::unique() %>%
                  base::setdiff(base::c("+","~")))

  comp[["variables"]] <- base::vector("list", base::length(all.vars))
  names(comp[["variables"]]) <- all.vars

  for (ii in 1:base::length(all.vars)) {
    var.class <- base::class(data[[all.vars[ii]]])
    comp[["variables"]][[all.vars[ii]]]$class <- var.class
    if (var.class=="factor") {
      comp[["variables"]][[all.vars[ii]]]$levels <- base::levels(data[[all.vars[ii]]])
    }
  }

  return(comp)

}

#' Calculate bootstrap confidence bands for the direct or indirect effect of interest
#'
#' @param object direct or indirect effect of interest (as obtained within a call to the function \code{effect}).
#' @param alpha the confidence level
#'
#' @return data.frame containing the unique event times, estimated effect, and lower and upper confidence bands
#' @export
#'
#' @keywords internal
add.ci <- function(object, alpha) {

  `%>%` <- dplyr::`%>%`

  # effect names (if factor then n.levels - 1 dummy variables, otherwise same as variable):
  effect.names <- base::names(object$coefs)[-1]

  # Add confidence bands:
  boot.coefs <- object$boot.coefs %>%
    dplyr::group_by(.data$times)

  # The cumulative sum may be NA especially for later time points due to rank deficient design matrices
  object$lower <- boot.coefs %>%
    dplyr::summarise_at(dplyr::vars(dplyr::one_of(effect.names)), function(x) stats::quantile(x, alpha/2, na.rm=TRUE))

  object$upper <- boot.coefs %>%
    dplyr::summarise_at(dplyr::vars(dplyr::one_of(effect.names)), function(x) stats::quantile(x, 1-alpha/2, na.rm=TRUE))

  return(object)

}

#' Wrapper function for survival::Surv
#'
#' @param ... parameters passed to survival::Surv
#'
#' @return object of class survival::Surv
#' @export
#'
#' @keywords internal
Surv <- function(...) {
  return(survival::Surv(...))
}
