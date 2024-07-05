#' Dynamic Path Analysis
#'
#' @param out.formula Survival formula for Aalen's additive hazards model.
#' @param mediator.formulas Mediator regression formula (in case of a single mediator), or a list of regression formulas (in case of multiple mediators).
#' The formulas must be ordered according to Directed Acyclic Graph Structure (see Details).
#' @param id character string indicating which column of 'data' corresponds to the subject ID. Bootstrapping will be performed on this id.
#' @param data Data set in counting process format. In particular the data should contain a "start", "stop" and "event" column along with
#' any mediators and baseline covariates.
#' @param boot.n Number of bootstrap samples.
#' @param method The underlying implementation of Aalen's additive regression model. Defaults to "aareg", which relies on the survival::aareg() implementation,
#' while method = "timereg" uses the timereg::aalen() implementation.
#' @param progress_bar Boolean. If TRUE, show progress bar. Defaults to FALSE.
#' @param ... other parameters passed to the Aalen's additive hazards model implementation. If method = "timereg", then ... will be passed to
#' timereg::aalen(), while if method = "aareg", then ... will be passed to survival::aareg(). If ... contains parameters that don't belong to the formalArgs of
#' the corresponding implementation then those parameters will be ignored.
#'
#' @return Object of class `dpa` with following fields:
#' \describe{
#'   \item{coefs}{list of estimated coefficients from each of the regressions listed in out.formula and mediator.formulas.}
#'   \item{boot.coefs}{list of bootstrap estimates corresponding to coefs. This stores all the bootstrap estimates to facilitate
#' calculation of direct, indirect and total effects along with bootstrap confidence intervals.}
#'   \item{meta}{a list keeping track of responses and covariates of each of the out.formula and mediator.formulas. Also keeps
#' track of all variable types and level names in case of factors.}
#' }
#' @export
#'
#' @details \code{dpa} performs Dynamic Path Analysis of a Directed Acyclic Graph (DAG). The out.formula
#' can have as covariates all mediators listed in mediator.formulas. The mediator.formulas must obey the
#' following DAG structure rule: The response of the k-th formula cannot appear as covariate in any of the formulas
#' k+1, ..., length(mediator.formulas).
#
#' @examples
#' library(dpasurv)
#'
#' data(simdata)
#'
#' set.seed(1)
#'
#' # Perform dynamic path analysis:
#' s <- dpa(survival::Surv(start,stop,event)~M+x, list(M~x), id="subject", data=simdata, boot.n=50)
#'
#' # Calculate cumulative direct, indirect, and total effects:
#' direct <- effect(x ~ outcome, s)
#' indirect <- effect(x ~ M ~ outcome, s)
#' total <- sum(direct, indirect)
#'
#' # Plot the effects using basic graphics:
#' par(mfrow=c(1,3))
#' plot(direct); abline(h=0, lty=2, col=2)
#' plot(indirect); abline(h=0, lty=2, col=2)
#' plot(total); abline(h=0, lty=2, col=2)
#'
#' # Plot the effects using ggplot2 graphics:
#' ggplot.effect(list(direct, indirect, total))
#'
dpa <- function(out.formula, mediator.formulas, id, data, boot.n=100, method = "aareg", progress_bar = FALSE, ...) {

  `%>%` <- dplyr::`%>%`

  # Make sure all categorical variables in data are set as factors:
  data <- dplyr::mutate_if(data, base::is.character, base::as.factor)

  # if mediator.formulas is a single formula (instead of a list of formulas), transform to list:
  if (base::class(mediator.formulas)=="formula")
    mediator.formulas <- base::list(mediator.formulas)

  # Parse all of the formulas into distinct (left & right hand side of ~) formula components:
  meta <- get.meta(out.formula, mediator.formulas, data)

  # Make sure that the formulas define a valid DAG:
  check.dag(meta, data)

  # How many mediator formulas?
  num.mediators <- base::length(mediator.formulas)

  # Arrange data by Subject and stop times:
  data <- data %>% dplyr::arrange(!!base::as.symbol(id), !!base::as.symbol(meta$outcome$stopt))

  # Keep track of observed times (of event). This gets passed to the function resolve.ties to
  # guarantee that (in case of ties) we don't generate death times that exceed other subject's death times
  obstimes <- unique(base::sort(data[[meta$outcome$stopt]][data[[meta$outcome$event]]==1]))

  # Resolve ties (if any) by adding random noise to duplicates:
  # data[[meta$outcome$stopt]] <- resolve.ties(data[[meta$outcome$stopt]], data[[meta$outcome$event]], obstimes)

  #########################
  # Obtain estimates
  #########################

  # List of coefficients from the different models (and corresponding bootstrap coefficients):
  coefs <- base::list()

  arguments <- list(out.formula = out.formula, data = data, id = id, method = method)

  dot.args <- list(...)

  # Retrieve and summarise coefs under "timereg" implementation
  if (method == "timereg") {

    # Aalen's additive hazard model:
    areg.obj <- do.call(Areg, c(arguments, dot.args[base::intersect(methods::formalArgs(timereg::aalen), names(dot.args))]))

    # Undo the random tie-breaking from timereg::aalen() and summarise the coefficients
    # at unique observed times. We sum up the coefficients across ties within unique times
    # since we are only interested in cumulative effects anyway. The mediator regression below
    # will only be applied at the unique event times and would be constant across tied
    # survival times. So this strategy also works for cumulative indirect effects.
    coefs[["outcome"]] <- areg.obj$coefs %>%
      dplyr::mutate(time.bins = base::cut(.data$times, breaks = c(obstimes, Inf), include.lowest=FALSE, labels = obstimes, right=FALSE)) %>%
      dplyr::group_by(.data$time.bins) %>%
      dplyr::mutate(times = .data$times[1L]) %>%
      dplyr::ungroup() %>%
      dplyr::select(-.data$time.bins) %>%
      dplyr::group_by(.data$times) %>%
      dplyr::summarise_all(sum)

  } else { # Retrieve and summarise coefs under "aareg" implementation

    # Aalen's additive hazard model:
    areg.obj <- do.call(Areg, c(arguments, dot.args[base::intersect(methods::formalArgs(survival::aareg), names(dot.args))]))

    coefs[["outcome"]] <- areg.obj$coefs %>%
      dplyr::group_by(.data$times) %>%
      dplyr::summarise_all(sum)

  }

  # Mediator models:
  for (ii in 1:num.mediators)
    coefs[[meta$mediator$y[ii]]] <- Mreg(regformula = stats::as.formula(meta$mediator$xreg[ii]),
                                         obstimes = coefs[["outcome"]]$times, startt = meta$outcome$startt,
                                         stopt = meta$outcome$stopt, event = meta$outcome$event,
                                         mediator = meta$mediator$y[ii], dataset = data, w=1)

  #############################
  # Obtain bootstrap estimates
  #############################

  # Create empty list of bootstrapped coefficients:
  boot.coefs <- base::vector("list", base::length(coefs))
  base::names(boot.coefs) <- base::names(coefs)

  if(progress_bar) {
    pb <- utils::txtProgressBar(min = 1,
                                max = boot.n,
                                initial = 1,
                                style = 3)
  }

  # Bootstrap:

  # Bootstrapping is performed by id:
  boot.by=id

  arguments[["id"]] <- "bootstrapID"

  if (method == "timereg") {
    # Set/Overwrite n.sim = 0 for all the bootstraps:
    dot.args[["n.sim"]] <- 0
  }

  for (b in 1:boot.n) {

    if(progress_bar) {
      utils::setTxtProgressBar(pb, b)
    }

    if (!base::is.null(boot.by)) {
      nested_by_ID <- data %>%
        dplyr::group_by(!!base::as.symbol(boot.by)) %>%
        tidyr::nest()

      # Override with bootstrap sample:
      nested_by_ID <- nested_by_ID[base::sample(dplyr::n_distinct(data[,boot.by]), replace = T),]

      # Create a unique bootstrapID for each sample that will be passed to timereg::aalen()
      nested_by_ID$bootstrapID <- 1:base::nrow(nested_by_ID)

      boot.data <- nested_by_ID %>%
        tidyr::unnest(data) %>%
        dplyr::ungroup()

    } else {
      boot.sample <- base::sample(1:base::nrow(data), replace=TRUE)
      boot.data <- data[boot.sample,]
      boot.data$bootstrapID <- 1:base::nrow(boot.data)
    }

    # Resolve ties by adding random noise to duplicates (by supplying all obstimes we guarantee against generating
    # death times of a tied subject that exceeds another subject's death time):
    # boot.data[[meta$outcome$stopt]] <- resolve.ties(boot.data[[meta$outcome$stopt]], boot.data[[meta$outcome$event]], obstimes)

    arguments[["data"]] <- boot.data

    # Retrieve and summarise coefs under "timereg" implementation
    if (method == "timereg") {

      # Aalen's additive hazard model:
      areg.obj.boot <- base::do.call(Areg, c(arguments, dot.args[base::intersect(methods::formalArgs(timereg::aalen), names(dot.args))]))

      # Undo the random tie-breaking from timereg::aalen() and summarise the coefficients
      # at unique observed times. We sum up the coefficients across ties within unique times
      # since we are only interested in cumulative effects anyway. The mediator regression below
      # will only be applied at the unique event times and would be constant across tied
      # survival times. So this strategy also works for cumulative indirect effects.
      boot.coefs[["outcome"]][[b]] <- areg.obj.boot$coefs %>%
        dplyr::mutate(time.bins = base::cut(.data$times, breaks = c(obstimes, Inf), include.lowest=FALSE, labels = obstimes, right=FALSE)) %>%
        dplyr::group_by(.data$time.bins) %>%
        dplyr::mutate(times = .data$times[1L]) %>%
        dplyr::ungroup() %>%
        dplyr::select(-.data$time.bins) %>%
        dplyr::group_by(.data$times) %>%
        dplyr::summarise_all(sum) %>%
        dplyr::mutate(boot.id = b)

    } else { # Retrieve and summarise coefs under "aareg" implementation

      # Aalen's additive hazard model:
      areg.obj.boot <- base::do.call(Areg, c(arguments, dot.args[base::intersect(methods::formalArgs(survival::aareg), names(dot.args))]))

      boot.coefs[["outcome"]][[b]] <- areg.obj.boot$coefs %>%
        dplyr::group_by(.data$times) %>%
        dplyr::summarise_all(sum) %>%
        dplyr::mutate(boot.id = b)

    }

    # Mediator models:
    for (ii in 1:num.mediators)
      boot.coefs[[meta$mediator$y[ii]]][[b]] <- Mreg(regformula = stats::as.formula(meta$mediator$xreg[ii]),
                                                      obstimes = boot.coefs[["outcome"]][[b]]$times, startt = meta$outcome$startt,
                                                      stopt = meta$outcome$stopt, event = meta$outcome$event,
                                                      mediator = meta$mediator$y[ii], dataset = boot.data, w=1) %>%
      dplyr::mutate(boot.id = b)

  }

  if(progress_bar) {
    close(pb)
  }

  # Bind the list of bootstrap coefficients together into single data.frames
  boot.coefs[["outcome"]] <- dplyr::bind_rows(boot.coefs[["outcome"]])
  for (ii in 1:num.mediators)
    boot.coefs[[meta$mediator$y[ii]]] <- dplyr::bind_rows(boot.coefs[[meta$mediator$y[ii]]])

  dpa.obj <- base::list(coefs = coefs, boot.coefs = boot.coefs, meta=meta, aalen=areg.obj$aalen)

  base::class(dpa.obj) <- "dpa"

  return(dpa.obj)

}
