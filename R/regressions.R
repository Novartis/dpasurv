#' Internal mediator regression function storing a linear regression estimate for each event time
#'
#' Do not call this function on its own
#'
#' @param regformula independent variable input as for standard lm fuction in R
#' @param obstimes observed event times at which Aalen's additive model is estimated via the function Areg
#' @param startt start time of the observation intervals
#' @param stopt end time of the observation intervals (actual event times)
#' @param event event indicator
#' @param mediator mediator of interest
#' @param dataset dataset (in counting process format)
#' @param w (optional) weights (not actually used in the default implementation of dynamic path analysis, set to 1)
#'
#' @return data.frame with observation times and estimated coefficients for independent variables in "regformula"
#' @export
#'
#' @keywords internal
Mreg <- function(regformula, obstimes, startt, stopt, event, mediator, dataset, w=1) {

  `%>%` <- dplyr::`%>%`

  # Design matrix:
  X <- stats::model.matrix(regformula, dataset)

  # Mediator:
  M <- dataset[[mediator]]

  # Starting times and stopping times:
  starts <- dataset[[startt]]
  stops <- dataset[[stopt]]

  #Ridge parameter epsilon
  epsilon = 1e-08

  subsets <- outer(starts, obstimes, function(x, y) x < y) &
    outer(stops, obstimes, function(x, y) x >= y)

  #Matrix to store b for each ii
  b.matrix = matrix(NA, ncol=ncol(X), nrow=length(obstimes))
  colnames(b.matrix) <- colnames(X)

  #Loop doing a linear regression in each time point:
  for(ii in 1:length(obstimes)) {

    # Design matrix on subset ii:
    Y <- X[subsets[,ii],,drop=FALSE]

    # define mediator on subset ii:
    E <- M[subsets[,ii]]

    # Least squares solution with ridge penalty for numerical stability
    b.matrix[ii,] <- base::solve((base::t(Y) %*% (w * Y) +
                                    epsilon * base::diag(ncol(X))), base::colSums((w * E) * Y))

  }

  #Output:
  base::cbind(times = obstimes, b.matrix) %>% dplyr::as_tibble()

}


#' Internal Aalen's Additive Hazards Model storing a linear effect estimate for each event time
#'
#' Do not call this function on its own
#'
#' @param out.formula Survival formula for Aalen's additive hazards model.
#' @param id character string indicating which column of 'data' corresponds to the subject ID.
#' @param data Data set in counting process format. In particular the data should contain a "start", "stop" and "event" column along with
#' any mediators and baseline covariates.
#' @param method passed from dpa call, defaults to "timereg", otherwise "aareg"
#' @param ... other parameters passed to Aalen's additive hazards regression function "timereg::aalen()"
#'
#' @return data.frame with observation times and estimated coefficients for independent variables in "regformula"
#' @importFrom rlang .data
#' @export
#'
#' @keywords internal
Areg = function(out.formula, id, data, method, ...) {

  `%>%` <- dplyr::`%>%`

  # Fit Aalen's additive model with timereg package
  if (method == "timereg") {

    # Aalen's additive hazard model:
    aalen.obj <- timereg::aalen(out.formula, data = data, id = data[[id]], ...)

    # Gather non-cumulative effect estimates from the timereg::aalen() regression
    aalen.coefs <- aalen.obj$cum[-1,] %>%
      dplyr::as_tibble() %>%
      dplyr::rename(times = .data$time) %>%
      dplyr::mutate_at(dplyr::vars(-c("times")), function(x) c(x[1],diff(x)))

    if (is.matrix(aalen.obj$gamma)) {

      const.params <- base::outer(c(aalen.coefs$times[1L], base::diff(aalen.coefs$times)) , base::as.numeric(aalen.obj$gamma), FUN="*")

      base::colnames(const.params) <- base::row.names(aalen.obj$gamma)
      base::colnames(const.params) <- base::gsub("const\\(","", base::colnames(const.params))
      base::colnames(const.params) <- base::substr(base::colnames(const.params), 1, base::nchar(base::colnames(const.params))-1)

      aalen.coefs <- aalen.coefs %>%
        dplyr::bind_cols(const.params)

    }

    names(aalen.coefs) <- base::gsub("\\(Intercept\\)", "Intercept", names(aalen.coefs))
    names(aalen.coefs) <- base::gsub("\\(Intercept\\)", "Intercept", names(aalen.coefs))

  } else { # Fit Aalen's model using the default method = "aareg"

    # Aalen's additive hazard model:
    aalen.obj <- survival::aareg(out.formula, data = data, ...)

    rownames(aalen.obj$coefficient) <- NULL

    # Gather non-cumulative effect estimates from the survival::aareg() regression
    aalen.coefs <- cbind(data.frame(times=aalen.obj$times), aalen.obj$coefficient) %>%
      dplyr::as_tibble()

  }

  return(list(aalen = aalen.obj, coefs = aalen.coefs))

}
