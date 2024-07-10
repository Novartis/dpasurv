#' Effect estimation in dynamic path analysis
#'
#' @description effect estimation method for class "dpa"
#'
#' @param formula the formula for the direct or indirect effect to be estimated. Should be of the form: covariate ~ outcome for direct effect of covariate on outcome, while it
#' should be of the form: covariate ~ mediator ~ outcome for indirect effect of covariate on outcome mediated through mediator.
#' Note that the word "outcome" is reserved for the survival outcome process, but the word "covariate" and "mediator" should match a corresponding variable name in the data input.
#' Alternatively the form can be: covariate ~ mediator for direct effects of covariate on mediator, or: covariate ~ mediator1 ~ mediator2
#' for indirect effects of covariate on mediator2 mediated through mediator1.
#' @param object object of class "dpa" (as obtained by calling the function \code{dpa}) from which the effect is to be estimated.
#' @param alpha The confidence level of the bootstrap intervals
#'
#' @return object of class "effect" with following fields:
#' \describe{
#' \item{coefs}{data.frame containing the unique event times along with the calculated effect coefficients. For effects corresponding to a continuous variable this results in a
#' single effect column. For factors with n.levels categories the data.frame contains n.levels-1 effect columns each representing the effect coefficients of a particular factor level (as compared to reference level).}
#' \item{lower}{data.frame of same dimension as coefs containing the lower confidence bands of the effects stored in coefs}
#' \item{upper}{data.frame of same dimension as coefs containing the upper confidence bands of the effects stored in coefs}
#' \item{boot.coefs}{data.frame with three columns: one column of bootstrap sample ID, a second column of unique event times (per bootstrap sample), and
#' a third column of the estimated effect coefficients (per bootstrap sample). The storing of the effects per bootstrap sample
#' facilitates calculation of bootstrap confidence intervals for sums of indirect and direct effects.}
#' \item{label}{effect label with path specification: "direct" for direct effect and "indirect" for indirect effect mediated through a path of mediator(s)}
#' \item{scale}{scale of effect coefficients in coefs, lower, upper: "cumulative" (for effects on outcome) or "identity" (for effects on mediators)}
#' \item{alpha}{confidence level of the bootstrap intervals}
#' }
#' @export
#'
#' @examples
#' library(dpasurv)
#'
#' data(simdata)
#'
#' set.seed(1)
#'
#' s <- dpa(Surv(start,stop,event)~M+x, list(M~x), id="subject", data=simdata, boot.n=50)
#'
#' direct <- effect(x ~ outcome, s)
#' indirect <- effect(x ~ M ~ outcome, s)
#' total <- sum(direct, indirect)
#'
effect <- function(formula, object, alpha=0.05) {

  # set up an empty output object (of class "effect"):
  output <- base::list(coefs = NULL, lower=NULL, upper=NULL, boot.coefs = NULL, label=NULL, formula=formula, scale=NULL, alpha=alpha)

  base::class(output) <- "effect"

  `%>%` <- dplyr::`%>%`

  # deparse formula to define path:
  path.variables <- find.variables(formula)

  if (base::length(path.variables) == 2) {
    output$label <- base::paste0("direct(",base::format(formula),")")
  } else {
    output$label <- base::paste0("indirect(",base::format(formula),")")
  }

  # Check the class of the covariate whose effect to calculate (and then find column names corresponding to associated variables.
  # If covariate is continuous then the variable name is trivial, while for factors there are n.levels-1 dummy variable names
  # corresponding to comparison with reference level[1]
  if (object$meta$variables[[path.variables[1]]]$class == "factor") {
    cols <- base::paste0(path.variables[1], object$meta$variables[[path.variables[1]]]$levels[-1])
  } else {
    cols <- path.variables[1]
  }

  # Initialize lists by calculating the direct effect of first variable on second variable (on path):
  output$coefs <- object$coefs[[path.variables[2]]] %>%
    dplyr::select(dplyr::one_of(base::c("times", cols)))

  output$boot.coefs <- object$boot.coefs[[path.variables[2]]] %>%
    dplyr::select(dplyr::one_of(base::c("boot.id", "times", cols)))

  # Calculate indirect effect if length(path.variables) > 2
  if (base::length(path.variables) > 2) {

    for (ii in 3:base::length(path.variables)) {
      output$coefs <- output$coefs %>%
        dplyr::mutate_at(dplyr::vars(-.data$times), function(x) x*(object$coefs[[path.variables[ii]]] %>%
                                                               dplyr::select(dplyr::one_of(path.variables[ii-1])) %>%
                                                               dplyr::pull()))

      output$boot.coefs <- output$boot.coefs %>%
        dplyr::mutate_at(dplyr::vars(-.data$boot.id, -.data$times), function(x) x*(object$boot.coefs[[path.variables[ii]]] %>%
                                                                         dplyr::select(dplyr::one_of(path.variables[ii-1])) %>%
                                                                         dplyr::pull()))
    }

  }

  # transform effects onto cumulative scale (if needed; otherwise identity transformation):
  # We calculate cumulative effects on "outcome", otherwise regular regression effects

  if(path.variables[base::length(path.variables)]=="outcome") {
    fun <- base::cumsum
    output$scale <- "cumulative"
  } else {
    fun <- base::identity
    output$scale <- "identity"
  }

  output$coefs <- output$coefs %>%
    dplyr::mutate_at(dplyr::vars(dplyr::one_of(cols)), fun)

  output$boot.coefs <- output$boot.coefs %>%
    dplyr::group_by(.data$boot.id) %>%
    dplyr::mutate_at(dplyr::vars(dplyr::one_of(cols)), fun) %>%
    dplyr::ungroup()

  # Add confidence bands:
  output <- add.ci(output)

  return(output)

}




#' Sum of direct and indirect effects from dynamic path analysis
#'
#' @description a sum method for class "effect"
#'
#' @param effect1 an object of class "effect" obtained from a call to the function effect()
#' @param effect2 a second object of class "effect" obtained from a call to the function effect()
#' @param ... additional objects of class "effect" (any number of effects allowed)
#'
#' @return an object of class "effect" containing the sum of effect1, effect2, ...
#' @export
#'
#' @examples
#' library(dpasurv)
#'
#' data(simdata)
#'
#' set.seed(1)
#'
#' s <- dpa(Surv(start,stop,event)~M+x, list(M~x), id="subject", data=simdata, boot.n=50)
#'
#' direct <- effect(x ~ outcome, s)
#' indirect <- effect(x ~ M ~ outcome, s)
#' total <- sum(direct, indirect)
#'
sum.effect <- function(effect1, effect2, ...) {

  # set up an empty output object (of class "effect"):
  output <- base::list(coefs = NULL, lower=NULL, upper=NULL, boot.coefs = NULL, label=NULL, formula=NULL, scale=NULL, alpha=NULL)

  base::class(output) <- "effect"

  `%>%` <- dplyr::`%>%`

  # Gather effects into list (and remove na.rm field, which occurs by default of sum function)
  effect.list <- base::list(effect1, effect2, ...); effect.list[["na.rm"]] <- NULL

  classes <- base::unlist(base::lapply(effect.list, class))
  alphas <- base::unlist(base::lapply(effect.list, function(x) x$alpha))
  scales <- base::unlist(base::lapply(effect.list, function(x) x$scale))
  formulas <- base::unlist(base::lapply(effect.list, function(x) x$formula))
  first.vars <-base::unlist(base::lapply(formulas, function(x) find.variables(x)[1]))
  last.vars <- base::unlist(base::lapply(formulas, function(x) find.variables(x)[base::length(find.variables(x))]))

  if (base::sum(classes!="effect") > 0)
    stop("input should only consist of objects of class 'effect' as obtained from calling the function 'dpasurv::effect'")

  if (base::sum(base::diff(alphas)!=0) > 0)
    stop("the confidence levels (alpha) are not the same across the effects that are being summed. Please recalculate effects with the same alpha levels.")

  if (!base::all(scales==scales[1]))
    stop("the effects that are being summed must be on the same scale (either 'cumulative' for effects on outcome or 'identity' for effects on mediators).")

  if (!base::all(first.vars==first.vars[1]))
    stop("the effects are not compatible: the first variable in the formulas effect1$formula, effect2$formula etc. do not match")

  if (!base::all(last.vars==last.vars[1]))
    stop("the effects are not compatible: the last variable in the formulas effect1$formula, effect2$formula etc. do not match")

  output$alpha <- effect.list[[1]]$alpha
  output$scale <- effect.list[[1]]$scale
  output$formula <- formulas
  output$label <- base::paste(base::unlist(base::lapply(effect.list, function(x) x$label)),collapse=" + ")

  output$coefs <- effect.list[[1]]$coefs
  output$boot.coefs <- effect.list[[1]]$boot.coefs

  for (ii in 2:base::length(effect.list)) {
    output$coefs[,-1] <- output$coefs[,-1] + effect.list[[ii]]$coefs[,-1]
    output$boot.coefs[,-(1:2)] <- output$boot.coefs[,-(1:2)] + effect.list[[ii]]$boot.coefs[,-(1:2)]
  }

  output$coefs <- output$coefs %>% dplyr::as_tibble()
  output$boot.coefs <- output$boot.coefs %>% dplyr::as_tibble()

  output <- add.ci(output)

  return(output)

}
