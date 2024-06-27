#' Plot effects from dynamic path analysis along with bootstrap confidence bands
#'
#' @description plotting method for class "effect"
#'
#' @param object object of class "effect", or list of objects of class "effect"
#' @param relative should the effect be plotted on a relative survival scale (i.e. `y=exp(-effect)`)?. Defaults to FALSE.
#' @param titles If NULL, function will automatically generate. Otherwise character vector of length equal to number of
#' elements in object list
#' @param yintercept y-intercept of horizontal line. Defaults to 0.
#' @param linetype type of horinzontal line to plot. Defaults to "dashed".
#' @param x_label Label for x-axis. Defaults to "Time"
#' @param y_label Label for y-axis. Default will be "Cumulative Effect" when object scale is "cumulative", "Effect" otherwise.
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' library(dpasurv)
#'
#' data(simdata)
#'
#' s <- dpa(survival::Surv(start,stop,event)~M+x, list(M~x), id="subject", data=simdata, boot.n=100)
#'
#' direct <- effect(x ~ outcome, s)
#' indirect <- effect(x ~ M ~ outcome, s)
#' total <- sum(direct, indirect)
#'
#' ggplot.effect(direct)
#' ggplot.effect(list(direct, indirect, total))
#'
ggplot.effect <- function(object,
                          relative=FALSE,
                          titles =NULL,
                          yintercept = 0,
                          linetype = "dashed",
                          x_label = "Time",
                          y_label = NULL) {

  # Set global variables (called in the pipes below)
  times <- y <- ymin <- ymax <- group <- effect_type <- NULL

  `%>%` <- dplyr::`%>%`
  all_plot_dat <- dplyr::tibble()

  # Make sure we can handle both effect objects and list of effect objects
  if(class(object) == "effect") {
    object_list <- list(object)
  } else{
    object_list <- object
  }

  # We generate a seperate tibble for each each effect-object in the object list,
  # and then bind them together
  for(iteration in seq_along(object_list)){
    object <- object_list[[iteration]]
    # effect names (if factor then n.levels - 1 dummy variables, otherwise same as variable):
    effect.names <- base::names(object$coefs)[-1]


    # Plot all of the effects:
    for (ii in 1:base::length(effect.names)) {
      # data to be plotted
      plot_dat <- object$coefs %>% dplyr::select(dplyr::one_of(c("times", effect.names[ii]))) %>%
        dplyr::inner_join(object$lower %>%
                            dplyr::select(dplyr::one_of(c("times", effect.names[ii]))) %>%
                            dplyr::rename(lower=dplyr::contains(effect.names[ii])), by="times") %>%
        dplyr::inner_join(object$upper %>%
                            dplyr::select(dplyr::one_of(c("times", effect.names[ii]))) %>%
                            dplyr::rename(upper=dplyr::contains(effect.names[ii])), by="times")

      plot_dat$group <- effect.names[ii]

      plot_dat$effect_type <- gsub("\\+", "+\n", object$label)

      if(!is.null(titles)) plot_dat$effect_type <- titles[iteration]


      # Plot cumulative effects on "outcome", otherwise regular regression effects
      if(object$scale=="cumulative") {
        plot_dat$ylab <- ifelse(!is.null(y_label), y_label, "Cumulative Effect")
      } else {
        plot_dat$ylab <- ifelse(!is.null(y_label), y_label, "Effect")
      }


      # Define effect_type y argument:
      if (relative) {
        plot_dat$y <- exp(-plot_dat[[effect.names[ii]]])
      } else {
        plot_dat$y <- plot_dat[[effect.names[ii]]]
      }

      # Add confidence band to data:

      if (relative) {
        plot_dat$ymin <- exp(-plot_dat$lower)
      } else {
        plot_dat$ymin <- plot_dat$lower
      }


      if (relative) {
        plot_dat$ymax <- exp(-plot_dat$upper)
      } else {
        plot_dat$ymax <- plot_dat$upper
      }

      all_plot_dat <- dplyr::bind_rows(all_plot_dat,plot_dat)
      all_plot_dat$effect_type <- factor(all_plot_dat$effect_type, levels = unique(all_plot_dat$effect_type))

    }
  }

  # Now that all data has been created, do the actual plotting

  plot_object <- ggplot2::ggplot(data = all_plot_dat,
                                 ggplot2::aes(x = times, y = y, ymin = ymin, ymax = ymax)) +
    ggplot2::geom_ribbon(fill = "azure3") +
    ggplot2::geom_line() +
    ggplot2::ylab(unique(all_plot_dat$ylab)) +
    ggplot2::xlab(x_label) +
    ggplot2::geom_hline(yintercept = yintercept, color = "red", linetype = linetype) +
    ggplot2::theme_bw()

  if(dplyr::n_distinct(all_plot_dat$group)>1){
    plot_object <- plot_object +
      ggplot2::facet_grid(rows = ggplot2::vars(group),
                          cols = ggplot2::vars(effect_type))
  }else{
    plot_object <- plot_object +
      ggplot2::facet_grid(cols = ggplot2::vars(effect_type))
  }

  plot_object



}

#' Plot effects from dynamic path analysis along with bootstrap confidence bands
#'
#' @description plotting method for class "effect"
#'
#' @param x object of class "effect"
#' @param relative should the effect be plotted on a relative survival scale (i.e. `y=exp(-effect)`)?. Defaults to FALSE.
#' @param ... other graphical parameters passed to the graphics::plot function.
#'
#' @export
#'
#' @examples
#' library(dpasurv)
#'
#' data(simdata)
#'
#' s <- dpa(survival::Surv(start,stop,event)~M+x, list(M~x), id="subject", data=simdata, boot.n=100)
#'
#' direct <- effect(x ~ outcome, s)
#' indirect <- effect(x ~ M ~ outcome, s)
#' total <- sum(direct, indirect)
#'
#' par(mfrow=c(1,3))
#' plot(direct); abline(h=0, lty=2, col=2)
#' plot(indirect); abline(h=0, lty=2, col=2)
#' plot(total); abline(h=0, lty=2, col=2)
#'
plot.effect <- function(x, relative=FALSE, ...) {

  `%>%` <- dplyr::`%>%`

  # effect names (if factor then n.levels - 1 dummy variables, otherwise same as variable):
  effect.names <- base::names(x$coefs)[-1]

  # Plot all of the effects:
  for (ii in 1:base::length(effect.names)) {

    # Request user input: "Hit <Return> for next plot"
    if (ii > 1) grDevices::devAskNewPage(ask = TRUE)

    args <- base::list(...)

    if (base::length(effect.names) > 1) {
      args$main <- ifelse("main" %in% base::names(args), args$main, paste0(gsub("\\+", "+\n", x$label), "\n$", effect.names[ii]))
    } else {
      args$main <- ifelse("main" %in% base::names(args), args$main, gsub("\\+", "+\n", x$label))
    }

    args$xlab <- ifelse("xlab" %in% base::names(args), args$xlab, "Time")

    args$type = ifelse("type" %in% base::names(args), args$type, "s")

    # Plot cumulative effects on "outcome", otherwise regular regression effects
    if(x$scale=="cumulative") {
      args$ylab <- ifelse("ylab" %in% base::names(args), args$ylab, "Cumulative Effect")
    } else {
      args$ylab <- ifelse("ylab" %in% base::names(args), args$ylab, "Effect")
    }

    # data to be plotted
    dat.plot <- x$coefs %>% dplyr::select(dplyr::one_of(c("times", effect.names[ii]))) %>%
      dplyr::inner_join(x$lower %>%
                          dplyr::select(dplyr::one_of(c("times"))) %>%
                          dplyr::bind_cols(x$lower %>%
                                             dplyr::select(dplyr::one_of(c(effect.names[ii]))) %>%
                                             dplyr::rename(lower=dplyr::contains(effect.names[ii]))), by="times") %>%
      dplyr::inner_join(x$upper %>%
                          dplyr::select(dplyr::one_of(c("times"))) %>%
                          dplyr::bind_cols(x$upper %>%
                                             dplyr::select(dplyr::one_of(c(effect.names[ii]))) %>%
                                             dplyr::rename(upper=dplyr::contains(effect.names[ii]))), by="times")

    if (!("ylim" %in% base::names(args))) {
      args$ylim <- base::range(dat.plot[,-1])
    }

    args$x <- dat.plot$times

    # Define main y argument:
    if (relative) {
      args$y <- exp(-dat.plot[[effect.names[ii]]])
    } else {
      args$y <- dat.plot[[effect.names[ii]]]
    }

    # plot estimated effect:
    base::do.call(plot, args)

    # Add confidence band to plot:
    args$col = "grey"

    if (relative) {
      args$y <- exp(-dat.plot$lower)
    } else {
      args$y <- dat.plot$lower
    }

    base::do.call(graphics::lines, args)

    if (relative) {
      args$y <- exp(-dat.plot$upper)
    } else {
      args$y <- dat.plot$upper
    }

    base::do.call(graphics::lines, args)

  }

  # Deactivate Hit Return for next plot:
  if (ii > 1)
    grDevices::devAskNewPage(ask = FALSE)

}
