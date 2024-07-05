##############################################################################
# Load packages
##############################################################################

library(dplyr, warn.conflicts = FALSE)
library(tidyr)

##############################################################################
# Data generating mechanism
##############################################################################
# Helper function that simulates trials

sim_trials <- function(id,
                       design = list(n_patients,
                                     measurement_weeks,
                                     time_grid),

                       parameters = list(
                         alpha_0,
                         alpha_1,
                         M_random_intercept_sd,
                         # Random intercept for M
                         M_noise_sd,
                         # Noise around measurements of M
                         beta_0,
                         beta_1,
                         M_eff_on_haz,
                         # constant effect that mediator has on hazard#
                         baseline_hazard,
                         censor_rate = .00005 # rate of the exp. distr. that simulates censoring
                       ),
                       parallel = FALSE)
{

  measurement_weeks_index <-
    which(!duplicated(findInterval(
      design$time_grid, design$measurement_weeks
    )))[-1]

  #############################################################################
  # Simulate
  #############################################################################
  # Randomise treatment/dose variable x
  x <- sample(x = c(0, 1, 2),
              size = design$n_patients,
              replace = TRUE)

  # Simulate individual patients' parameter values

  # Generate mediator matrix M
  M_true <-
    rep(1, length(x)) %*% t(parameters$alpha_0) + x %*% t(parameters$alpha_1)
  if (!is.null(parameters$M_random_intercept_sd)) {
    M_intercept <-
      rnorm(design$n_patients, 0, parameters$M_random_intercept_sd)
  } else{
    M_intercept <- 0
  }
  M_true <- M_true + M_intercept
  # Add noise around M_true
  # So far, assume independence across observations
  Sigma <-
    matrix(0, length(design$time_grid), length(design$time_grid))
  diag(Sigma) <- (parameters$M_noise_sd) ^ 2
  noise <-
    MASS::mvrnorm(n = design$n_patients, rep(0, length(design$time_grid)), Sigma = Sigma)
  M_noise <- (M_true + noise)
  # Not all the values we generate for M will actually be measured
  M_measured <- M_noise[, measurement_weeks_index] %>% data.frame() %>% as_tibble()
  # M_measured <- M_noise
  names(M_measured) <- paste0("M", design$measurement_weeks)
  # combine into a data set
  sim_dat <- bind_cols(
    subject = 1:design$n_patients,
    x = x,
    M0 = M_intercept + parameters$alpha_0[1],
    M = M_measured
  )

  M_eff_on_haz <-
    matrix(
      rep(parameters$M_eff_on_haz, design$n_patients),
      byrow = T,
      nrow = design$n_patients
    )
  M_eff_on_haz_mat <- (M_noise * M_eff_on_haz)
  M_eff_on_haz_mat <- as_tibble(data.frame(M_noise * parameters$M_eff_on_haz))

  # We now have everything to define the hazard matrix
  hazard_mat <-
    rep(1, length(x)) %*% t(parameters$beta_0) + sim_dat$x %*% t(parameters$beta_1) + M_eff_on_haz_mat
  bh <- 0
  if (min(hazard_mat) < 0) {
    bh <- abs(min(hazard_mat))
    hazard_mat <- hazard_mat + bh
  }
  names(hazard_mat) <- paste0("week_", design$time_grid)
  M_at_study_end <- paste0("M", tail(design$measurement_weeks, 1))
  sim_dat_long <-
    pivot_longer(sim_dat, cols = starts_with("M"), names_to = "M") %>%
    filter(M != M_at_study_end) %>% # last mediator measurement is irrelevant
    mutate(week = rep(c(0, design$measurement_weeks[-length(design$measurement_weeks)]), design$n_patients)) %>%
    select(-M) %>%
    rename(M = value)

  # For each patient, we now have one hazard value for each week
  # This approximates the probability of having an event in that week
  event_table <- t(apply(hazard_mat*diff(design$time_grid)[1], 1, rbinom, n = length(design$time_grid), size = 1))
  # The last column is irrelevant for us, because this depicts the probability of having an event in the week after the study finishes
  # The first column is also ignored, because by definition everyone is alive at baseline

  # Which patients have an event, and when?
  identify_event <- function(x){
    if(sum(x)==0) return(tibble(event = 0, time = tail(design$time_grid,1))) # no event
    return(tibble(event = 1, time = design$time_grid[min(which(x==1))]))
  }

  event_times_per_patient <-
    apply(event_table, 1, identify_event) %>%
    bind_rows() %>%
    mutate(subject = 1:design$n_patients)

  ############################################################################
  # Censoring
  ############################################################################

  # Each patient's observations have a chance of being censored before the
  # study finishes
  censoring_times <- round(rexp(design$n_patients, parameters$censor_rate), 2)
  # Every patient gets censored when the study finishes
  censoring_times <-
    pmin(censoring_times, design$measurement_weeks[length(design$measurement_weeks)])

  event_times_per_patient <- event_times_per_patient %>%
    mutate(
      event = ifelse(censoring_times <= time, 0, event),
      time = ifelse(censoring_times <= time, censoring_times, time)
    )


  temp <-
    tmerge(
      event_times_per_patient,
      event_times_per_patient,
      id = subject,
      event = event(time, event)
    )
  sim_dat_long <-
    tmerge(
      temp,
      sim_dat_long,
      id = subject,
      M = tdc(week, M),
      x = tdc(week, x)
    ) %>% data.frame() %>%
    as_tibble()

  list(sim_dat = sim_dat_long %>% select(subject, x, tstart, tstop, M, event) %>% rename(start=tstart,stop=tstop),
       baseline_hazard = bh)
}

# Helper functions for the Emax model:
f_constant <- function(t, c) {
  return(rep(c, length(t)))
}

f_emax <- function(t, e_max, et_50, hill=1, e_0 = 0) {
  return(t^hill * e_max/(t^hill + et_50^hill) + e_0)
}

#############################################################################
# Set up design parameters
#############################################################################
# number of patients
n_patients <- 200

# The weeks measurements are made after treatment
measurement_weeks <- c(4,8, 12, 26 , 52, 78, 104, 156, 260) # the last number depicts not a measurement, but the end of the study

time_grid <- 1:260

#############################################################################
# Set up parameter values
#############################################################################

#############################################################################
# Define mediator value under treatment and control
#############################################################################

# Set model parameters that models the impact of treatment on the hazard over time
e_0_M <- 6.8 # mean value of M for control patients (and at time = 0 for patients on treatment)
e_max_M <- 3.7/2 # max treatment effect, reached as time approaches infinity
et_50_M <- 30 # time at which half of the treatment effect is achieved

M_eff_on_haz <- rep(0.00045, length(time_grid))

alpha_0 <- f_constant(t=time_grid, c=e_0_M) # Mediator value for patients on control
alpha_1 <- f_emax(t=time_grid, e_max=e_max_M, et_50=et_50_M, hill=1) # alpha_0 + alpha_1 yields mediator value for patients on experimental

M_random_intercept_sd <- 1.5

M_noise_sd <- sqrt(0.05)

#############################################################################
# Define Treatment effect on hazard (and baseline hazard)
#############################################################################

e_0_X <- 0 # Impact of treatment on hazard at time 0
e_max_X <- -.007/2 # treatment effect, reached as time approaches infinity
et_50_X <- 30 # time at which half of the treatment effect is achieved
h_X <- 2

baseline_hazard <- 0.0055

beta_0 <- f_constant(t=time_grid, c=baseline_hazard)
beta_1 <- f_emax(t=time_grid, e_max=e_max_X, et_50 = et_50_X, hill = h_X, e_0 = e_0_X)

#############################################################################
# Define baseline hazard
#############################################################################

# Save parameters values

design = list(n_patients = n_patients,
              measurement_weeks = measurement_weeks,
              time_grid = time_grid)

parameters = list(alpha_0 = alpha_0,
                  alpha_1 = alpha_1,
                  M_random_intercept_sd = M_random_intercept_sd, # Random intercept for M
                  M_noise_sd = M_noise_sd, # Noise around measurements of M
                  beta_0 = beta_0,
                  beta_1 = beta_1,
                  M_eff_on_haz = M_eff_on_haz, # constant effect that mediator has on hazard#
                  baseline_hazard = baseline_hazard,
                  censor_rate = .0006 # rate of the exp. distr. that simulates censoring
)

######################################################################################################
# Simulation and analysis
######################################################################################################

# Create dataset:
set.seed(254)
simdata <- sim_trials(1, design=design, parameters=parameters)$sim_dat %>%
  mutate(dose = factor(ifelse(x==0, "ctrl", ifelse(x==1, "low", "high")), levels=c("ctrl", "low", "high")),
         x = factor(ifelse(x==0, "ctrl", "trt")),
         stop = ifelse(event==1, stop - round(runif(n(),min=0, max=1), 2), stop),
         subject = factor(subject)) %>%
  relocate(dose, M, .after=x)

usethis::use_data(simdata, overwrite = TRUE)
