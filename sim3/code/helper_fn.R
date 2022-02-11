#############################################################################
#                               HELPER FUNCTIONS                            #
#############################################################################

#' generate data in a setting similar to what might be expected 
#' in a correlates study of a COVID-19 vaccine
#' 
#' @param n Sample size (40,000 for JnJ)
#' @param cens_rate The intercept related to the hazard of right censoring. 
#' @param covid_rate The intercept in the outcome regression. Make larger to 
#' add more events to mimic what might happen with further follow-up time
#' @study_stop The stop date of the study
#' @return A data.frame with self-explanatory column names

make_ows_data_survival <- function(
  n = 40000, cens_rate = -8, covid_rate = -6.5,
  study_stop = 67
){
  # baseline covariates
  age <- rbinom(n, 1, 0.4)
  race <- rbinom(n, 1, 0.25)
  risk <- rbinom(n, 1, 0.25)
  
  # vaccine distribution
  vax <- rbinom(n, 1, 0.5)
  
  # probability of getting 0 immune response among participants
  ab_prob <- plogis(2 - 3 * age)
  ab_binom <- rbinom(n, 1, ab_prob)
  # generate continuous nonzero immune marker truncated at 0
  mean_ab_vax_cov <- 3 - 0.5 * age
  sd_ab_vax_cov <- 1
  ab = numeric(n)
  ab[ab_binom == 1] <- rnorm(sum(ab_binom==1), mean = mean_ab_vax_cov, sd = sd_ab_vax_cov)
  ab[vax == 0] <- 0
  ab[vax == 1 & ab < 0] <- 0
  
  # hazard function of getting infected
  haz_covid_vax_ab_cov <- plogis(covid_rate - 4 * ab - 2.2 * vax + 0.5 * risk + 0.5 * age + 0.1 * race)
  # time to covid infection
  time_to_covid <- rgeom(n, haz_covid_vax_ab_cov) + 1

  # hazard function of right censoring
  haz_right_cens_vax_ab_cov <- plogis(cens_rate + 1.5*age + risk)
  # time to right censoring
  time_to_right_cens <- rgeom(n, haz_right_cens_vax_ab_cov) + 1
  
  # observed failure time
  observed_ftime <- pmin(time_to_covid, time_to_right_cens)
  # observed indicator
  observed_covid_ind <- as.numeric(observed_ftime == time_to_covid)
  # add censoring due to study stop
  observed_ftime[observed_ftime > study_stop] <- study_stop
  observed_covid_ind[observed_ftime > study_stop] <- 0

  # ows sampling design
  strata_matrix <- expand.grid(age = c(0,1), race = c(0,1), risk = c(0,1))
  vax_strata_size <- 113
  placebo_strata_size <- 15
  random_subcohort <- rep(0, n)
  for(i in seq_len(nrow(strata_matrix))){
    random_subcohort[age == strata_matrix$age[i] & 
                       race == strata_matrix$race[i] & 
                       risk == strata_matrix$risk[i] & 
                       vax == 1][1:vax_strata_size] <- 1
    random_subcohort[age == strata_matrix$age[i] & 
                       race == strata_matrix$race[i] & 
                       risk == strata_matrix$risk[i] & 
                       vax == 0][1:placebo_strata_size] <- 1				   
  }
  
  measure_ab <- random_subcohort
  
  # sample all cases
  measure_ab[observed_covid_ind == 1] <- 1
  
  # "sample" abs
  ab[measure_ab == 0] <- NA
  
  dat <- data.frame(age = age, race = race, risk = risk,
                    vax = vax, ab = ab, 
                    ftime = observed_ftime,
                    ftype = observed_covid_ind,
                    measure_ab = measure_ab, 
                    random_subcohort = random_subcohort)
  return(dat)
}


#' Compute the true value of four risk estimators based on the data
#' generating process encoded by \code{make_ows_data}.
#' 
#' @param n Sample size (should be very large to have accurate approximation)
#' @param covid_rate The intercept in the outcome regression. Make larger to 
#' add more events to mimic what might happen with further follow-up time
#' @param t0 the time for evaluation
#' @return the true value of four risk estimators
get_ows_truth1_survival <- function(
  n = 1e6, 
  covid_rate = -6.5, parm_ab = -4, parm_vax = -2.2,
  t0 = 66, study_stop = 67
){
  stopifnot(t0 <= study_stop)
  
  # baseline covariates
  age <- rbinom(n, 1, 0.4)
  race <- rbinom(n, 1, 0.25)
  risk <- rbinom(n, 1, 0.25)
  
  # vaccine distribution
  vax <- rbinom(n, 1, 0.5)
  
  # probability of getting 0 immune response among participants
  ab_prob <- plogis(2 - 3 * age)
  ab_binom <- rbinom(n, 1, ab_prob)
  # generate continuous nonzero immune marker truncated at 0
  mean_ab_vax_cov <- 3 - 0.5 * age
  sd_ab_vax_cov <- 1
  ab = numeric(n)
  ab[ab_binom == 1] <- rnorm(sum(ab_binom==1), mean = mean_ab_vax_cov, sd = sd_ab_vax_cov)
  ab[vax == 0] <- 0
  ab[vax == 1 & ab < 0] <- 0
  
  # immune response if everyone is not vaccinated
  ab_vax_0 <- rep(0, n)
  # immune response if everyone get vaccinated
  ab_vax_1 <- numeric(n)
  ab_vax_1[ab_binom == 1] <- rnorm(sum(ab_binom==1), mean = mean_ab_vax_cov, sd = sd_ab_vax_cov)
  ab_vax_1[ab_vax_1 < 0] <- 0
  
  # compute hazard function of covid infection
  # haz_covid_vax_ab_cov <- plogis(covid_rate - 4 * ab - 2.2 * vax + 0.5 * risk + 0.5 * age + 0.1 * race)
  haz_covid_vax0_ab0_cov <- plogis(covid_rate + parm_ab * 0 + parm_vax * 0 + 0.5 * risk + 0.5 * age + 0.1 * race)
  haz_covid_vax1_ab0_cov <- plogis(covid_rate + parm_ab * 0 + parm_vax * 1 + 0.5 * risk + 0.5 * age + 0.1 * race)
  haz_covid_vax0_ab1_cov <- plogis(covid_rate + parm_ab * ab_vax_1 + parm_vax * 0 + 0.5 * risk + 0.5 * age + 0.1 * race)
  haz_covid_vax1_ab1_cov <- plogis(covid_rate + parm_ab * ab_vax_1 + parm_vax * 1 + 0.5 * risk + 0.5 * age + 0.1 * race)
  
  # compute time to covid infection
  time_to_covid_vax0_ab0 <- rgeom(n, haz_covid_vax0_ab0_cov) + 1
  time_to_covid_vax1_ab0 <- rgeom(n, haz_covid_vax1_ab0_cov) + 1
  time_to_covid_vax0_ab1 <- rgeom(n, haz_covid_vax0_ab1_cov) + 1
  time_to_covid_vax1_ab1 <- rgeom(n, haz_covid_vax1_ab1_cov) + 1

  ey11 <- mean(time_to_covid_vax1_ab1 <= t0)
  ey10 <- mean(time_to_covid_vax1_ab0 <= t0)
  ey01 <- mean(time_to_covid_vax0_ab1 <= t0)
  ey00 <- mean(time_to_covid_vax0_ab0 <= t0)
  
  VE <- 1-ey11/ey00
  
  return(VE)
  # return(c(ey11, ey00, ey10, ey01))
}


covid_rate_choices = seq(-10, -1, by = 1)
parm_ab_choices = seq(-6, -0.1, by = 0.1)
parm_vax_choices = seq(-6, -0.1, by = 0.1)
for(parm_ab in parm_ab_choices){
  for(parm_vax in parm_vax_choices){
    
  } 
}
