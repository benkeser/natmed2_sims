#############################################################################
#                               HELPER FUNCTIONS                            #
#############################################################################

#' generate data in a setting similar to what might be expected 
#' in a correlates study of a COVID-19 vaccine
#' 
#' @param n Sample size (30,000 for Moderna)
#' @param covid_rate The intercept in the outcome regression. Make larger to 
#' add more events to mimic what might happen with further follow-up time
#' @return A data.frame with self-explanatory column names

make_ows_data <- function(n = 30000, 
                          covid_rate = -5){
  # baseline covariates
  age <- rbinom(n, 1, 0.4)
  race <- rbinom(n, 1, 0.25)
  risk <- rbinom(n, 1, 0.25)
  
  vax <- rbinom(n, 1, 0.5)
  
  # continuous immune marker truncated at 0
  mean_ab_vax_cov <- 2 - 0.5 * age
  sd_ab_vax_cov <- 1
  ab <- rnorm(n, mean = mean_ab_vax_cov, sd = sd_ab_vax_cov)
  ab[vax == 0] <- 0
  ab[vax == 1 & ab < 0] <- 0
  
  risk_covid_vax_ab_cov <- plogis(covid_rate - 0.5 * ab - 1.8 * vax + 0.7 * risk + 0.2 * age + 0.1 * race)
  covid <- rbinom(n, 1, risk_covid_vax_ab_cov)
  
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
  
  measure_ab = random_subcohort
  
  # sample all cases
  measure_ab[covid == 1] <- 1
  
  # "sample" abs
  ab[measure_ab == 0] <- NA
  
  dat <- data.frame(age = age, race = race, risk = risk,
                    vax = vax, ab = ab, covid = covid,
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
#' @param return the true value of four risk estimators
get_ows_truth1 <- function(n = 1e6, 
                           covid_rate = -5){
  # baseline covariates
  age <- rbinom(n, 1, 0.4)
  race <- rbinom(n, 1, 0.25)
  risk <- rbinom(n, 1, 0.25)
  
  vax <- rbinom(n, 1, 0.5)
  
  # continuous immune marker truncated at 0
  mean_ab_vax_cov <- 2 - 0.5 * age
  sd_ab_vax_cov <- 1
  ab <- rnorm(n, mean = mean_ab_vax_cov, sd = sd_ab_vax_cov)
  ab[vax == 0] <- 0
  ab[vax == 1 & ab < 0] <- 0
  
  ab_vax_0 <- rep(0, n)
  ab_vax_1 <- rnorm(n, mean = mean_ab_vax_cov, sd = sd_ab_vax_cov)
  ab_vax_1[ab_vax_1 < 0] <- 0
  
  risk_covid_vax0_ab0_cov <- plogis(covid_rate - 0.5 * 0 - 1.8 * 0 + 0.7 * risk + 0.2 * age + 0.1 * race)
  risk_covid_vax1_ab0_cov <- plogis(covid_rate - 0.5 * 0 - 1.8 * 1 + 0.7 * risk + 0.2 * age + 0.1 * race)
  risk_covid_vax0_ab1_cov <- plogis(covid_rate - 0.5 * ab_vax_1 - 1.8 * 0 + 0.7 * risk + 0.2 * age + 0.1 * race)
  risk_covid_vax1_ab1_cov <- plogis(covid_rate - 0.5 * ab_vax_1 - 1.8 * 1 + 0.7 * risk + 0.2 * age + 0.1 * race)
  
  ey11 <- mean(risk_covid_vax1_ab1_cov)
  ey10 <- mean(risk_covid_vax1_ab0_cov)
  ey01 <- mean(risk_covid_vax0_ab1_cov)
  ey00 <- mean(risk_covid_vax0_ab0_cov)
  
  return(c(ey11, ey00, ey10, ey01))
}

#' Compute the true value of risk estimators, indirect, direct and total effect
#' based on the data generating process encoded by \code{make_ows_data}.
#' 
#' @param n Sample size (should be very large to have accurate approximation)
#' @param covid_rate The intercept in the outcome regression. Make larger to 
#' add more events to mimic what might happen with further follow-up time
#' @param return a list of the targeted results
get_ows_truth2 <- function(n = 1e6, 
                           covid_rate = -5){
  # baseline covariates
  age <- rbinom(n, 1, 0.4)
  race <- rbinom(n, 1, 0.25)
  risk <- rbinom(n, 1, 0.25)
  
  vax <- rbinom(n, 1, 0.5)
  
  # continuous immune marker truncated at 0
  mean_ab_vax_cov <- 2 - 0.5 * age
  sd_ab_vax_cov <- 1
  ab <- rnorm(n, mean = mean_ab_vax_cov, sd = sd_ab_vax_cov)
  ab[vax == 0] <- 0
  ab[vax == 1 & ab < 0] <- 0
  
  ab_vax_0 <- rep(0, n)
  ab_vax_1 <- rnorm(n, mean = mean_ab_vax_cov, sd = sd_ab_vax_cov)
  ab_vax_1[ab_vax_1 < 0] <- 0
  
  risk_covid_vax0_ab0_cov <- plogis(covid_rate - 0.5 * 0 - 1.8 * 0 + 0.7 * risk + 0.2 * age + 0.1 * race)
  risk_covid_vax1_ab0_cov <- plogis(covid_rate - 0.5 * 0 - 1.8 * 1 + 0.7 * risk + 0.2 * age + 0.1 * race)
  risk_covid_vax0_ab1_cov <- plogis(covid_rate - 0.5 * ab_vax_1 - 1.8 * 0 + 0.7 * risk + 0.2 * age + 0.1 * race)
  risk_covid_vax1_ab1_cov <- plogis(covid_rate - 0.5 * ab_vax_1 - 1.8 * 1 + 0.7 * risk + 0.2 * age + 0.1 * race)
  
  ey11 <- mean(risk_covid_vax1_ab1_cov)
  ey10 <- mean(risk_covid_vax1_ab0_cov)
  ey01 <- mean(risk_covid_vax0_ab1_cov)
  ey00 <- mean(risk_covid_vax0_ab0_cov)
  
  return(list(risk = c(ey11, ey10, ey00),
              indirect = ey11 / ey10, 
              direct = ey10 / ey00,
              total = ey11 / ey00))
}