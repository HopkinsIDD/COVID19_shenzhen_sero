library(tidyverse)
library(rstan)
library(ggpubr)
library(loo)
library(mgcv)

make_analysis_data <- function(dat,
                               stan_model,
                               T_max,
                               poly_est,
                               poly_pred,
                               save_stan=F,
                               ...){
  ## sample from Stan model
  stan_sample <- sampling(stan_model,
                          data=list(N=nrow(dat),
                                    T_max=T_max,
                                    test_n=dat$n_adj,
                                    test_pos=dat$test_pos_adj,
                                    t_ort=as.matrix(poly_est),
                                    t_new_ort=poly_pred
                                    ),
                          ...)
  ## get Stan likelihood
  stan_ll <- suppressWarnings(loo::extract_log_lik(stan_sample) %>% loo::loo())
  
  ## extract parameters
  ## sensitivity (sens) of PCR: P(PCR+ | covid+)
  ## false negative rate (fnr) of PCR: P(PCR- | covid+) = 1 - sens
  sens <- extract(stan_sample, pars="sens")[[1]]
  
  plot_dat <- as_tibble(sens) %>%
    gather("days", "sens") %>%
    mutate(days_since_exposure=gsub(pattern="V", "", days) %>% as.numeric) %>%
    group_by(days_since_exposure) %>%
    summarise(fnr_med=median(1-sens),
              fnr_lb=quantile(1-sens,probs=.025),
              fnr_ub=quantile(1-sens,probs=.975)
              )
  
  if(save_stan){
    return(list(plot_dat=plot_dat,
                stan_ll=stan_ll,
                stan_sample=stan_sample))
  } else{
    return(list(plot_dat=plot_dat,
                stan_ll=stan_ll))
  }
}


make_analysis_data_steve <- function(dat,
                               stan_model,
                               T_max,
                               poly_est,
                               poly_pred,
                               save_stan=F,
                               ...){
  ## sample from Stan model
  stan_sample <- sampling(stan_model,
                          data=list(N=nrow(dat),
                                    J=max(dat$study_idx),
                                    T_max=T_max,
                                    test_n=dat$n_adj,
                                    test_pos=dat$test_pos_adj,
                                    study_idx=dat$study_idx,
                                    t_ort=as.matrix(poly_est),
                                    t_new_ort=poly_pred
                          ),
                          ...)
  ## get Stan likelihood
  stan_ll <- suppressWarnings(loo::extract_log_lik(stan_sample) %>% loo::loo())
  
  ## extract parameters
  ## sensitivity (sens) of PCR: P(PCR+ | covid+)
  ## false negative rate (fnr) of PCR: P(PCR- | covid+) = 1 - sens
  sens <- extract(stan_sample, pars="sens")[[1]]

  plot_dat <- as_tibble(sens) %>%
    gather("days", "sens") %>%
    mutate(days_since_exposure=gsub(pattern="V", "", days) %>% as.numeric) %>%
    group_by(days_since_exposure) %>%
    summarise(fnr_med=median(1-sens),
              fnr_lb=quantile(1-sens,probs=.025),
              fnr_ub=quantile(1-sens,probs=.975)
    )
  
  if(save_stan){
    return(list(plot_dat=plot_dat,
                stan_ll=stan_ll,
                stan_sample=stan_sample))
  } else{
    return(list(plot_dat=plot_dat,
                stan_ll=stan_ll))
  }
}