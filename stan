
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

data {
    int<lower=1> N;//number of rows in our input data, which has a row for each number of days since symptom onset.
    int<lower=1> T_max; // number of days that we would like to predict into the future. 
    int<lower=1> test_n[N]; // total number of pcr tests
    int<lower=0> test_pos[N]; // total number of positive pcr tests
    matrix[N,3] t_ort;
    matrix[T_max, 3] t_new_ort;
}

// the beta terms are the coefficients for the cubic polynomial for time.
parameters{
    real beta_0;
    real beta_1;
    real beta_2;
    real beta_3;
}

transformed parameters{
    vector[N] mu;

    for(i in 1:N){
        mu[i] = beta_0+beta_1*t_ort[i,1]+beta_2*t_ort[i,2]+beta_3*t_ort[i,3]; 
    }
}

model {
    target += binomial_logit_lpmf(test_pos | test_n, mu);
}

// 'sens' is the sensitivity of the RT-PCR over time for the predicted values.
generated quantities{
    vector<lower=0, upper=1>[T_max] sens;
    vector[N] log_lik;

    for(i in 1:T_max){
        sens[i] = inv_logit(beta_0+beta_1*t_new_ort[i,1]+beta_2*t_new_ort[i,2]+beta_3*t_new_ort[i,3]);
    }

    for(i in 1:N){
        log_lik[i] = binomial_logit_lpmf(test_pos[i] | test_n[i], mu[i]);
    }
}
