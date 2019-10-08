data {
  
  int<lower=0> n; // number of years
  int<lower=0> p; // number of parishes
  
  matrix[p, n] births_obs; // parish records for baptims, zero means missing value
  matrix[p, n] deaths_obs; // parish records for burials, zero means missing value
  int<lower=0> n_miss_pop; // number of missing observations in population series
  int<lower=0> n_obs_pop;  // number of observations in population series
  int<lower=0> i_miss_pop[n_miss_pop]; // time points of missing population totals
  int<lower=0> i_obs_pop[n_obs_pop]; // time points of observed population totals
  vector<lower=0>[n_obs_pop] data_obs_pop; // observed population totals
  vector<lower=0>[n] soldiers; // additional deaths (fallen in combat)
 
}
transformed data{
  
  // various helper variables for looping through the observations
  
  vector[n] timeseq; // sequence 1:n needed for vectorization of gamma-coefficents
  int<lower=0> t_miss_b[n] = rep_array(0, n); //how many missing values at each time point
  int<lower=0> t_miss_d[n] = rep_array(0, n); 
  int<lower=0> t_obs_b[n]  = rep_array(0, n); //how many observed values at each time point
  int<lower=0> t_obs_d[n] = rep_array(0, n); 
  int<lower=0> i_miss_b[n] = rep_array(0, n); // indices of time points which contain missing values
  int<lower=0> i_miss_d[n] = rep_array(0, n); // indices of time points which contain missing values
  int<lower=0> cumt_miss_b[n+1]; // cumulative sums of missing values
  int<lower=0> cumt_obs_b[n+1];
  int<lower=0> cumt_miss_d[n+1];
  int<lower=0> cumt_obs_d[n+1];
  int<lower=0> pn_miss_b[p * n] = rep_array(0, p * n); //location of each missing value
  int<lower=0> pn_miss_d[p * n] = rep_array(0, p * n);
  int<lower=0> pn_obs_b[p * n]  = rep_array(0, p * n); //location of each observed value
  int<lower=0> pn_obs_d[p * n]  = rep_array(0, p * n);
  
  int n_miss_b = 0;  // total number of missing time points
  int n_miss_d = 0;  // total number of missing time points
  int n_miss_total_b = 0; // total number of missing values
  int n_miss_total_d = 0;
  int n_obs_total_b = 0; // total number of observed values
  int n_obs_total_d = 0;
  
  for(t in 1:n) {
    timeseq[t] = t;
    t_miss_b[t] = 0;
    t_miss_d[t] = 0;
    t_obs_b[t] = 0;
    t_obs_d[t] = 0;
    for (i in 1:p) {
      if (births_obs[i, t] > 0) {
        t_obs_b[t] += 1; // increase number of observed values at time t
        n_obs_total_b += 1;
        pn_obs_b[n_obs_total_b] = i; // add index of new observed value
      } else {
        t_miss_b[t] += 1; // increase number of missing values at time t
        n_miss_total_b += 1;
        pn_miss_b[n_miss_total_b] = i; // add index of new missing value
      }
      if (deaths_obs[i, t] > 0) {
        t_obs_d[t] += 1; // increase number of observed values at time t
        n_obs_total_d += 1;
        pn_obs_d[n_obs_total_d] = i; // add index of new observed value
      } else {
        t_miss_d[t] += 1;
        n_miss_total_d += 1;
        pn_miss_d[n_miss_total_d] = i;
       
      }
    }
    if(t_miss_b[t] > 0) {
      n_miss_b += 1; // total number of missing time points
      i_miss_b[n_miss_b] = t;
    }
    if(t_miss_d[t] > 0) {
      n_miss_d += 1;
      i_miss_d[n_miss_d] = t;
    }
  }
  cumt_miss_b[1] = 0;
  cumt_obs_b[1] = 0;
  cumt_miss_d[1] = 0;
  cumt_obs_d[1] = 0;
  for(t in 2:(n+1)) {
    cumt_miss_b[t] = cumt_miss_b[t-1] + t_miss_b[t-1];
    cumt_obs_b[t] = cumt_obs_b[t-1] + t_obs_b[t-1];
    cumt_miss_d[t] = cumt_miss_d[t-1] + t_miss_d[t-1];
    cumt_obs_d[t] = cumt_obs_d[t-1] + t_obs_d[t-1];
  }

}

parameters {

  // random walks
  real<lower=0> sigma_drift_b; // sd of random walks
  real<lower=0> sigma_drift_d;
  vector[n] drift_b; // drifts of random walks
  vector[n] drift_d;
  real<lower=0> sigma_rw_b; // sd of random walks
  real<lower=0> sigma_rw_d;
  matrix[p, n] rw_b; 
  matrix[p, n] rw_d;
  
  // observations
  real<lower=0> rate_bd; // rates of the Gamma distributions
  vector<lower=0>[n_miss_b] miss_b; // estimates of total missing records per year 
  vector<lower=0>[n_miss_d] miss_d;
 
  // population estimates of the population mean
  // not including 1697
  vector<lower=0>[n] mu;
  real<lower=0> rate_mu; // noise on mu-process due to crude shape of gammas
  real<lower=0> mu_1647; // initial population size
  
  real<lower=0> sigma_pop;
  real famine; // parish level famine effect, multiplicative increase of exp(famine) to deaths
   // population level famine effect, due to prior knowledge about the drop in population from 1695 to 1697
  real<lower=0, upper=1> famine_ratio;
  simplex[3] lc; // for lambda_1, proportion of lambda_n
  real<lower=0, upper=1> lambda_n; // gamma in 1850
  real<lower=0> r_b; // growth rate of lambda_b
  real<lower=0> r_d; // growth rate of lambda_d
  real<lower=0, upper=1> mid_b_uc;
  real<lower=0, upper=1> mid_d_uc;
}

transformed parameters {
 
  // quality of baptisms is better in 1648
  real<lower=0, upper=1> lbc = lc[1] + lc[2];
  real<lower=0, upper=1> ldc = lc[1];
  real<lower=0,upper=1> l1_b = lbc * lambda_n; // lambda_b in 1648
  real<lower=0,upper=1> l1_d = ldc * lambda_n; // lambda_b in 1648
  real<lower=0, upper=183> mid_b = 183 * mid_b_uc;
  real<lower=0, upper=183> mid_d = 183 * mid_d_uc;
  vector<lower=0>[n] births; // total births (lambda-corrected)
  vector<lower=0>[n] deaths; // total deaths (lambda-corrected)
  vector<lower=0,upper=1>[n] lambda_b; 
  vector<lower=0,upper=1>[n] lambda_d;


  // create lambda curves
  {
     // for births
     vector[n] uc_lambda = 1.0 ./ (1.0 + exp(-r_b * (timeseq - mid_b)));
     real minl = uc_lambda[1];
     real maxl = uc_lambda[n];
     // scale
     lambda_b = (uc_lambda - minl) / (maxl - minl) * (lambda_n - l1_b) + l1_b;
     // and for deaths
     uc_lambda = 1.0 ./ (1.0 + exp(-r_d * (timeseq - mid_d)));
     minl = uc_lambda[1];
     maxl = uc_lambda[n];
     // scale
     // note, we will modify lambda[50] later due to the Great Famine effect
     lambda_d = (uc_lambda - minl) / (maxl - minl) * (lambda_n - l1_d) + l1_d;
  }
  
  for (t in 1:n) {
    births[t] = sum(births_obs[, t]);
    deaths[t] = sum(deaths_obs[, t]);
  }
  
  births[i_miss_b[1:n_miss_b]] += miss_b;
  deaths[i_miss_d[1:n_miss_d]] += miss_d;
  births ./= lambda_b;
  deaths ./= lambda_d;
}

model {
  
  int nmb = 0;
  int nmd = 0;

  // random walk components
  sigma_rw_b ~ gamma(2, 10);
  sigma_rw_d ~ gamma(2, 10);
  sigma_drift_b ~ gamma(2, 20);
  sigma_drift_d ~ gamma(2, 20);
  
  drift_b[1] ~ normal(0, 0.01);
  drift_d[1] ~ normal(0, 0.01);
  rw_b[, 1] ~ normal(3, 1);
  rw_d[, 1] ~ normal(1, 1);
  for(t in 2:n) {
    drift_b[t] ~ normal(drift_b[t-1], sigma_drift_b);
    drift_d[t] ~ normal(drift_d[t-1], sigma_drift_d);
    rw_b[, t] ~ normal(drift_b[t-1] + rw_b[, t - 1], sigma_rw_b);
    rw_d[, t] ~ normal(drift_d[t-1] + rw_d[, t - 1], sigma_rw_d);
  }

  rate_bd ~ gamma(2, 4);
  
  // parish records
  for (t in 1:n) {
    vector[p] exprw = exp(rw_b[, t]);
    if(t_miss_b[t] > 0){
      int m[t_miss_b[t]] = pn_miss_b[(cumt_miss_b[t] + 1):cumt_miss_b[t+1]];
      nmb += 1;
      miss_b[nmb] ~ gamma(rate_bd * sum(exprw[m]), rate_bd);
    }
    if(t_obs_b[t] > 0){
      int o[t_obs_b[t]] = pn_obs_b[(cumt_obs_b[t] + 1):cumt_obs_b[t+1]];
      births_obs[o, t] ~ gamma(rate_bd * exprw[o], rate_bd);
    }
  }

  for (t in 1:49) {
     vector[p] exprw = exp(rw_d[, t]); 
    if(t_miss_d[t] > 0){
      int m[t_miss_d[t]] = pn_miss_d[(cumt_miss_d[t] + 1):cumt_miss_d[t+1]];
      nmd += 1;
      miss_d[nmd] ~ gamma(rate_bd * sum(exprw[m]), rate_bd);
    }
    if(t_obs_d[t] > 0){
      int o[t_obs_d[t]] = pn_obs_d[(cumt_obs_d[t] + 1):cumt_obs_d[t+1]];
      deaths_obs[o, t] ~ gamma(rate_bd * exprw[o], rate_bd);
    }
  }

  // extra famine effect as that is a clear outlier year
  // crude estimate for total burials before and after famine
  // by mean imputation about 5000-10000, in 1697 50000
  // => at least 5-10 times more burials assuming no additional
  // missingness
  // log-link so multiplicative effect (exp(famine)*exp(rw_t))
  // prior for famine as N(2, 0.25)
  famine ~ normal(2, 0.25);
  {
    int t = 50;
    vector[p] exprw = exp(rw_d[, t]); 
    if(t_miss_d[t] > 0){
      int m[t_miss_d[t]] = pn_miss_d[(cumt_miss_d[t] + 1):cumt_miss_d[t+1]];
      nmd += 1;
      miss_d[nmd] ~ gamma(rate_bd * sum(exprw[m] * exp(famine)), rate_bd);
    }
    if(t_obs_d[t] > 0){
      int o[t_obs_d[t]] = pn_obs_d[(cumt_obs_d[t] + 1):cumt_obs_d[t+1]];
      deaths_obs[o, t] ~ gamma(rate_bd * (exprw[o] * exp(famine)), rate_bd);
    }
  }
  for (t in 51:n) {
     vector[p] exprw = exp(rw_d[, t]); 
    if(t_miss_d[t] > 0){
      int m[t_miss_d[t]] = pn_miss_d[(cumt_miss_d[t] + 1):cumt_miss_d[t+1]];
      nmd += 1;
      miss_d[nmd] ~ gamma(rate_bd * sum(exprw[m]), rate_bd);
    }
    if(t_obs_d[t] > 0){
      int o[t_obs_d[t]] = pn_obs_d[(cumt_obs_d[t] + 1):cumt_obs_d[t+1]];
      deaths_obs[o, t] ~ gamma(rate_bd * exprw[o], rate_bd);
    }
  }
  
  // lambda
  
  // means 0.5*ln and 0.75*ln for l1, 
  // 99% prior CIs approx [0.23,0.77] and [0.47, 0.94]
  lc ~ dirichlet([10, 5, 5]'); 
  // 99% prior CI [0.04, 0.2] for growth rate
  r_b ~ gamma(10, 100);
  r_d ~ gamma(10, 100);
  mid_b_uc ~ beta(5,5);
  mid_d_uc ~ beta(5,5);
  // last lambda, mean p/197 where p is number of parishes and 197 is total number
  lambda_n ~ beta(p / 197.0 * 10, (1 - p / 197.0) * 10);
  
  // there should be about 15-30% drop in population from 1695 to 1697
  // prior probability that drop is between 15.5-30.5% is about 99%
  famine_ratio ~ beta(0.775 * 200, (1 - 0.775) * 200);
  
  // 99% prior CI for mu_1647 [358220, 509289]
  mu_1647 ~ gamma(430000*0.0005, 0.0005);
  // 99% prior CI for sd of mu with shape=rate*430,000 [400, 13000] and with  shape=rate*1,500,000 [750,24400]
  rate_mu ~ gamma(2, 4);
  mu[1] ~ gamma(rate_mu * (mu_1647 + births[1] - deaths[1] - soldiers[1]), rate_mu);
  mu[2:49] ~ gamma(rate_mu * (mu[1:48] + births[2:49] - deaths[2:49] - soldiers[2:49]), rate_mu);
  mu[50] ~ gamma(rate_mu * famine_ratio * mu[48], rate_mu);
  mu[51:n] ~ gamma(rate_mu * (mu[50:(n-1)] + births[51:n] - deaths[51:n] - soldiers[51:n]), rate_mu);
  sigma_pop  ~ gamma(5, 0.001);
  data_obs_pop ~ normal(mu[i_obs_pop], sigma_pop);
}

generated quantities {
  real<lower=1647, upper=1850> m_b = 1647 + 183 * mid_b_uc;
  real<lower=1647, upper=1850> m_d = 1647 + 183 * mid_d_uc;
  real famine_correction = (mu[49] + births[50] - soldiers[50] - famine_ratio * mu[48]) / deaths[50];
}
