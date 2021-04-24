
functions{
  real[] SIR(real t,  // time
            real[] ln_init,// system state {infected,cases,susceptible}
            real[] theta,
            real[] x_r,
            int[] x_i)

    {
        real dy_dt[3];

        real kI = theta[1];
        real kR = theta[2];
        real lS = ln_init[1];  # susceptible
        real lI = ln_init[2];  # infected
        real lR = ln_init[3];  # recovered

        dy_dt[1] = -kI*lS*lI; //dS
        dy_dt[2] = kI*lS*lI - kR*lI; //dI
        dy_dt[3] = kR*lR; //dR
        return dy_dt;
    }
}

data {
    int<lower=1> Nsusc; //number that are susceptible (the effective population of interest)
    int<lower = 1> nobs; // number of days observed
    real t[nobs];
    real kI;
    real kR;
    real<lower=0, upper=1> S0; //normalized 0-1
    real phi; //inflation factor for Poisson variance
}
transformed data {
  real x_r[0];
  int x_i[0];
  }
generated quantities{
  real theta[2]; 
  real ln_init[3]; // 
  real ln[nobs,3]; // Nsusc normalized rates
  real I[nobs];
 
  theta[1] = kI;
  theta[2] = kR;
  
  ln_init[1] = S0;//S0, initial susceptible, S>0 and I>0, and S + I = 1 
  ln_init[2] = 1-S0;//I0, initial infected
  ln_init[3] = 0.0;//R0, no one recovered
  
  ln = integrate_ode_rk45(SIR, ln_init, t[1]-1, t, theta, x_r, x_i);
  
#   for (i in 1:nobs) {
#     I[i] = poisson_rng(ln[i, 2] * Nsusc); // lnS + lnI = 1 (norm fractions of pop), need to rescale
#   }

  for (i in 1:nobs) {
    I[i] = neg_binomial_2_rng(ln[i, 2] * Nsusc, phi); // lnS + lnI = 1 (norm fractions of pop), need to rescale
  }
    }
