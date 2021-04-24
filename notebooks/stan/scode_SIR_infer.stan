
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
    int I[nobs];
}

transformed data {
  real x_r[0];
  int x_i[0];
  }

parameters{
    real<lower = 0> theta[2]; // model parameters 
    real<lower = 0, upper = 1> S0;  // initial fraction of susceptible individuals
    real phi;
  }

transformed parameters{
    real ln[nobs,3]; // solution from the ODE solver
    real ln_init[3];     // initial conditions for fractions

    ln_init[1] = S0; //S, ODE system is Nsusc normalized
    ln_init[2] = 1-S0;//1/Nsusc; //I
    ln_init[3] = 0; //R

    ln = integrate_ode_rk45(SIR, ln_init, t[1]-1, t, theta, x_r, x_i);
  }

model{
    //priors
    S0 ~ beta(2, 2); //some prior for between 0 and 1 fraction of the population
    for (i in 1:2){
        theta[i] ~ lognormal(0,1);
    }

    //likelihood
    for (i in 1:nobs){
        # target += poisson_lpmf(I[i]|ln[i,2]*Nsusc);  # I
        target += neg_binomial_2_lpmf(I[i]|ln[i, 2] * Nsusc, phi);
        }
  }

