############################################################
# Datafile S1 - R code for statistical analyses
# File accompanies:
#   Aged soils contribute little to contemporary carbon cycling downstream of thawing permafrost peatlands
#   Andrew J Tanentzap, Katheryn Burd, McKenzie Kuhn, Cristian Estop-Aragonés, Suzanne E. Tank, David Olefeldt
# Please cite if you reuse / find helpful
# File prepared by AJ Tanentzap (ajt65@cam.ac.uk) in Jan 2021

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# The file contains the following sections:
# 0) load packages and create useful functions
# 1) does burning or thaw have stronger have younger C?
# 2) run age attribution model for groundwaters
# 3) run age attribution model in outlets
# 4) fit mixing model to apportion carbon in outlets
# 5) does DOC concentration and age differ between catchments?
# 6) mixing model fitted to consumers
# 7) re-fit mixing model to consumers without 2H
# 8) fit mixing model to algae
# 9) compare DIC and DOC between sites
# Fig 1 – F14C in soils and rivers
# Fig 2 - age-depth profiles of porewater DOC 
# Fig 3 - Age profiles of outlets from both source approportion models (a/b = Raymond, c = Wild)
# Fig 4 – a) consumer ages; b) contributions
# Fig S1 – posterior comparison between consumers and DOC pool
# Fig S2 – isotope biplots
# Fig S3 – histogram of algal 13C compared to Ishikawa dataset  https://link.springer.com/article/10.1007/s11284-018-1623-z
# Fig S4 – compare estimates of resource use with and without H
  
                              '
############################################################
### 0) load data and packages
############################################################
### read data
a14c <- read.csv("./data/atm_14C_series.csv")
o14c <- read.csv("./data/outlet_14C.csv")
p14c <- read.csv("./data/porewater_14C.csv")
f14c <- read.csv("./data/foodweb_14C.csv")
odoc <- read.csv("./data/outlet_DOC.csv")

# create dataframe with Sept DIC data
o14c_DIC <- as.data.frame(cbind(rep(c('burned','unburned'),each=2),rep(c('R1','R2'),2)))
colnames(o14c_DIC) <- c('Site','Sample')
o14c_DIC$d14C <- c(16.29,13.37,-7.65,0.57) 
o14c_DIC$F14C <- c(1.0246,1.0216,1.0004,1.0087)

library(rstan)
library(gplots)
library(rbacon)
library(car)
library(lme4)
library(kdensity)
library(loo)
library(coin) 


############################################################
### 1) does burning or thaw have stronger have younger C?
############################################################
p14c$age <- with(p14c,-8033*log(F14C))
p14c$Site <- relevel(p14c$Site,'UPP')
m1 <- lm(sqrt(age+67) ~ Site, data = p14c)

waterdepth_bpp <- c(40,24,42)
waterdepth_upp <- c(34,24,24)
waterdepth_bog <- c(17.5,19,21)

wilcox_test(c(waterdepth_bog,waterdepth_upp)~as.factor(c(rep('B',3),rep('P',3))))


############################################################
### 2) run age attribution model for groundwaters
############################################################
# find annual mean fraction modern F14C and fit spline
a14c$year2 <- ceiling(a14c$CalibBP)
permod <- with(a14c, tapply(F14C*100,year2,mean))
permod <- as.data.frame(cbind(as.numeric(names(permod)),as.numeric(permod)))
colnames(permod) <- c('year','permod')
ss1 <- with(permod, smooth.spline(year,permod,nknots=1000))

# predict values for last few years from spline
permod2 <- cbind(-60:(1950-2017), predict(ss1,-60:(1950-2017))$y)
colnames(permod2) <- c('year','permod')
permod <- rbind(permod2[nrow(permod2):1,], permod)
permod$t <- as.numeric(rownames(permod))

# concatenate all together and clean-up time-series
permod3 <- as.data.frame(cbind(10000:(1950-2017), predict(ss1,10000:(1950-2017))$y))
permod3 <- permod3[nrow(permod3):1,]
colnames(permod3) <- c('year','permod')
rownames(permod3) <- seq(1:nrow(permod3))
permod3$t <- as.numeric(rownames(permod3))

# fit eq. 1 from Main Text and eq. S1 from Methods S1
kSTAN_sim <- "
 data{
    int<lower=1> N;
    int<lower=1> Nyrs;
    real F14C[N];
    vector[Nyrs] time;
    vector[Nyrs] time_1;
    vector[Nyrs] atmF14C;
   }
 parameters{
    real<lower=0> sdF14C;
    real k;
   }
 transformed parameters{
    real muF14C = dot_product((exp(-k*(time_1))-exp(-k*time)), atmF14C);
   }
 model{
      sdF14C ~ normal(0,10);
      k ~ beta(1,6);
      F14C ~ normal(muF14C, sdF14C);
   }
 generated quantities{
    vector[Nyrs] cont;
     cont = (exp(-k*(time_1))-exp(-k*time)) .*  atmF14C;
     cont = cont/sum(cont)*100;
 }"
its_sim = function(){ list( k = 0.0024, sdF14C = 3 )}
 
# fit for bogs
STANdata_sim_bog <- list(   N = nrow(p14c[p14c$Site=='Bog',]),
                    F14C = p14c[p14c$Site=='Bog',]$F14C*100,
                    time = 1:10000, time_1 = 1:10000-1, Nyrs = 10000,
                    atmF14C = permod3[match((1:10000),permod3$t),'permod']
                  )
k_model_simple_bog <- stan(model_code=kSTAN_sim, data=STANdata_sim_bog, init=its_sim, chains=4, iter=5000, warmup=2500, thin=10, refresh=100, cores=4, open_progress=T, control = list(adapt_delta=0.95))

# fit for peat plateaus - don't split burned/unburned as no difference in ages
STANdata_simple_PP <- list(   N = nrow(p14c[p14c$Site%in%c('BPP','UPP'),]),
                              F14C = p14c[p14c$Site%in%c('BPP','UPP'),]$F14C*100,
                              time = 1:10000, time_1 = 1:10000-1, Nyrs = 10000,
                              atmF14C = permod3[match((1:10000),permod3$t),'permod']
                  )
k_model_simple_PP <- stan(model_code=kSTAN_sim, data=STANdata_simple_PP, init=its_sim, chains=4, iter=5000, warmup=2500, thin=10, refresh=100, cores=4, open_progress=T, control = list(adapt_delta=0.99))

# classify age contributions into bins
prop_cut_bog <- apply(extract(k_model_simple_bog,"cont")[[1]], 1, function(x){ tapply(x, cut(1:10000, breaks = c(0,9,49,99,299,499,999,9999)),sum) })
prop_cut_pp <-  apply(extract(k_model_simple_PP,"cont")[[1]],  1, function(x){ tapply(x, cut(1:10000, breaks = c(0,9,49,99,299,499,999,9999)),sum) })


############################################################
### 3) run age attribution model in outlets (eq. 1 from Main Text)
############################################################
kSTAN <- "
 data{
    int<lower=1> N;
    int<lower=1> Nyrs;
    real F14C[N];
    vector[Nyrs] time;
    vector[Nyrs] time_1;
    vector[Nyrs] atmF14C;
    int site[N];
    int season[N];
   }
 parameters{
    real<lower=0> sdmu;
    real k[N];
    real<lower=0> varb[2];
    real<upper=0> lammy1;
    real<upper=-lammy1> lammy2;
    matrix[2,3] ran_mon;
    real<lower=0,upper=.01> beta_var;
   }
 transformed parameters{
    real muF14C[N];
    real mub[N];
    real<lower=0> alphab[N];
    real<lower=0> betab[N];
    for (i in 1:N){
      mub[i] = inv_logit(lammy1 + lammy2*(site[i]-1) + ran_mon[site[i],season[i]]*varb[site[i]]);
      alphab[i] = ( ((1 - mub[i]) / beta_var) - (1 / mub[i]) ) * pow(mub[i],2);
      betab[i] = alphab[i] * (1 / mub[i] - 1);
    }
    for (i in 1:N)
     muF14C[i] = dot_product((exp(-k[i]*time_1)-exp(-k[i]*time)), atmF14C);
   }
 model{
      lammy1 ~ normal(0,1);
      lammy2 ~ normal(0,1);
      beta_var ~ gamma(1,1);
      varb ~ normal(0,10);
      sdmu ~ normal(0,10);
      to_vector(ran_mon) ~ normal(0,1);
      k ~ beta(alphab,betab);
      F14C ~ normal(muF14C, sdmu);
 }"
o14c$Month <- relevel(o14c$Month,'May')
o14c$Site <- relevel(o14c$Site,'Scott')
STANdata <- list(  N = nrow(o14c),
                    F14C = o14c$F14C*100,
                    time = 1:10000, time_1 = 1:10000-1, Nyrs = 10000,
                    atmF14C = permod3[match((1:10000),permod3$t),'permod'],
                    site = as.numeric(o14c$Site), season = as.numeric(o14c$Month)
                  )
its = function(){ list( k = rep(1e-4,12), betavar = 1e-6, sdmu = 1, lammy1 = -4, lammy2 = 0, varb = rep(1,2), ran_mon = matrix(rep(0,6),nrow=2) )}
k_model <- stan(model_code = kSTAN, data = STANdata, init = its, chains = 4, iter = 5000, warmup = 2500, thin = 10, refresh = 100, cores = 4, open_progress = T)

# find contribution of each year over the past 10k to F14C in each site
cont_s1 <- sapply(1:300,function(x){ (exp(-plogis(extract(k_model,'lammy1')[[1]][101:400][x])*STANdata$time_1)-exp(-plogis(extract(k_model,'lammy1')[[1]][101:400][x])*STANdata$time)) * STANdata$atmF14C })
cont_s1 <- t(apply(cont_s1,2,function(x){x/sum(x)*100}))
cont_s2 <- sapply(1:300,function(x){ (exp(-plogis(extract(k_model,'lammy1')[[1]][101:400][x]+extract(k_model,'lammy2')[[1]][101:400][x])*STANdata$time_1)-exp(-plogis(extract(k_model,'lammy1')[[1]][101:400][x]+extract(k_model,'lammy2')[[1]][101:400][x])*STANdata$time)) * STANdata$atmF14C })
cont_s2 <- t(apply(cont_s2,2,function(x){x/sum(x)*100}))

prop_c_n <- apply(cont_s1, 1, function(x){ tapply(x, cut(1:10000, breaks = c(0,9,49,99,299,499,999,9999)),sum) })
prop_c_s <- apply(cont_s2, 1, function(x){ tapply(x, cut(1:10000, breaks = c(0,9,49,99,299,499,999,9999)),sum) })

            
############################################################
### 4) fit mixing model to apportion carbon in outlets
############################################################
stanMIX1 <- "
data{
   int<lower=1> N;
   real y14C[N];
   real deltaPPmn14C[2];
   real deltaPPsd14C[2];
   real deltaATmn14C[2];
   real deltaATsd14C[2];
   int site[N];
   int season[N];
}
  parameters{
   real<lower=0> sigma[2];
   matrix[3,2] ran_seas;
   real betab[2];
   real<lower=0,upper=0.25> beta_var;
   real<lower=0,upper=1> p_pp[N];
}
transformed parameters{
   vector[N] muy14C;
   vector[N] sdy14C;
   vector<lower=0,upper=1>[N] mu_p;
   vector<lower=0>[N] alpha_p;
   vector<lower=0>[N] beta_p;
   for (i in 1:N){
    mu_p[i] = inv_logit(betab[1] + betab[2]*(site[i]-1) + ran_seas[season[i],1]*sigma[2]);
    alpha_p[i] = (((1 - mu_p[i]) / beta_var) - (1 / mu_p[i])) * pow(mu_p[i],2);
    beta_p[i] = alpha_p[i] * ((1 / mu_p[i]) - 1); 
   }
  for (i in 1:N){
    muy14C[i] = p_pp[i]*deltaPPmn14C[site[i]] + (1-p_pp[i])*deltaATmn14C[site[i]];
    sdy14C[i] = sqrt(pow(p_pp[i]*deltaPPsd14C[site[i]],2) + pow((1-p_pp[i])*deltaATsd14C[site[i]],2) + pow(sigma[1],2));
  }
}
model{
     to_vector(ran_seas) ~ normal(0,1);
     beta_var ~ normal(0,1);
     sigma ~ normal(0,1);
     betab ~ normal(0,10);
     p_pp ~ beta(alpha_p,beta_p);
     y14C ~ normal(muy14C, sdy14C);
}"

# generate atmospheric estimates from extrapolating observed data
# need to convert FMC to d14C so that estimate of error in spline is on same scale
# and we can then compare that with actual measurement error, which is much lower (i.e. we are propagating more uncertainty)
permod_tmp <- with(a14c, tapply((F14C*exp((1950-2017)/8267) - 1)*1000,year2,mean))
permod_tmp <- as.data.frame(cbind(as.numeric(names(permod_tmp)),as.numeric(permod_tmp)))
colnames(permod_tmp) <- c('year','permod')
ss2 <- with(permod_tmp, smooth.spline(year,permod,nknots=1000))

o14c$Month <- relevel(o14c$Month,'July')
o14c$Site <- relevel(o14c$Site,'Noto')
STANdata_a <- list(  N = nrow(o14c),
                   y14C = o14c$d14C, site = as.numeric(o14c$Site), season = as.numeric(o14c$Month),
                   deltaPPmn14C = rev(as.numeric( with(p14c[which(p14c$Site != 'Bog'),], tapply(d14C, Site[drop=T], mean)))),
                   deltaPPsd14C = rev(as.numeric( with(p14c[which(p14c$Site != 'Bog'),], tapply(d14C, Site[drop=T], sd)))),
                   deltaATmn14C = as.vector(tapply(c(f14c[f14c$Type=='peat_permafrost_litter','d14C'],rep(f14c[f14c$Sample=='Bog_Litter','d14C'],2)),
                                                   c(as.character(f14c[f14c$Type=='peat_permafrost_litter','Site']),rep(levels(f14c[f14c$Type=='peat_permafrost_litter','Site']), each = length(f14c[f14c$Sample=='Bog_Litter','d14C']))),
                                                   mean,na.rm=T)),
                   deltaATsd14C = as.vector(tapply(c(f14c[f14c$Type=='peat_permafrost_litter','d14C'],rep(f14c[f14c$Sample=='Bog_Litter','d14C'],2)),
                                                   c(as.character(f14c[f14c$Type=='peat_permafrost_litter','Site']),rep(levels(f14c[f14c$Type=='peat_permafrost_litter','Site']), each = length(f14c[f14c$Sample=='Bog_Litter','d14C']))),
                                                   sd,na.rm=T))
                 )
its = function(){ list( betab = c(1,1), sigma = c(1,1), p_pp = runif(12,0,1), ran_seas = matrix(rep(0,6),nrow=3), beta_var = 1e-4 ) }

# runs using the groundwater DOC values as end members
stanMIX1_a_d <- stan(model_code = stanMIX1, data = STANdata_a, init = its, chains = 4, iter = 8000, warmup = 3000, thin = 20, refresh = 100, cores = 4, open_progress = T, control = list(adapt_delta = 0.95, max_treedepth = 20))

# re-run with solid peat at the bottom of the active layer contributing to outlet DOC
# modify model to allow for active layer contribution to vary among 3 seasons in each of the 2 catchments
stanMIX1_al <- "
data{
   int<lower=1> N;
   real y14C[N];
   real deltaPPmn14C[2,3];
   real deltaPPsd14C[2,3];
   real deltaATmn14C[2];
   real deltaATsd14C[2];
   int site[N];
   int season[N];
}
  parameters{
   real<lower=0> sigma[2];
   real ran_seas[3];
   real betab[2];
   real<lower=0,upper=0.25> beta_var;
   real<lower=0,upper=1> p_pp[N];
}
transformed parameters{
   vector[N] muy14C;
   vector[N] sdy14C;
   vector<lower=0,upper=1>[N] mu_p;
   vector<lower=0>[N] alpha_p;
   vector<lower=0>[N] beta_p;
   for (i in 1:N){
    mu_p[i] = inv_logit(betab[1] + betab[2]*(site[i]-1) + ran_seas[season[i]]*sigma[2]);
    alpha_p[i] = (((1 - mu_p[i]) / beta_var) - (1 / mu_p[i])) * pow(mu_p[i],2);
    beta_p[i] = alpha_p[i] * ((1 / mu_p[i]) - 1);
   }
  for (i in 1:N){
    muy14C[i] = p_pp[i]*deltaPPmn14C[site[i],season[i]] + (1-p_pp[i])*deltaATmn14C[site[i]];
    sdy14C[i] = sqrt(pow(p_pp[i]*deltaPPsd14C[site[i],season[i]],2) + pow((1-p_pp[i])*deltaATsd14C[site[i]],2) + pow(sigma[1],2));
  }
}
model{
     ran_seas ~ normal(0,10);
     beta_var ~ normal(0,1);
     sigma ~ normal(0,1);
     betab ~ normal(0,10);
     p_pp ~ beta(alpha_p,beta_p);
     y14C ~ normal(muy14C, sdy14C);
}"

o14c$Month <- relevel(o14c$Month,'May')
STANdata_a$season = as.numeric(o14c$Month)
STANdata_a$deltaPPmn14C <- rbind(c(-203.6752, -353.3284, -368.1946),c(-192.1408, -222.552, -251.4155))
STANdata_a$deltaPPsd14C <- rbind(c(7.347400, 9.905079, 9.588409),c(10.289166, 6.809884, 8.030875))
its = function(){ list( betab = c(1,1), sigma = c(1,1), p_pp = runif(12,0,1), ran_seas = rep(0,3), beta_var = 1e-4 ) }
stanMIX1_a_l <- stan(model_code = stanMIX1_al, data = STANdata_a, init = its, chains = 4, iter = 8000, warmup = 3000, thin = 20, refresh = 100, cores = 4, open_progress = T, control = list(adapt_delta = 0.95, max_treedepth = 20))

 
############################################################
### 5) does DOC concentration and age differ between catchments?
############################################################
shapiro.test(odoc$DOC_mg_L)
leveneTest(DOC_mg_L ~ Site, data = odoc)
t.test(DOC_mg_L ~ Site, paired=T, var.equal=T, data = odoc)

# fit outlet mixing model
outlet_age_m <- "
 data{
    int<lower=1> N;
    real F14C[N];
    int Site[N];
    int Month[N];
   }
 parameters{
    real<lower=0> sdF14C;
    real<lower=0> month_sd[2];
    matrix[2,3] ranef_m;
    real beta[2];
   }
 transformed parameters{
    real muF14C[N];
    for (i in 1:N)
     muF14C[i] = beta[1] + beta[2]*(Site[i]-1) + ranef_m[Site[i],Month[i]]*month_sd[Site[i]];
   }
 model{
      sdF14C ~ normal(0,1);
      to_vector(ranef_m) ~ normal(0,1);
      month_sd ~ normal(0,1);
      beta ~ normal(0,10);
      F14C ~ normal(muF14C, sdF14C);
   }"
   
o14c$Site <- relevel(o14c$Site,'Scott')
o14c$Month <- relevel(o14c$Month,'May')
outlet_age_STANdata <- list( N = nrow(o14c),F14C = o14c$F14C, Site = as.numeric(o14c$Site), Month = as.numeric(o14c$Month))
it_outlet_age = function(){ list( sdF14C = 0.01, ranef_m = matrix(rep(0.1,6),nrow=2), month_sd = rep(1,2), beta = c(1,0.01))}
m1_outlet_age <- stan(model_code = outlet_age_m, data = outlet_age_STANdata, init = it_outlet_age, chains = 4, iter = 7500, warmup = 5000, thin = 10, refresh = 100, cores = 4, open_progress = T, control = list(adapt_delta = 0.90))


############################################################
### 6) mixing model fitted to consumers
############################################################
stanFWmix <- "data{                                                             // HERE WE DEFINE THE DATA THAT WE INPUT INTO STAN
   int<lower=1> N;                                                              // total number of consumer observations
   int<lower=1> Ncon;                                                           // total number of consumer types
   matrix[N,4] yCNHR;                                                           // matrix of consumer isotope values
   int<lower=1> consumer[N];                                                    // consumer ID
   matrix[4,3] dResourcemu[2];                                                  // matrix of means for each isotope * corresponding resource
   matrix[4,3] dResourcesd[2];                                                  // matrix of means for each isotope * corresponding resource
   real deltaWmn[2];                                                            // mean 2H for water corresponding with each consumer observation
   real deltaWsd[2];                                                            // SD 2H for water corresponding with each consumer observation
   real ptlOmegamu[Ncon];                                                       // mean per-trophic-level consumption of dietary water for each consumer type
   real ptlOmegasd[Ncon];                                                       // SD per-trophic-level consumption of dietary water for each consumer type
   vector[Ncon] tposmu;                                                         // mean trophic level for each consumer type
   real tpossd;                                                                 // SD in trophic level
   int<lower=1> site[N];                                                        // site ID
 }
  parameters{                                                                   // HERE WE DEFINE THE PARAMETERS THAT WE ESTIMATE IN STAN
   real<lower=0> sigma[6];                                                      // estimate 6 SD coefficients and ensure that they are positive
   vector<lower=0>[Ncon] tpos;                                                  // trophic position of each consumer
   vector[Ncon] ptlf;                                                           // per-trophic-level fractionation in N for each consumer
   real<lower=0,upper=1> OmegaDIET[Ncon];                                       // dietary water contribution towards 2H for each consumer
   corr_matrix[4] OmegaCNHR;                                                    // correlation matrix among isotopes
   real betab[4];
   matrix[Ncon,2] ran_eff_c;
   real<lower=1> lambdaA;
   real<lower=0,upper=1> phiA[N]; 
   real<lower=1> lambdaP;
   real<lower=0,upper=1> phiP[N]; 
  }
transformed parameters{                                                         // HERE WE DEFINE THE PARAMETERS THAT WE CALCULATE IN STAN
   simplex[3] p[N];
   vector[N] linregA;
   vector<lower=0>[N] phiAalpha;                                                  
   vector<lower=0>[N] phiAbeta;                                                  
   vector[N] linregP;   
   vector<lower=0>[N] phiPalpha;                                                  
   vector<lower=0>[N] phiPbeta;                                                   
   matrix[N,4] muyCNHR;                                                         // matrix of mean consumer isotope values
   vector[4] sdyCNHR[N];                                                        // SD consumer isotope values
   matrix[4,4] SigmaYCNHR[N];                                                   // variance-covariance matrix of consumer isotope values
   real Omegamu[Ncon];                                                          // mean dietary water contribution towards 2H values
   real Omegavar[Ncon];                                                         // SD dietary water contribution towards 2H values
   real BETAalpha[Ncon];                                                        // alpha parameter for Beta prior distribution of dietary water contribution towards 2H values
   real BETAbeta[Ncon];                                                         // beta parameter for Beta prior distribution of dietary water contribution towards 2H values
   vector[Ncon] ptlfmu = tpos * 2.52;                                           // calculate mean of per-trophic-level N fractionation
   vector[Ncon] ptlfsd;                                                         // sd in per-trophic-level fractionation for each consumer type
     for (i in 1:Ncon)                                                          // calculate SD of per-trophic-level N fractionation
      ptlfsd[i] = sqrt( pow(2.52*tpossd,2) + pow(tpos[i]*1.46,2));
     for (i in 1:Ncon){
      Omegamu[i]   = 1 - pow((1 - ptlOmegamu[i]),tpos[i]);                      // calculate mean contribution of dietary water towards 2H
      Omegavar[i]  = ( pow(tpos[i]*pow((1-ptlOmegamu[i]),tpos[i]-1),2)*         // calculate SD for contribution of dietary water towards 2H
                        pow(ptlOmegasd[i],2) +
                        pow(-pow((1-ptlOmegamu[i]),tpos[i])*
                        log(1-ptlOmegamu[i]),2)*(tpossd*tpossd)
                      );
      BETAalpha[i] = ( ((1-Omegamu[i])/Omegavar[i]) - (1/Omegamu[i]) ) *        // convert mean and variance into alpha parameter of Beta distribution
                        pow(Omegamu[i],2);
      BETAbeta[i]  = BETAalpha[i] * ((1/Omegamu[i]) - 1);                       // convert mean and variance into beta parameter of Beta distribution
     }


     for (i in 1:N)
      linregA[i] = inv_logit(betab[1] + betab[2]*(site[i]-1) + ran_eff_c[consumer[i],1]*sigma[5]);      // eq. 4 in Main Text
     for (i in 1:N)
      linregP[i] = inv_logit(betab[3] + betab[4]*(site[i]-1) + ran_eff_c[consumer[i],2]*sigma[6]);
     phiAalpha = linregA * lambdaA;                                       
     phiAbeta  = (1-linregA) * lambdaA;                                                                //lambdaA here is effective sample size of beta distribution as per https://en.wikipedia.org/wiki/Beta_distribution#Mean_and_sample_size
     phiPalpha = linregP * lambdaP;                                       
     phiPbeta  = (1-linregP) * lambdaP; 
      
     for (i in 1:N){
      p[i,1] = phiA[i];      
      p[i,3] = (1-p[i,1])*phiP[i];
      p[i,2] = (1-p[i,1])*(1-phiP[i]);
     }

     for (i in 1:N){                                                            // calculate  mean and SD for each consumer observation with 4 isotopes and package SDs into Cholesky-decomposed VCV matrix
      muyCNHR[i,1] = dot_product(p[i],dResourcemu[site[i]][1]);
      muyCNHR[i,2] = dot_product(p[i],dResourcemu[site[i]][2]) + ptlf[consumer[i]];
      muyCNHR[i,3] = OmegaDIET[consumer[i]]*deltaWmn[site[i]] + (1-OmegaDIET[consumer[i]])*dot_product(p[i],dResourcemu[site[i]][3]);
      muyCNHR[i,4] = dot_product(p[i],dResourcemu[site[i]][4]);
      sdyCNHR[i,1] = sqrt(pow(p[i,1]*dResourcesd[site[i],1,1],2) + pow(p[i,2]*dResourcesd[site[i],1,2],2) + pow(p[i,3]*dResourcesd[site[i],1,3],2) + pow(sigma[1],2));
      sdyCNHR[i,2] = sqrt(pow(p[i,1]*dResourcesd[site[i],2,1],2) + pow(p[i,2]*dResourcesd[site[i],2,2],2) + pow(p[i,3]*dResourcesd[site[i],2,3],2) + pow(sigma[2],2));
      sdyCNHR[i,3] = sqrt(pow(OmegaDIET[consumer[i]]*deltaWsd[site[i]],2) +
                        pow( (p[i,1]-p[i,1]*OmegaDIET[consumer[i]])*dResourcesd[site[i],3,1],2) +
                        pow( (p[i,2]-p[i,2]*OmegaDIET[consumer[i]])*dResourcesd[site[i],3,2],2) +
                        pow( (p[i,3]-p[i,3]*OmegaDIET[consumer[i]])*dResourcesd[site[i],3,3],2) +
                        pow(sigma[3]*10,2));
      sdyCNHR[i,4] = sqrt(pow(p[i,1]*dResourcesd[site[i],4,1],2) + pow(p[i,2]*dResourcesd[site[i],4,2],2) + pow(p[i,3]*dResourcesd[site[i],4,3],2) + pow(sigma[4]*5,2));
      SigmaYCNHR[i] = quad_form_diag(OmegaCNHR,sdyCNHR[i]);
     }
  }
  model{                                                                        // HERE WE DEFINE THE PRIOR LIKELIHOODS THAT WE ESTIMATE IN STAN
     matrix[4,4] LLCNHR[N];                                                     // Cholesky decomposition of variance-covariance matrix of consumer isotope values
     for (i in 1:N)
      LLCNHR[i] = cholesky_decompose(SigmaYCNHR[i]);
     tpos ~ normal(tposmu, tpossd);
     ptlf ~ normal(ptlfmu, ptlfsd);
     sigma ~ normal(0,1);
     lambdaA ~ pareto(1,2);
     lambdaP ~ pareto(1,2);
     betab ~ normal(0,10);
     to_vector(ran_eff_c) ~ normal(0,10);
     phiA  ~ beta(phiAalpha,phiAbeta);
     phiP  ~ beta(phiPalpha,phiPbeta);      
     OmegaCNHR ~ lkj_corr(2);
     OmegaDIET ~ beta(BETAalpha,BETAbeta);
     for (i in 1:N)
      yCNHR[i] ~ multi_normal_cholesky(muyCNHR[i], LLCNHR[i]);
  }"


STANdata_f <- list(  yCNHR = as.matrix(f14c[f14c$Type %in% c('fish','fish','dragonfly','hyallela','mayfly'),c('d13C','d15N','d2H','d14C')]),
                     dResourcemu = t(as.matrix(apply(f14c[f14c$Type %in% c('oldpeat','modern_veg_litter','peat_permafrost_litter','algae'),c('d13C','d15N','d2H','d14C')],2,function(x){
                                                tapply(x,f14c[f14c$Type %in% c('oldpeat','modern_veg_litter','peat_permafrost_litter','algae'),'Type'][drop=T],mean,na.rm=T) }) )),
                     dResourcesd = t(as.matrix(apply(f14c[f14c$Type %in% c('oldpeat','modern_veg_litter','peat_permafrost_litter','algae'),c('d13C','d15N','d2H','d14C')],2,function(x){
                                                tapply(x,f14c[f14c$Type %in% c('oldpeat','modern_veg_litter','peat_permafrost_litter','algae'),'Type'][drop=T],sd,na.rm=T) }) )),
                     deltaWmn = as.numeric(with(f14c[f14c$Type=='water',],tapply(d2H,Site,mean))), deltaWsd = as.numeric(with(f14c[f14c$Type=='water',],tapply(d2H,Site,sd))), tpossd = 0.1
                 )
STANdata_f$N <- nrow(STANdata_f$yCNHR)
STANdata_f$consumer <- f14c[f14c$Type %in% c('fish','fish','dragonfly','hyallela','mayfly'),'Type']
levels(STANdata_f$consumer)[levels(STANdata_f$consumer)=="mayfly"] <- 'dragonfly'
STANdata_f$tposmu <- levels(STANdata_f$consumer[drop=T])
STANdata_f$ptlOmegamu <- levels(STANdata_f$consumer[drop=T])
STANdata_f$ptlOmegasd <- levels(STANdata_f$consumer[drop=T])
STANdata_f$consumer <- as.numeric(STANdata_f$consumer[drop=T])
STANdata_f$Ncon <- max(STANdata_f$consumer)
STANdata_f$tposmu <- c(2,2.5,1.5)                                          #"aquatic insects"   "fish"    "amphipods"
STANdata_f$ptlOmegamu <- c(0.09,0.12,0.12)                                 #"aquatic insects"   "fish"    "amphipods"
STANdata_f$ptlOmegasd <- c(sqrt(0.5*0.12^2+0.5*0.06^2),0.02,0.12)          #"aquatic insects"   "fish"    "amphipods"
STANdata_f$site <- as.numeric(f14c[f14c$Type %in% c('fish','dragonfly','hyallela','mayfly'),'Site'])

# here split up resources by site
burned_res_df <- f14c[-which( f14c$Type %in% c('oldpeat','modern_veg_litter','peat_permafrost_litter','algae') == F | f14c$Site == 'unburned'),]
unburned_res_df <- f14c[-which( f14c$Type %in% c('oldpeat','modern_veg_litter','peat_permafrost_litter','algae') == F | f14c$Site == 'burned'),]
levels(burned_res_df$Type)[levels(burned_res_df$Type)=='peat_permafrost_litter'] <- 'modern_veg_litter'
levels(unburned_res_df$Type)[levels(unburned_res_df$Type)=='peat_permafrost_litter'] <- 'modern_veg_litter'
dResourcesd_raw <- dResourcemu_raw <- array(NA,c(2,4,3))
dResourcemu_raw[1,,] <- t(as.matrix(apply(burned_res_df[,c('d13C','d15N','d2H','d14C')],2,function(x){tapply(x,burned_res_df$Type[drop=T],mean,na.rm=T)})))
dResourcemu_raw[2,,] <- t(as.matrix(apply(unburned_res_df[,c('d13C','d15N','d2H','d14C')],2,function(x){tapply(x,unburned_res_df$Type[drop=T],mean,na.rm=T)})))
dResourcesd_raw[1,,] <- t(as.matrix(apply(burned_res_df[,c('d13C','d15N','d2H','d14C')],2,function(x){tapply(x,burned_res_df$Type[drop=T],sd,na.rm=T)})))
dResourcesd_raw[2,,] <- t(as.matrix(apply(unburned_res_df[,c('d13C','d15N','d2H','d14C')],2,function(x){tapply(x,unburned_res_df$Type[drop=T],sd,na.rm=T)})))

dResourcemu_raw[1,,2] <- (colMeans(burned_res_df[burned_res_df$Sample %in% c('Bog_Litter','BPP_Litter'),c('d13C','d15N','d2H','d14C')])*48/(48+31)) + 
                          colMeans(burned_res_df[burned_res_df$Sample =='For_Litter',c('d13C','d15N','d2H','d14C')])*31/(48+31)
dResourcemu_raw[2,,2] <- (colMeans(unburned_res_df[unburned_res_df$Sample %in% c('Bog_Litter','UPP_Litter'),c('d13C','d15N','d2H','d14C')])*42/(41+42)) + 
                          colMeans(unburned_res_df[unburned_res_df$Sample =='For_Litter',c('d13C','d15N','d2H','d14C')])*41/(41+42)
dResourcesd_raw[1,,2] <- sqrt( (apply(burned_res_df[burned_res_df$Sample %in% c('Bog_Litter','BPP_Litter'),c('d13C','d15N','d2H','d14C')],2,sd)^2*(48/(48+31))^2) + 
                               (apply(burned_res_df[burned_res_df$Sample =='For_Litter',c('d13C','d15N','d2H','d14C')],2,sd)^2)*(31/(48+31)^2) )
dResourcesd_raw[2,,2] <- sqrt( (apply(unburned_res_df[unburned_res_df$Sample %in% c('Bog_Litter','UPP_Litter'),c('d13C','d15N','d2H','d14C')],2,sd)^2*(48/(48+31))^2) + 
                               (apply(unburned_res_df[unburned_res_df$Sample =='For_Litter',c('d13C','d15N','d2H','d14C')],2,sd)^2)*(31/(48+31)^2) )
STANdata_f$dResourcemu <- dResourcemu_raw
STANdata_f$dResourcesd <- dResourcesd_raw
                             
# set porewater DOC as end-member
STANdata_f$dResourcemu[,4,3] <- rev(as.numeric( with(p14c[which(p14c$Site != 'Bog'),], tapply(d14C, Site[drop=T], mean))))
STANdata_f$dResourcesd[,4,3] <- rev(as.numeric( with(p14c[which(p14c$Site != 'Bog'),], tapply(d14C, Site[drop=T], sd))))
its = function(){ list( sigma = rep(1,6), tpos = c(2,2.5,1.5), ptlf = c(2,2.5,1.5)*2.52, OmegaDIET = c(0.09,0.12,0.12), OmegaCNHR = diag(4),
                        betab = rep(1,4), ran_eff_c = matrix(rep(1,6),3), phiA = rep(0.8,16), phiP = rep(0.5,16), lambdaA = 10, lambdaP = 10 )}
stanFWmix1_p <- stan(model_code=stanFWmix, data=STANdata_f, init=its, chains=4, iter=8000, warmup=3000, thin=20, refresh=100, cores=4, open_progress=T, control = list(adapt_delta=0.95, max_treedepth=14))

# refit using solid peat at base of active layer
burned_res_df <- f14c[-which( f14c$Type %in% c('oldestpeat','modern_veg_litter','peat_permafrost_litter','algae') == F | f14c$Site == 'unburned'),]
unburned_res_df <- f14c[-which( f14c$Type %in% c('oldpeat','modern_veg_litter','peat_permafrost_litter','algae') == F | f14c$Site == 'burned'),]
levels(burned_res_df$Type)[levels(burned_res_df$Type)=='peat_permafrost_litter'] <- 'modern_veg_litter'
levels(unburned_res_df$Type)[levels(unburned_res_df$Type)=='peat_permafrost_litter'] <- 'modern_veg_litter'
dResourcesd_raw <- dResourcemu_raw <- array(NA,c(2,4,3))
dResourcemu_raw[1,,] <- t(as.matrix(apply(burned_res_df[,c('d13C','d15N','d2H','d14C')],2,function(x){tapply(x,burned_res_df$Type[drop=T],mean,na.rm=T)})))
dResourcemu_raw[2,,] <- t(as.matrix(apply(unburned_res_df[,c('d13C','d15N','d2H','d14C')],2,function(x){tapply(x,unburned_res_df$Type[drop=T],mean,na.rm=T)})))
dResourcesd_raw[1,,] <- t(as.matrix(apply(burned_res_df[,c('d13C','d15N','d2H','d14C')],2,function(x){tapply(x,burned_res_df$Type[drop=T],sd,na.rm=T)})))
dResourcesd_raw[2,,] <- t(as.matrix(apply(unburned_res_df[,c('d13C','d15N','d2H','d14C')],2,function(x){tapply(x,unburned_res_df$Type[drop=T],sd,na.rm=T)})))
dResourcemu_raw[1,,2] <- (colMeans(burned_res_df[burned_res_df$Sample %in% c('Bog_Litter','BPP_Litter'),c('d13C','d15N','d2H','d14C')])*48/(48+31)) + 
                          colMeans(burned_res_df[burned_res_df$Sample =='For_Litter',c('d13C','d15N','d2H','d14C')])*31/(48+31)
dResourcemu_raw[2,,2] <- (colMeans(unburned_res_df[unburned_res_df$Sample %in% c('Bog_Litter','UPP_Litter'),c('d13C','d15N','d2H','d14C')])*42/(41+42)) + 
                          colMeans(unburned_res_df[unburned_res_df$Sample =='For_Litter',c('d13C','d15N','d2H','d14C')])*41/(41+42)
dResourcesd_raw[1,,2] <- sqrt( (apply(burned_res_df[burned_res_df$Sample %in% c('Bog_Litter','BPP_Litter'),c('d13C','d15N','d2H','d14C')],2,sd)^2*(48/(48+31))^2) + 
                               (apply(burned_res_df[burned_res_df$Sample =='For_Litter',c('d13C','d15N','d2H','d14C')],2,sd)^2)*(31/(48+31)^2) )
dResourcesd_raw[2,,2] <- sqrt( (apply(unburned_res_df[unburned_res_df$Sample %in% c('Bog_Litter','UPP_Litter'),c('d13C','d15N','d2H','d14C')],2,sd)^2*(48/(48+31))^2) + 
                               (apply(unburned_res_df[unburned_res_df$Sample =='For_Litter',c('d13C','d15N','d2H','d14C')],2,sd)^2)*(31/(48+31)^2) ) 
STANdata_f_al <- STANdata_f                               
STANdata_f_al$dResourcemu <- dResourcemu_raw
STANdata_f_al$dResourcesd <- dResourcesd_raw    
STANdata_f_al$dResourcemu[,4,3] <- c(-368.1946, -251.4155)
STANdata_f_al$dResourcesd[,4,3]<- c(9.588409, 8.030875)
stanFWmix1_p_al <- stan(model_code=stanFWmix, data=STANdata_f_al, init=its, chains=4, iter=8000, warmup=3000, thin=20, refresh=100, cores=4, open_progress=T, control = list(adapt_delta=0.95, max_treedepth=14))


# check it is reasonable to assume 13C at 50 to 60 cm is like that at 27 and 35 cm where DO14C measured
#plot(c(with(f14c[which(f14c$Sample == 'UPP_Litter'),], tapply(d13C, Type[drop=T], mean)),f14c[which(f14c$Type %in% c('oldpeat','oldestpeat')),'d13C'])[order(c(0, f14c[which(f14c$Type %in% c('oldpeat','oldestpeat')),'Rep']))],
#     c(0, f14c[which(f14c$Type %in% c('oldpeat','oldestpeat')),'Rep'])[order(c(0, f14c[which(f14c$Type %in% c('oldpeat','oldestpeat')),'Rep']))],
#     type = 'l', xlab=expression(delta*'13C (\u2030)'), ylab = 'Depth (cm)', ylim=c(140,0), xlim = c(-30,-24), bty='n',las=1)
summary(with(f14c[which(f14c$Type %in% c('oldpeat','oldestpeat')),], lm(d13C~Rep)))
summary(with(f14c[which(f14c$Type %in% c('oldpeat','oldestpeat')),], lm(d15N~Rep)))
summary(with(f14c[which(f14c$Type %in% c('oldpeat','oldestpeat')),], lm(d2H~Rep)))


############################################################
### 7) re-fit mixing model to consumers without 2H
############################################################
stanFWmix_noH <- "data{                                                         // HERE WE DEFINE THE DATA THAT WE INPUT INTO STAN
   int<lower=1> N;                                                              // total number of consumer observations
   int<lower=1> Ncon;                                                           // total number of consumer types
   matrix[N,3] yCNR;                                                            // matrix of consumer isotope values
   int<lower=1> consumer[N];                                                    // consumer ID
   matrix[3,3] dResourcemu[2];                                                  // matrix of means for each isotope * corresponding resource
   matrix[3,3] dResourcesd[2];                                                  // matrix of means for each isotope * corresponding resource
   real ptlOmegamu[Ncon];                                                       // mean per-trophic-level consumption of dietary water for each consumer type
   real ptlOmegasd[Ncon];                                                       // SD per-trophic-level consumption of dietary water for each consumer type
   vector[Ncon] tposmu;                                                         // mean trophic level for each consumer type
   real tpossd;                                                                 // SD in trophic level
   int<lower=1> site[N];                                                        // site ID
 }
  parameters{                                                                   // HERE WE DEFINE THE PARAMETERS THAT WE ESTIMATE IN STAN
   real<lower=0> sigma[5];                                                      
   vector<lower=0>[Ncon] tpos;                                                  // trophic position of each consumer
   vector[Ncon] ptlf;                                                           // per-trophic-level fractionation in N for each consumer
   corr_matrix[3] OmegaCNR;                                                     // correlation matrix among isotopes
   real betab[4];
   matrix[Ncon,2] ran_eff_c;
   real<lower=1> lambdaA;
   real<lower=0,upper=1> phiA[N]; 
   real<lower=1> lambdaP;
   real<lower=0,upper=1> phiP[N]; 
  }
transformed parameters{                                                         // HERE WE DEFINE THE PARAMETERS THAT WE CALCULATE IN STAN
   simplex[3] p[N];
   vector[N] linregA;
   vector<lower=0>[N] phiAalpha;                                                  
   vector<lower=0>[N] phiAbeta;                                                  
   vector[N] linregP;   
   vector<lower=0>[N] phiPalpha;                                                  
   vector<lower=0>[N] phiPbeta;                                                   
   matrix[N,3] muyCNR;                                                          // matrix of mean consumer isotope values
   vector[3] sdyCNR[N];                                                         // SD consumer isotope values
   matrix[3,3] SigmaYCNR[N];                                                    // variance-covariance matrix of consumer isotope values
   matrix[3,3] LLCNR[N];                                                        // Cholesky decomposition of variance-covariance matrix of consumer isotope values
   vector[Ncon] ptlfmu = tpos * 2.52;                                           // calculate mean of per-trophic-level N fractionation
   vector[Ncon] ptlfsd;                                                         // sd in per-trophic-level fractionation for each consumer type
     for (i in 1:Ncon)                                                          // calculate SD of per-trophic-level N fractionation
      ptlfsd[i] = sqrt( pow(2.52*tpossd,2) + pow(tpos[i]*1.46,2));

     for (i in 1:N)
      linregA[i] = inv_logit(betab[1] + betab[2]*(site[i]-1) + ran_eff_c[consumer[i],1]*sigma[4]);
     for (i in 1:N)
      linregP[i] = inv_logit(betab[3] + betab[4]*(site[i]-1) + ran_eff_c[consumer[i],2]*sigma[5]);
     phiAalpha = linregA * lambdaA;                                       
     phiAbeta  = (1-linregA) * lambdaA; 
     phiPalpha = linregP * lambdaP;                                       
     phiPbeta  = (1-linregP) * lambdaP; 
      
     for (i in 1:N){
      p[i,1] = phiA[i];      
      p[i,3] = (1-p[i,1])*phiP[i];
      p[i,2] = (1-p[i,1])*(1-phiP[i]);
     }                        

     for (i in 1:N){                                                            // calculate  mean and SD for each consumer observation with 4 isotopes and package SDs into Cholesky-decomposed VCV matrix
      muyCNR[i,1] = dot_product(p[i],dResourcemu[site[i]][1]);
      muyCNR[i,2] = dot_product(p[i],dResourcemu[site[i]][2]) + ptlf[consumer[i]];
      muyCNR[i,3] = dot_product(p[i],dResourcemu[site[i]][3]);
      sdyCNR[i,1] = sqrt(pow(p[i,1]*dResourcesd[site[i],1,1],2) + pow(p[i,2]*dResourcesd[site[i],1,2],2) + pow(p[i,3]*dResourcesd[site[i],1,3],2) + pow(sigma[1],2));
      sdyCNR[i,2] = sqrt(pow(p[i,1]*dResourcesd[site[i],2,1],2) + pow(p[i,2]*dResourcesd[site[i],2,2],2) + pow(p[i,3]*dResourcesd[site[i],2,3],2) + pow(sigma[2],2));
      sdyCNR[i,3] = sqrt(pow(p[i,1]*dResourcesd[site[i],3,1],2) + pow(p[i,2]*dResourcesd[site[i],3,2],2) + pow(p[i,3]*dResourcesd[site[i],3,3],2) + pow(sigma[3]*5,2));
      SigmaYCNR[i] = quad_form_diag(OmegaCNR,sdyCNR[i]);
      LLCNR[i] = cholesky_decompose(SigmaYCNR[i]);
     }
  }
  model{                                                                        // HERE WE DEFINE THE PRIOR LIKELIHOODS THAT WE ESTIMATE IN STAN
     tpos ~ normal(tposmu, tpossd);
     ptlf ~ normal(ptlfmu, ptlfsd);
     sigma ~ normal(0,1);
     lambdaA ~ pareto(1,2);
     lambdaP ~ pareto(1,2);
     betab ~ normal(0,10);
     to_vector(ran_eff_c) ~ normal(0,10);
     phiA  ~ beta(phiAalpha,phiAbeta);
     phiP  ~ beta(phiPalpha,phiPbeta);      
     OmegaCNR ~ lkj_corr(2);
     for (i in 1:N)
      yCNR[i] ~ multi_normal_cholesky(muyCNR[i], LLCNR[i]);
  }"

STANdata_f_noH <- STANdata_f
STANdata_f_noH$yCNR <- STANdata_f_noH$yCNHR[,-3]
STANdata_f_noH$yCNHR <- NULL
STANdata_f_noH$deltaWmn   <- NULL
STANdata_f_noH$deltaWsd   <- NULL
STANdata_f_noH$dResourcemu <- STANdata_f_noH$dResourcemu[,-3,]
STANdata_f_noH$dResourcesd <- STANdata_f_noH$dResourcesd[,-3,]
 
its = function(){ list( sigma = rep(1,5), tpos = c(2,2.5,1.5), ptlf = c(2,2.5,1.5)*2.52, OmegaCNR = diag(3), betab = rep(1,4), lambdaA = 10, lambdaP = 10, 
                        ran_eff_c = matrix(rep(1,6),3), phiA = rep(0.8,16), phiP = rep(0.5,16) )}
                        
stanFWmix1_noH <- stan(model_code = stanFWmix_noH, data = STANdata_f_noH, init = its, chains = 4, iter = 8000, warmup = 3000, thin = 20, refresh = 100, cores = 4,
                              open_progress = T, control = list(adapt_delta = 0.95, max_treedepth = 14))


############################################################
### 8) fit mixing model to algae
############################################################
# read in fractionation data from Finlay's 2004 L&O paper
alg_frac <- read.csv("finlay_fig7.csv", header = T)
# fit distribution to the data to see what serves it best
library(MASS)
AIC(fitdistr(alg_frac[,1], "lognormal"))
AIC(fitdistr(alg_frac[,1], "normal"))
AIC(fitdistr(alg_frac[,1], "exponential"))
sum(dunif(alg_frac[,1], min(alg_frac[,1]), max(alg_frac[,1]),log=T))*-2+4


stanFWmixA <- "data{                                                            // HERE WE DEFINE THE DATA THAT WE INPUT INTO STAN
   int<lower=1> N;                                                              // total number of observations
   matrix[N,2] yCR;                                                             // matrix of consumer isotope values
   matrix[2,3] dResourcemu[2];                                                     // matrix of means for each isotope * corresponding resource
   matrix[2,3] dResourcesd[2];                                                     // matrix of means for each isotope * corresponding resource
   int<lower=1> site[N];                                                        // site ID
 }  
 transformed data{
    vector[3] alpha;
     for (i in 1:3)
      alpha[i] = 1;
   }
  parameters{                                                                   // HERE WE DEFINE THE PARAMETERS THAT WE ESTIMATE IN STAN
   real<lower=0> sigma[2];                                                      
   corr_matrix[2] OmegaCR; 
   real<lower=0,upper=20> frac;
   simplex[3] p[2];   
  }
transformed parameters{                                                         // HERE WE DEFINE THE PARAMETERS THAT WE CALCULATE IN STAN
   matrix[N,2] muyCR;                                                           // matrix of mean consumer isotope values
   vector[2] sdyCR[N];                                                          // SD consumer isotope values
   matrix[2,2] SigmaYCR[N];                                                     // variance-covariance matrix of consumer isotope values
   matrix[2,2] LLCR[N];                                                         // Cholesky decomposition of variance-covariance matrix of consumer isotope values
     for (i in 1:N){                                                            // calculate the mean and SD for each consumer observation with 2 isotopes and package SDs into Cholesky-decomposed VCV matrix
      muyCR[i,1] = dot_product(p[site[i]],dResourcemu[site[i]][1]) - frac;
      muyCR[i,2] = dot_product(p[site[i]],dResourcemu[site[i]][2]);      
      sdyCR[i,1] = sqrt(pow(p[site[i],1]*dResourcesd[site[i],1,1],2) + pow(p[site[i],2]*dResourcesd[site[i],1,2],2) + pow(p[site[i],3]*dResourcesd[site[i],1,3],2) + pow(sigma[1],2));
      sdyCR[i,2] = sqrt(pow(p[site[i],1]*dResourcesd[site[i],2,1],2) + pow(p[site[i],2]*dResourcesd[site[i],2,2],2) + pow(p[site[i],3]*dResourcesd[site[i],2,3],2) + pow(sigma[2],2));
      SigmaYCR[i] = quad_form_diag(OmegaCR,sdyCR[i]);
      LLCR[i] = cholesky_decompose(SigmaYCR[i]);
     }
  }                                                                     
  model{                                                                    
     sigma ~ normal(0,5);                                        
     for (n in 1:2)
      p[n] ~ dirichlet(alpha);
     OmegaCR ~ lkj_corr(2);
     for (i in 1:N)
      yCR[i] ~ multi_normal_cholesky(muyCR[i], LLCR[i]);
  }
  generated quantities {
  vector[N] log_lik;
  for (i in 1:N)
    log_lik[i] = multi_normal_cholesky_lpdf(yCR[i] | muyCR[i], LLCR[i]);
  }"
                            
f14c_v2 <- f14c
f14c_v2$Type <- factor(f14c_v2$Type, levels = c("modern_soilC",levels(f14c_v2$Type)))
f14c_v2[f14c_v2$Sample=='Bog_Litter','Type'] <- levels(f14c_v2$Type)[1]

STANdata_f_A <- list(  yCR = as.matrix(f14c_v2[f14c_v2$Type =='algae',c('d13C','d14C')]),
                     dResourcemu = t(as.matrix(apply(f14c_v2[f14c_v2$Type %in% c('oldpeat','modern_soilC','water'),c('d13C','d14C')],2,function(x){
                                                tapply(x,f14c_v2[f14c_v2$Type %in% c('oldpeat','modern_soilC','water'),'Type'][drop=T],mean,na.rm=T) }) )),
                     dResourcesd = t(as.matrix(apply(f14c_v2[f14c_v2$Type %in% c('oldpeat','modern_soilC','water'),c('d13C','d14C')],2,function(x){
                                                tapply(x,f14c_v2[f14c_v2$Type %in% c('oldpeat','modern_soilC','water'),'Type'][drop=T],sd,na.rm=T) }) ))
                   )
STANdata_f_A$N <- nrow(STANdata_f_A$yCR)
STANdata_f_A$site <- as.numeric(f14c[f14c$Type =='algae','Site'])

# allow for fractionation of soil OM to soil CO2 by averaging 1-4‰ estimate (https://www.sciencedirect.com/science/article/pii/S0009254199000376?via%3Dihub#BIB37   https://www.nature.com/articles/s41598-017-09049-9)
STANdata_f_A$dResourcemu[1,1:2] <- STANdata_f_A$dResourcemu[1,1:2] + 2.5 
# estimated values from active layer                                           
STANdata_f_A$dResourcemu[2,2] <- mean(c(-295.3064,-237.4975))
STANdata_f_A$dResourcesd[2,2] <- sqrt(.25*(8.456323^2) + .25*(5.126971^2))

# values from Point Barrow http://scrippsco2.ucsd.edu/data/atmospheric_co2/ptb
# and account for fractionation of 7% from dissolution of gaseous CO2 in water (file:///C:/Users/Andrew/Downloads/1994_Pawellek__Veizer_Carbon_cycle_upper_Danube_and_its_tributaries.pdf)
STANdata_f_A$dResourcemu[1,3] <- -8.5 + 7
STANdata_f_A$dResourcesd[1,3] <- 0.382557185
colnames(STANdata_f_A$dResourcemu)[3] <- 'atm_CO2'
STANdata_f_A$dResourcemu[2,3] <- permod3[1,2]
STANdata_f_A$dResourcesd[2,3] <- sqrt(var((ss2$yin - ss2$y)/(1-ss2$lev)))

# values come from Whitcar (1996) review for 13C
colnames(STANdata_f_A$dResourcemu)[2] <- 'geogenic'
STANdata_f_A$dResourcemu[1,2] <- -1.25
STANdata_f_A$dResourcesd[1,2] <- 0.4
# assume rocks Devoniain https://ags.aer.ca/document/Atlas/chapter_1_low_res.pdf
STANdata_f_A$dResourcemu[2,2] <- -1000
STANdata_f_A$dResourcesd[2,2] <- 1
 
# duplicate site resource data
dResourcesd_raw <- dResourcemu_raw <- array(NA,c(2,2,3))
dResourcemu_raw[1,,] <- STANdata_f_A$dResourcemu
dResourcesd_raw[1,,] <- STANdata_f_A$dResourcesd
dResourcemu_raw[2,,] <- STANdata_f_A$dResourcemu
dResourcesd_raw[2,,] <- STANdata_f_A$dResourcesd
STANdata_f_A$dResourcemu <- dResourcemu_raw
STANdata_f_A$dResourcesd <- dResourcesd_raw

# here go and replace values of modern_veg_litter from bog with catchment specific PP litter
STANdata_f_A$dResourcemu[1,1,1] <- mean(f14c_v2[which(f14c_v2$Type == 'peat_permafrost_litter' & f14c_v2$Site =='burned'),'d13C']) + 2.5 
STANdata_f_A$dResourcesd[1,1,1] <- sd(f14c_v2[which(f14c_v2$Type == 'peat_permafrost_litter' & f14c_v2$Site =='burned'),'d13C'])
STANdata_f_A$dResourcemu[2,1,1] <- mean(f14c_v2[which(f14c_v2$Type == 'peat_permafrost_litter' & f14c_v2$Site =='unburned'),'d13C']) + 2.5 
STANdata_f_A$dResourcesd[2,1,1] <- sd(f14c_v2[which(f14c_v2$Type == 'peat_permafrost_litter' & f14c_v2$Site =='unburned'),'d13C'])
STANdata_f_A$dResourcemu[1,2,1] <- mean(f14c_v2[which(f14c_v2$Type == 'peat_permafrost_litter' & f14c_v2$Site =='burned'),'d14C']) 
STANdata_f_A$dResourcesd[1,2,1] <- sd(f14c_v2[which(f14c_v2$Type == 'peat_permafrost_litter' & f14c_v2$Site =='burned'),'d14C'])
STANdata_f_A$dResourcemu[2,2,1] <- mean(f14c_v2[which(f14c_v2$Type == 'peat_permafrost_litter' & f14c_v2$Site =='unburned'),'d14C']) 
STANdata_f_A$dResourcesd[2,2,1] <- sd(f14c_v2[which(f14c_v2$Type == 'peat_permafrost_litter' & f14c_v2$Site =='unburned'),'d14C'])

# fit model
its = function(){ list( sigma = rep(1,2), OmegaCR = diag(2), p = matrix(rep(c(.6,.1,.3),2),nrow=2,byrow=T), frac = 10 )}
stanFWmix1A <- stan(model_code = stanFWmixA, data = STANdata_f_A, init = its, chains = 4, iter = 8000, warmup = 3000, thin = 20, refresh = 100, cores = 4, open_progress = T, control=list(adapt_delta=0.99))
stanFWmix1A_log_lik <- extract_log_lik(stanFWmix1A)
stanFWmix1A_loo <- loo(stanFWmix1A_log_lik,cores=6)
  
         
# compare with scenario where use porewater DOC values
STANdata_f_A_s2 <- STANdata_f_A 
STANdata_f_A_s2$dResourcemu[1:2,2,1] <- rev(as.numeric( with(p14c[which(p14c$Site != 'Bog'),], tapply(d14C, Site[drop=T], mean))))
STANdata_f_A_s2$dResourcesd[1:2,2,1] <- rev(as.numeric( with(p14c[which(p14c$Site != 'Bog'),], tapply(d14C, Site[drop=T], sd))))

stanFWmix1A_s2 <- stan(model_code = stanFWmixA, data = STANdata_f_A_s2, init = its, chains=4, iter=8000, warmup=3000, thin=20, refresh=100, cores=4, open_progress = T, control = list(adapt_delta=0.9))
stanFWmix1A_s2_log_lik <- extract_log_lik(stanFWmix1A_s2)
stanFWmix1A_s2_loo <- loo(stanFWmix1A_s2_log_lik,cores=6)        
                 
                                 
# and scenario where use solid peat
STANdata_f_A_s3 <- STANdata_f_A        
STANdata_f_A_s3$dResourcemu[1,1,1] <- mean(f14c_v2[f14c_v2$Type == 'oldestpeat','d13C']) + 2.5 
STANdata_f_A_s3$dResourcesd[1,1,1] <- sd(f14c_v2[f14c_v2$Type == 'oldestpeat','d13C'])
STANdata_f_A_s3$dResourcemu[2,1,1] <- mean(f14c_v2[f14c_v2$Type == 'oldpeat','d13C']) + 2.5 
STANdata_f_A_s3$dResourcesd[2,1,1] <- sd(f14c_v2[f14c_v2$Type == 'oldpeat','d13C'])

STANdata_f_A_s3$dResourcemu[1,2,1] <- -368.1946
STANdata_f_A_s3$dResourcesd[1,2,1] <- 9.588409
STANdata_f_A_s3$dResourcemu[2,2,1] <- -251.4155
STANdata_f_A_s3$dResourcesd[2,2,1] <- 8.030875
                 
stanFWmix1A_s3 <- stan(model_code = stanFWmixA, data = STANdata_f_A_s3, init = its, chains=4, iter=8000, warmup=3000, thin=20, refresh=100, cores=4, open_progress = T, control = list(adapt_delta=0.9))
stanFWmix1A_s3_log_lik <- extract_log_lik(stanFWmix1A_s3)
stanFWmix1A_s3_loo <- loo(stanFWmix1A_s3_log_lik,cores=6)        
                 
                                            
# compare with scenario where methanotrophy is contributing
STANdata_f_MOB <- STANdata_f_A_s2     
# values come from Whitcar (1999) review
# are consistent with Table 3 here (https://www.biogeosciences-discuss.net/11/16447/2014/bgd-11-16447-2014-print.pdf), if assume mean +/- 3SD = range 
# also consistent with paper says 13C values of respired CO2 from methanotrophs goes to -60 (https://cpb-us-e1.wpmucdn.com/blogs.uoregon.edu/dist/d/3735/files/2013/07/depth_profile_paper-z30i0w.pdf)
STANdata_f_MOB$dResourcemu[1:2,1,1] <- -50
STANdata_f_MOB$dResourcesd[1:2,1,1] <- 5

stanFWmix1A_MOB <- stan(model_code = stanFWmixA, data = STANdata_f_MOB, init = its, chains=4, iter=8000, warmup=3000, thin=20, refresh=100, cores=4, open_progress = T, control = list(adapt_delta=0.9))
stanFWmix1A_MOB_log_lik <- extract_log_lik(stanFWmix1A_MOB)
stanFWmix1A_MOB_loo <- loo(stanFWmix1A_MOB_log_lik,cores=6)           
         
              
# extract model weights and finalise calculations
# want to use pseudo-BMA rather than stacking because we expect some models to have similar predictive pformance and we want them to count uniquely, i.e. don't want those to share weights
# https://cran.r-project.org/web/packages/loo/vignettes/loo2-weights.html 
mod_wt <- loo_model_weights(list(stanFWmix1A_loo,stanFWmix1A_s2_loo,stanFWmix1A_s3_loo,stanFWmix1A_MOB_loo),method="pseudobma")



# are biota older in burned site?
with(f14c[f14c$Type %in% c('fish','dragonfly','algae','hyallela','mayfly'),], leveneTest(d14C ~ Site))
shapiro.test(f14c[f14c$Type %in% c('fish','dragonfly','algae','hyallela','mayfly'),]$d14C)
with(f14c[f14c$Type %in% c('fish','dragonfly','algae','hyallela','mayfly'),], t.test(d14C ~ Site,var.equal=T))
with(f14c[which(f14c$Type=='peat_permafrost_litter'),], t.test(d14C~Site))

# Discussion stats
# approximation of Hutchins et al. 2020 distribution
quantile(rlnorm(1e6,log(0.63),1.2),c(0.025,0.975))
quantile(rtnorm(1e6,0.63,2.5,lower=0),c(0.025,0.975))
quantile(rtnorm(800,0.63,2.5,lower=0)*(extract(stanFWmix1A_s2,'p')[[1]][,,1]*mod_wt[2]+extract(stanFWmix1A_s3,'p')[[1]][,,1]*mod_wt[3]),c(.5,0.025,0.975,1))


############################################################
### 9) compare DIC and DOC between sites
############################################################
shapiro.test(o14c_DIC$d14C)
with(o14c_DIC,leveneTest(d14C~Site))
wilcox.test(o14c[o14c$Month=='Sep',]$d14C,o14c_DIC$d14C[c(3:4,1:2)],paired=T)


############################################################
### Fig 1 – F14C in soils and rivers
############################################################
with(p14c, plot(tapply(F14C,Site,mean), ylab='Fraction modern C (%)', las=1, xaxt='n', xlab='', bty='n', ylim=c(0.7,1.1), xlim=c(1,9), pch=15, cex=1.2) )
axis(1,1:3,lab=c('UPP','Bog','BPP'))
sems_p14c <- with(p14c, tapply(F14C,Site,function(x){sem <- sd(x)/sqrt(length(x)); c(mean(x)+sem,mean(x)-sem) }))
sapply(1:3, function(y){ lines(rep(y,2),sems_p14c[[y]] ) })

points(as.numeric(o14c$Site)+4, o14c$F14C, pch = 19, col = c('gray80','gray50','gray20')[as.numeric(o14c$Month)], cex=0.9 )
points(5:6, with(o14c, tapply(F14C,Site,mean)), pch=19, cex=1.2 )
sems_o14c <- with(o14c, tapply(F14C,Site,function(x){sem <- sd(x)/sqrt(length(x)); c(mean(x)+sem,mean(x)-sem) }))
sapply(1:2, function(y){ lines(rep(y+4,2),sems_o14c[[y]] ) })
axis(1,5:6,lab=c('Unburned','Burned'), las=2)

with(o14c_DIC, points(9:8, tapply(F14C,Site,mean), pch=19, cex=1.2) )
axis(1,8:9,lab=c('Unburned','Burned'), las=2)

sems_o14c_DIC <- with(o14c_DIC, tapply(F14C,Site,function(x){sem <- sd(x)/sqrt(length(x)); c(mean(x)+sem,mean(x)-sem) }))
sapply(1:2, function(y){ lines(rep(c(9:8)[y],2),sems_o14c_DIC[[y]] ) })


############################################################
### Fig 2 - age-depth profiles of porewater DOC 
############################################################
par(mfrow=c(1,3))
barplot2(rowMeans(prop_cut_bog),
                   xlim = c(0, 60), cex.names = 1, plot.ci = TRUE, horiz=T,
                   ci.l = apply(prop_cut_bog,1,quantile,prob=0.025),
                   ci.u = apply(prop_cut_bog,1,quantile,prob=0.975))
              
barplot2(rowMeans(prop_cut_pp),
                   xlim = c(0, 60), cex.names = 1, plot.ci = TRUE, horiz=T,
                   ci.l = apply(prop_cut_pp,1,quantile,prob=0.025),
                   ci.u = apply(prop_cut_pp,1,quantile,prob=0.975))


############################################################
### Fig 3 - Age profiles of outlets from both source approportion models (a/b = Raymond, c = Wild)
############################################################
par(mfrow=c(1,3))

barplot2(rowMeans(prop_c_s),
                   xlim = c(0, 50), cex.names = 1.5, plot.ci = TRUE, horiz=T,
                   ci.l = apply(prop_c_s,1,quantile,prob=0.025),
                   ci.u = apply(prop_c_s,1,quantile,prob=0.975))
barplot2(rowMeans(prop_c_n),
                   xlim = c(0, 50), cex.names = 1.5, plot.ci = TRUE, horiz=T,
                   ci.l = apply(prop_c_n,1,quantile,prob=0.025),
                   ci.u = apply(prop_c_n,1,quantile,prob=0.975))

plot(1:2, c( mean(extract(stanMIX1_a_d, 'p_pp')[[1]][,which(STANdata_a$site==2)])*100,
                 mean(1-extract(stanMIX1_a_d, 'p_pp')[[1]][,which(STANdata_a$site==2)])*100 ),
               xlim = c(0,6), ylim = c(0,100),
               bty = 'n', las = 1, xaxt = 'n',
               cex.axis = 1.3, cex.lab = 1.3,
               pch = 22, bg = 'gray30', cex = 1.2,
               xlab = 'Source', ylab = '% contribution')
sapply(1:2, function(x){lines(rep((1:2)[x],2), rbind(quantile(extract(stanMIX1_a_d, 'p_pp')[[1]][,which(STANdata_a$site==2)],prob=c(0.025,0.975))*100,
                                                         quantile(1-extract(stanMIX1_a_d, 'p_pp')[[1]][,which(STANdata_a$site==2)],prob=c(0.025,0.975))*100)[x,], col = 'gray30') })
points(3:4, c( mean(extract(stanMIX1_a_d, 'p_pp')[[1]][,which(STANdata_a$site==1)])*100,
                 mean(1-extract(stanMIX1_a_d, 'p_pp')[[1]][,which(STANdata_a$site==1)])*100 ), pch = 21, bg = 'white', cex = 1.3)
sapply(1:2, function(x){lines(rep((3:4)[x],2), rbind(quantile(extract(stanMIX1_a_d, 'p_pp')[[1]][,which(STANdata_a$site==1)],prob=c(0.025,0.975))*100,
                                                         quantile(1-extract(stanMIX1_a_d, 'p_pp')[[1]][,which(STANdata_a$site==1)],prob=c(0.025,0.975))*100)[x,], col = 'gray30') })

# add legend and axes
points(0.6, 99, pch = 22, bg = 'gray30', cex = 1.2)
points(0.6, 90, pch = 21, bg = 'white', cex = 1.3)
text(.98, 100, 'Burned', col = 'gray30', cex = 1.1)
text(1.1, 91, 'Unburned', col = 'black', cex = 1.1)
axis(1, at = 1:2, lab = c('Old peat','Surface litterfall'), cex.lab = 1.3)
#dev.off()

# report model fit
quantile(rstantools::bayes_R2(extract(k_model,'muF14C')[[1]],STANdata$F14C))
quantile(rstantools::bayes_R2(extract(stanMIX1_a_l, "muy14C")[[1]],STANdata_a$y14C))
quantile(rstantools::bayes_R2(extract(stanMIX1_a_d, "muy14C")[[1]],STANdata_a$y14C))


############################################################
### Fig 4 – a) consumer ages; b) contributions
############################################################
out_quants <-  sapply(order(STANdata_f$tposmu), function(x) {
                c( quantile(extract(stanFWmix1_p, 'p')[[1]][,which(STANdata_f$consumer == x & STANdata_f$site == 2),3], prob = c(0.5,0.025,0.975)),
                   quantile(extract(stanFWmix1_p, 'p')[[1]][,which(STANdata_f$consumer == x & STANdata_f$site == 1),3], prob = c(0.5,0.025,0.975)),
                   quantile(extract(stanFWmix1_p, 'p')[[1]][,which(STANdata_f$consumer == x & STANdata_f$site == 2),1], prob = c(0.5,0.025,0.975)),
                   quantile(extract(stanFWmix1_p, 'p')[[1]][,which(STANdata_f$consumer == x & STANdata_f$site == 1),1], prob = c(0.5,0.025,0.975)) )
                })
out_quants[1,] <- sapply(order(STANdata_f$tposmu), function(x){mean(extract(stanFWmix1_p, 'p')[[1]][,which(STANdata_f$consumer == x & STANdata_f$site == 2),3])})
out_quants[4,] <- sapply(order(STANdata_f$tposmu), function(x){mean(extract(stanFWmix1_p, 'p')[[1]][,which(STANdata_f$consumer == x & STANdata_f$site == 1),3])})
out_quants[7,] <- sapply(order(STANdata_f$tposmu), function(x){mean(extract(stanFWmix1_p, 'p')[[1]][,which(STANdata_f$consumer == x & STANdata_f$site == 2),1])})
out_quants[10,] <- sapply(order(STANdata_f$tposmu), function(x){mean(extract(stanFWmix1_p, 'p')[[1]][,which(STANdata_f$consumer == x & STANdata_f$site == 1),1])})
       
par(mfrow=c(1,2))
f14c2 <- f14c
levels(f14c2$Type)[levels(f14c2$Type)=="mayfly"] <- 'dragonfly'
levels(f14c2$Type)[levels(f14c2$Type)=="stickleback"] <- 'greyling'
with(f14c2[which(f14c2$Type %in% c('algae','hyallela','dragonfly','greyling')),], boxplot(d14C ~ Site + Type[drop=T], las = 1,
          ylab = expression(Delta * '14C (‰)'),
          lty=1, ylim = c(-40,40)))
axis(4,  at= (exp(c(0,50,200,400)/-8033) * exp((1950-2017)/8267) - 1)*1000, lab = c(0,50,200,400), las = 1)

plot(1:12-0.1, as.numeric(out_quants[seq(1,nrow(out_quants),by=3),])*100,
               xlim = c(0.5,12.5), ylim = c(0,100),
               bty = 'n', las = 1, xaxt = 'n', #log='y',
               cex.axis = 1.3, cex.lab = 1.3, cex = 1.2,
               pch = rep(c(22,21),times=6), bg = rep(c('#FFF9AA','#807815','#B3E095','#387013'),times=4),
               xlab = 'consumer', ylab = '% contribution')
sapply(seq(1,12,by=2), function(x){lines(rep((1:12-0.1)[x],2),c(as.numeric(out_quants[seq(2,nrow(out_quants),3),])[x],as.numeric(out_quants[seq(3,nrow(out_quants),3),])[x])*100, col='#387013') })
sapply(seq(2,12,by=2), function(x){lines(rep((1:12-0.1)[x],2),c(as.numeric(out_quants[seq(2,nrow(out_quants),3),])[x],as.numeric(out_quants[seq(3,nrow(out_quants),3),])[x])*100, col='#807815') })
axis(1, at = seq(2.5,12,by=4), lab = c("amphipods", "aquatic insects", "fish"), cex.lab = 1.3)
sapply(seq(4.5,12,by=4), function(x) { abline(v = x, lty = 3) })
points(0.6, 99, pch = 21, cex = 1.2)
points(0.6, 90, pch = 22, cex = 1.3)
text(2, 100, 'Burned', cex = 1.1)
text(2, 91, 'Unburned',cex = 1.1)

quantile(rstantools::bayes_R2(t(apply(extract(stanFWmix1_p,'muyCNHR')[[1]],1,as.vector)),as.vector(STANdata_f$yCNHR)))
  
    
############################################################
### Fig S1 – posterior comparison between consumers and DOC pool
############################################################
t1 <- extract(stanMIX1_a, 'p_pp')[[1]]
t1 <- extract(stanMIX1_a_d, 'p_pp')[[1]]
t2 <- extract(stanFWmix1, 'p')[[1]][,,3]
t2 <- extract(stanFWmix1_p, 'p')[[1]][,,3]                 
quantile(sapply(1:800, function(x){ as.numeric(ks.test(t1[x,],t2[x,])$p.value) }),c(.5, .025,.975))
plot(kdensity(t1[1,],kernel='beta'), ylim = c(0,60),xlim=c(0,0.4),col=adjustcolor('blue', alpha.f = 0.1), main='',las=1,xlab='contribution')
apply(t1[-1,], 1, function(x) { lines(kdensity(x,kernel='beta'),col=adjustcolor('blue', alpha.f = 0.1)) })
apply(t2, 1, function(x) { lines(kdensity(x,kernel='beta'),col=adjustcolor('red', alpha.f = 0.1)) })


############################################################
### Fig S2 – isotope biplots
############################################################
par(mfrow=c(3,2), mar = c(4.5,4.5,.5,.5), oma = c(1,1,1,1))
#plot(STANdata_f$yCNHR[,c(1,4)], pch = 19, xlim = c(-60,0), ylim = c(-1000,100))
colnames(STANdata_f$yCNHR)[4] <- 'D14C'
plot(STANdata_f$yCNHR[,c(1,4)], pch = c(19,17)[STANdata_f$site], xlim = c(-60,-25), ylim = c(-350,100), cex=1.4, col=rgb(0,0,0,.3), cex.lab=1.1)
plotCI(STANdata_f$dResourcemu[1,1,], STANdata_f$dResourcemu[1,4,], col = c('lightgreen','darkgreen','tan4'), err='y', uiw=STANdata_f$dResourcesd[1,4,], pch=19, add=T, cex=1.5)
plotCI(STANdata_f$dResourcemu[1,1,], STANdata_f$dResourcemu[1,4,], col = c('lightgreen','darkgreen','tan4'), err='x', uiw=STANdata_f$dResourcesd[1,1,], pch=NA, add=T)
plotCI(STANdata_f$dResourcemu[2,1,], STANdata_f$dResourcemu[2,4,], col = c('lightgreen','darkgreen','tan4'), err='y', uiw=STANdata_f$dResourcesd[2,4,], pch=17, add=T, cex=1.5)
plotCI(STANdata_f$dResourcemu[2,1,], STANdata_f$dResourcemu[2,4,], col = c('lightgreen','darkgreen','tan4'), err='x', uiw=STANdata_f$dResourcesd[2,1,], pch=NA, add=T)
# add MOB
plotCI(STANdata_f_MOB$dResourcemu[1,1,1],STANdata_f_MOB$dResourcemu[1,2,1], err='y', uiw=STANdata_f_MOB$dResourcesd[1,2,1], col='red', pch=19, add=T)
plotCI(STANdata_f_MOB$dResourcemu[1,1,1],STANdata_f_MOB$dResourcemu[1,2,1], err='x', uiw=STANdata_f_MOB$dResourcesd[1,1,1], col='red', pch=NA, add=T)
plotCI(STANdata_f_MOB$dResourcemu[2,1,1],STANdata_f_MOB$dResourcemu[2,2,1], err='y', uiw=STANdata_f_MOB$dResourcesd[2,2,1], col='red', pch=17, add=T)
plotCI(STANdata_f_MOB$dResourcemu[2,1,1],STANdata_f_MOB$dResourcemu[2,2,1], err='x', uiw=STANdata_f_MOB$dResourcesd[2,1,1], col='red', pch=NA, add=T)

# add ROCKS
#plotCI(STANdata_f_rocks$dResourcemu[1,2],STANdata_f_rocks$dResourcemu[2,2], err='y', uiw=STANdata_f_rocks$dResourcesd[2,2], col='gray80', pch=19, add=T)
#plotCI(STANdata_f_rocks$dResourcemu[1,2],STANdata_f_rocks$dResourcemu[2,2], err='x', uiw=STANdata_f_rocks$dResourcesd[1,2], col='gray80', pch=19, add=T)

plot(STANdata_f$yCNHR[,3:4], pch = c(19,17)[STANdata_f$site], xlim = c(-300,-150), ylim = c(-350,100), cex=1.4, col=rgb(0,0,0,.3), cex.lab=1.1)
plotCI(STANdata_f$dResourcemu[1,3,], STANdata_f$dResourcemu[1,4,], col = c('lightgreen','darkgreen','tan4'), err='y', uiw=STANdata_f$dResourcesd[1,4,], pch=19, add=T, cex=1.5)
plotCI(STANdata_f$dResourcemu[1,3,], STANdata_f$dResourcemu[1,4,], col = c('lightgreen','darkgreen','tan4'), err='x', uiw=STANdata_f$dResourcesd[1,3,], pch=NA, add=T)
plotCI(STANdata_f$dResourcemu[2,3,], STANdata_f$dResourcemu[2,4,], col = c('lightgreen','darkgreen','tan4'), err='y', uiw=STANdata_f$dResourcesd[2,4,], pch=17, add=T, cex=1.5)
plotCI(STANdata_f$dResourcemu[2,3,], STANdata_f$dResourcemu[2,4,], col = c('lightgreen','darkgreen','tan4'), err='x', uiw=STANdata_f$dResourcesd[2,3,], pch=NA, add=T)

plot(STANdata_f$yCNHR[,1:2], pch = c(19,17)[STANdata_f$site], xlim = c(-60,-25), ylim = c(-5,10), cex=1.4, col=rgb(0,0,0,.3), cex.lab=1.1)
plotCI(STANdata_f$dResourcemu[1,1,], STANdata_f$dResourcemu[1,2,], col = c('lightgreen','darkgreen','tan4'), err='y', uiw=STANdata_f$dResourcesd[1,2,], pch=19, add=T, cex=1.5)
plotCI(STANdata_f$dResourcemu[1,1,], STANdata_f$dResourcemu[1,2,], col = c('lightgreen','darkgreen','tan4'), err='x', uiw=STANdata_f$dResourcesd[1,1,], pch=NA, add=T)
plotCI(STANdata_f$dResourcemu[2,1,], STANdata_f$dResourcemu[2,2,], col = c('lightgreen','darkgreen','tan4'), err='y', uiw=STANdata_f$dResourcesd[2,2,], pch=17, add=T, cex=1.5)
plotCI(STANdata_f$dResourcemu[2,1,], STANdata_f$dResourcemu[2,2,], col = c('lightgreen','darkgreen','tan4'), err='x', uiw=STANdata_f$dResourcesd[2,1,], pch=NA, add=T)

plot(STANdata_f$yCNHR[,3:2], pch = c(19,17)[STANdata_f$site], xlim = c(-300,-150), ylim = c(-5,10), cex=1.4, col=rgb(0,0,0,.3), cex.lab=1.1)
plotCI(STANdata_f$dResourcemu[1,3,], STANdata_f$dResourcemu[1,2,], col = c('lightgreen','darkgreen','tan4'), err='y', uiw=STANdata_f$dResourcesd[1,2,], pch=19, add=T, cex=1.5)
plotCI(STANdata_f$dResourcemu[1,3,], STANdata_f$dResourcemu[1,2,], col = c('lightgreen','darkgreen','tan4'), err='x', uiw=STANdata_f$dResourcesd[1,3,], pch=NA, add=T)
plotCI(STANdata_f$dResourcemu[2,3,], STANdata_f$dResourcemu[2,2,], col = c('lightgreen','darkgreen','tan4'), err='y', uiw=STANdata_f$dResourcesd[2,2,], pch=17, add=T, cex=1.5)
plotCI(STANdata_f$dResourcemu[2,3,], STANdata_f$dResourcemu[2,2,], col = c('lightgreen','darkgreen','tan4'), err='x', uiw=STANdata_f$dResourcesd[2,3,], pch=NA, add=T)

plot(STANdata_f$yCNHR[,c(4,2)], pch = c(19,17)[STANdata_f$site], xlim = c(-350,100), ylim = c(-5,10), cex=1.4, col=rgb(0,0,0,.3), cex.lab=1.1)
plotCI(STANdata_f$dResourcemu[1,4,], STANdata_f$dResourcemu[1,2,], col = c('lightgreen','darkgreen','tan4'), err='y', uiw=STANdata_f$dResourcesd[1,2,], pch=19, add=T, cex=1.5)
plotCI(STANdata_f$dResourcemu[1,4,], STANdata_f$dResourcemu[1,2,], col = c('lightgreen','darkgreen','tan4'), err='x', uiw=STANdata_f$dResourcesd[1,4,], pch=NA, add=T)
plotCI(STANdata_f$dResourcemu[2,4,], STANdata_f$dResourcemu[2,2,], col = c('lightgreen','darkgreen','tan4'), err='y', uiw=STANdata_f$dResourcesd[2,2,], pch=17, add=T, cex=1.5)
plotCI(STANdata_f$dResourcemu[2,4,], STANdata_f$dResourcemu[2,2,], col = c('lightgreen','darkgreen','tan4'), err='x', uiw=STANdata_f$dResourcesd[2,4,], pch=NA, add=T)

plot(0,0,axes=F,xlab='',ylab='',bty='n',pch=NA,xlim=c(0,10),ylim=c(0,10))
points( rep(0,5), c(9,7,5,3,1), col = c('lightgreen','darkgreen','tan4','red','gray70'), pch=19, cex=1.8)
sapply(1:5, function(x){ text(1,c(9,7,5,3,1)[x],c('litter','algae','porewater DOC','methanotroph CO2','consumers')[x],cex=1.2,adj=0) })
points( rep(5,2), c(9,7), pch=c(21,24), cex=1.8,bg='white')
sapply(1:2, function(x){ text(6,c(9,7)[x],c('burned','unburned')[x],cex=1.2,adj=0) })


############################################################
### Fig S3 – histogram of algal 13C compared to Ishikawa dataset  https://link.springer.com/article/10.1007/s11284-018-1623-z
############################################################
globperi <- read.table('global_SI_periphyton.csv',header=T,sep=',')
hist(globperi[, ]$d13C_permil,main='',las=1,xlab=expression(Delta13C*' (%)'),xlim=c(-50,0))  #which(globperi$dominant_taxa %in% c('macrophytes','bryophyte') == F) makes no difference
par(new=TRUE)
hist(f14c[f14c$Type=='algae','d13C'],xlim=c(-50,0),col='black', ylim = c(0,30),axes=F,main='',xlab='',ylab='')

leveneTest(c(f14c[f14c$Type=='algae','d13C'],globperi$d13C_permil[!is.na(globperi$d13C_permil)]),
           c(rep(1,sum(f14c$Type=='algae')),rep(2,sum(globperi$d13C_permil != 0,na.rm=T))))
t.test(f14c[f14c$Type=='algae','d13C'],globperi$d13C_permil[!is.na(globperi$d13C_permil)],var.equal=T)


############################################################
### Fig S4 – compare estimates of resource use with and without H
############################################################
cor(as.vector(sapply(1:3, function(x){ colMeans(extract(stanFWmix1_p, 'p')[[1]][,,x]) })), as.vector(sapply(1:3, function(x){ colMeans(extract(stanFWmix1_noH, 'p')[[1]][,,x]) })), method = "spearman")

par(mfrow=c(1,3))
plot(colMeans(extract(stanFWmix1_p, 'p')[[1]][,,1]), colMeans(extract(stanFWmix1_noH, 'p')[[1]][,,1]), ylab='p_algae (without H)', xlab='p_algae (with H)', las = 1, xlim = c(0.6,1), ylim = c(0.6,1), pch=19); abline(0,1)
plot(colMeans(extract(stanFWmix1_p, 'p')[[1]][,,2]), colMeans(extract(stanFWmix1_noH, 'p')[[1]][,,2]), ylab='p_litter (without H)', xlab='p_litter (with H)', las = 1, xlim = c(0,.4), ylim = c(0,.4), pch=19); abline(0,1)
plot(colMeans(extract(stanFWmix1_p, 'p')[[1]][,,3]), colMeans(extract(stanFWmix1_noH, 'p')[[1]][,,3]), ylab='p_soil_DOC (without H)', xlab='p_soil_DOC (with H)', las = 1, xlim = c(0,0.08), ylim = c(0,0.08), pch=19); abline(0,1)