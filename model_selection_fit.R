### Trajectory Project: model selection and goodness of fit
### Chenyang Gu
### Jul.2, 2014

library(MASS)
library(coda)
library(doBy)


## 1. read data 
## 2. read MCMC output, select one chain


mcmc_output = trajectory.chains[[1]]

Nsim = nrow(mcmc_output)





## facility-level random effects, d.est
# standard deviations
sigma.d1 = mcmc_output[,9]
sigma.d2 = mcmc_output[,10]
sigma.d3 = mcmc_output[,11]
sigma.d4 = mcmc_output[,12]

# correlations
rho.d11 = mcmc_output[,13]
rho.d21 = mcmc_output[,14]
rho.d31 = mcmc_output[,15]
rho.d41 = mcmc_output[,16]

rho.d12 = mcmc_output[,17]
rho.d22 = mcmc_output[,18]
rho.d32 = mcmc_output[,19]
rho.d42 = mcmc_output[,20]

rho.d13 = mcmc_output[,21]
rho.d23 = mcmc_output[,22]
rho.d33 = mcmc_output[,23]
rho.d43 = mcmc_output[,24]

rho.d14 = mcmc_output[,25]
rho.d24 = mcmc_output[,26]
rho.d34 = mcmc_output[,27]
rho.d44 = mcmc_output[,28]

# covariances
sigma.d21 = rho.d21*sigma.d1*sigma.d2
sigma.d31 = rho.d31*sigma.d1*sigma.d3
sigma.d41 = rho.d41*sigma.d1*sigma.d4

sigma.d12 = rho.d12*sigma.d2*sigma.d1
sigma.d32 = rho.d32*sigma.d2*sigma.d3
sigma.d42 = rho.d42*sigma.d2*sigma.d4

sigma.d13 = rho.d13*sigma.d3*sigma.d1
sigma.d23 = rho.d23*sigma.d3*sigma.d2
sigma.d43 = rho.d43*sigma.d3*sigma.d4

sigma.d14 = rho.d14*sigma.d4*sigma.d1
sigma.d24 = rho.d24*sigma.d4*sigma.d2
sigma.d34 = rho.d34*sigma.d4*sigma.d3

# combine MCMC outputs for facility-level variance parameters
Sigma.d = cbind(sigma.d1^2,sigma.d12, sigma.d13, sigma.d14, 
                sigma.d21, sigma.d2^2,sigma.d23, sigma.d24, 
                sigma.d31, sigma.d32, sigma.d3^2,sigma.d34, 
                sigma.d41, sigma.d42, sigma.d43, sigma.d4^2)
#dim(Sigma.d)
#[1] 200  16

d.samp = function(x) {
	Sigma = matrix(x, nrow = 4, ncol = 4, byrow = TRUE)
	d = mvrnorm(nFac, mu = rep(0,4), Sigma = Sigma)
	return(d)
}

# estimated facility-level random effects
# save into array 
set.seed(0704)
D.est = array(NA, c(Nsim,nFac,4))
for(i in 1:Nsim){
	D.est[i,,] = d.samp(Sigma.d[i,])
}
#dim(D.est)
#[1]  200 2659    4

# Extend facility-level random effects with facility index "facility"
d1.est = D.est[ , facility, 1]
d2.est = D.est[ , facility, 2]
d3.est = D.est[ , facility, 3]
d4.est = D.est[ , facility, 4]




# effect of mixture component
g1 = mcmc_output[,33:34]
g2 = mcmc_output[,35:36]
g3 = mcmc_output[,37:38]
g4 = mcmc_output[,39:40]

C = mcmc_output[,59:2717]
#dim(C)
#[1]  200 2659

g1.est = g2.est = g3.est = g4.est = matrix(0, nrow = Nsim, ncol = ncol(C))
for(i in 1:nrow(C)){
	g1.est[i,] = g1[i,C[i,]]
	g2.est[i,] = g2[i,C[i,]]
	g3.est[i,] = g3[i,C[i,]]
	g4.est[i,] = g4[i,C[i,]] 
}
#dim(g1.est)
#[1]  200 2659

# estimated mean of mixture component
g1.est2 = g1.est[, facility]
g2.est2 = g2.est[, facility]
g3.est2 = g3.est[, facility]
g4.est2 = g4.est[, facility]

#dim(g1.est2)
#[1]  200 3106



# fixed-effects
beta11 = mcmc_output[,41]
beta12 = mcmc_output[,42]
beta13 = mcmc_output[,43]
beta14 = mcmc_output[,44]
beta15 = mcmc_output[,45]
beta16 = mcmc_output[,46]

beta21 = mcmc_output[,47]
beta22 = mcmc_output[,48]
beta23 = mcmc_output[,49]
beta24 = mcmc_output[,50]

beta31 = mcmc_output[,51]
beta32 = mcmc_output[,52]
beta33 = mcmc_output[,53]
beta34 = mcmc_output[,54]
beta35 = mcmc_output[,55]

beta41 = mcmc_output[,56]
beta42 = mcmc_output[,57]
beta43 = mcmc_output[,58]

Beta1 = mcmc_output[,41:46]
Beta2 = mcmc_output[,47:50]
Beta3 = mcmc_output[,51:55]
Beta4 = mcmc_output[,56:58]

fix.eff1 = Beta1 %*% t(X1)
fix.eff2 = Beta2 %*% t(X2)
fix.eff3 = Beta3 %*% t(X3)
fix.eff4 = Beta4 %*% t(X4)

#dim(fix.eff1)
#[1]  200 3106


## calculate mu.B1 and mu.B2
mu.B11 = g1.est2 + d1.est + fix.eff1
mu.B12 = g2.est2 + d2.est + fix.eff2
mu.B21 = g3.est2 + d3.est + fix.eff3
mu.B22 = g4.est2 + d4.est + fix.eff4




## subject-level random effects
# pre- and post-hospitalization random effect covariance matrix
sigma.a1 = mcmc_output[,3]
sigma.b1 = mcmc_output[,4]
rho1     = mcmc_output[,5]

sigma.a2 = mcmc_output[,6]
sigma.b2 = mcmc_output[,7]
rho2     = mcmc_output[,8]

sigma.cov1 = rho1*sigma.a1*sigma.b1
sigma.cov2 = rho2*sigma.a2*sigma.b2

Sigma.ab1.comb = cbind(sigma.a1^2, sigma.cov1, sigma.cov1, sigma.b1^2)
Sigma.ab2.comb = cbind(sigma.a2^2, sigma.cov2, sigma.cov2, sigma.b2^2)



ab.samp = function(mu.vec1, mu.vec2, sigma.vec) {
	mu.vec = cbind(mu.vec1, mu.vec2)
	Sigma = matrix(sigma.vec, nrow = 2, ncol = 2, byrow = TRUE)
	ab = matrix(0, nrow = nrow(mu.vec), ncol = 2)
	for(i in 1:nrow(mu.vec)){
		ab[i,] = mvrnorm(1, mu.vec[i,], Sigma = Sigma)
	}
	return(ab)
}


# pre-hospitalization
set.seed(0704)
ab1.est = array(NA, c(Nsim,nSubject,2))
for(i in 1:Nsim){
	ab1.est[i,,] = ab.samp(mu.B11[i,], mu.B12[i,], Sigma.ab1.comb[i,])
}

# post-hospitalization
set.seed(0704)
ab2.est = array(NA, c(Nsim,nSubject,2))
for(i in 1:Nsim){
	ab2.est[i,,] = ab.samp(mu.B21[i,], mu.B22[i,], Sigma.ab2.comb[i,])
}

#dim(ab1.est)
#[1]  200 3106    2

a1.est = ab1.est[,subject1,1]
b1.est = ab1.est[,subject1,2]
a2.est = ab2.est[,subject2,1]
b2.est = ab2.est[,subject2,2]


## measuremen-level 
sigma.y1 = mcmc_output[,1]
sigma.y2 = mcmc_output[,2]


mu1 = matrix(0, nrow = Nsim, ncol = length(time1))
mu2 = matrix(0, nrow = Nsim, ncol = length(time2))

for(i in 1:Nsim){
	mu1[i,] = a1.est[i,] + b1.est[i,] * time1
	mu2[i,] = a2.est[i,] + b2.est[i,] * time2
}


obs1 = matrix(0, nrow = Nsim, ncol = length(time1))
obs2 = matrix(0, nrow = Nsim, ncol = length(time2))

set.seed(0704)
for(i in 1:Nsim){
	obs1[i,] = rnorm(length(time1), mu1[i,], sigma.y1[i])
	obs2[i,] = rnorm(length(time2), mu2[i,], sigma.y2[i])
}

#dim(obs1)
#[1]  200 8298
#dim(obs2)
#[1]   200 13032


# truncated predicted observations
# if obs <= 0, obs = 0; if obs >= 28, obs = 28
obs1.tr = obs1
obs2.tr = obs2
obs1.tr[obs1.tr <= 0] = 0
obs1.tr[obs1.tr >= 28] = 28
obs2.tr[obs2.tr <= 0] = 0
obs2.tr[obs2.tr >= 28] = 28


obs1.round = round(obs1)
obs2.round = round(obs2)



### posterior predictive checking
## Checking the goodness of fit of a model using Chi-square type statistic
chi_rep1 = chi_obs1 = NULL
for(i in 1:Nsim){
    chi_rep1[i] = sum( (obs1[i,] - mu1[i,])^2 / (sigma.y1[i]^2) )
    chi_obs1[i] = sum( (y1 - mu1[i,])^2 / (sigma.y1[i]^2) )
}

bpvalue1 = sum(chi_rep1 > chi_obs1) / Nsim


chi_rep2 = chi_obs2 = NULL
for(i in 1:Nsim){
    chi_rep2[i] = sum( (obs2[i,] - mu2[i,])^2 / (sigma.y2[i]^2) )
    chi_obs2[i] = sum( (y2 - mu2[i,])^2 / (sigma.y2[i]^2) )
}

bpvalue2 = sum(chi_rep2 > chi_obs2) / Nsim



## Checking structural assumptions of a model using Skewness and Kurtosis of standardized residuals
skew_rep1 = skew_obs1 = NULL
kurt_rep1 = kurt_obs1 = NULL
for(i in 1:Nsim){
    skew_rep1[i] = mean( ((obs1[i,] - mu1[i,]) / sigma.y1[i])^3 )
    skew_obs1[i] = mean( ((y1 - mu1[i,]) / sigma.y1[i])^3 )
    
    kurt_rep1[i] = mean( ((obs1[i,] - mu1[i,]) / sigma.y1[i])^4 ) - 3
    kurt_obs1[i] = mean( ((y1 - mu1[i,]) / sigma.y1[i])^4 ) - 3
}

sum(skew_rep1 > skew_obs1) / Nsim
sum(kurt_rep1 > kurt_obs1) / Nsim


## check the frequency of ADL scores 
apply(obs1, 1, FUN = function(x) sum(x <= 0))
apply(obs1, 1, FUN = function(x) sum(x >= 28))

apply(obs1.round, 1, FUN = function(x) sum(x == 27))


apply(obs2, 1, FUN = function(x) sum(x <= 0))
apply(obs2, 1, FUN = function(x) sum(x >= 28))

apply(obs2.round, 1, FUN = function(x) sum(x == 15))



par(mfrow = c(4,4))

for(i in 0:15) {
    hist(apply(obs1.round, 1, FUN = function(x) sum(x == i)), main = sum(y1 == i))
    abline(v = sum(y1 == i), col = "red")
}

for(i in 16:28) {
    hist(apply(obs1.round, 1, FUN = function(x) sum(x == i)), main = sum(y1 == i))
    abline(v = sum(y1 == i), col = "red")
}



