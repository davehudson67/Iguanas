library(jagsUI)
library(rjags)
library(coda)

#load data

load("data.Rdata")
source("functions.R")

#Initial values
zinit <- cjs.init.z(jags.data$y,jags.data$f)
inits <- function(){ list(mu_Hs=29.3,sd_Hs=2.57,
                          mu_y2=runif(1,65,68),
                          maleEff_y2=runif(1,7,8),sd_y2=runif(1,1,2),
                          mu_k=runif(1,-7,-6.5),maleEff_k=runif(1,0.45,0.50),sd_k=runif(1,0.5,0.6),
                          sd_Lr=runif(1,0.1,0.5),psi=runif(1,0.45,0.55),
                          z=zinit,
                          betaS.sex=runif(1,-5,5),betaS.svl=runif(1,-5,5),
                          betaP.sex=runif(1,-5,5),betaP.svl=runif(1,-5,5),
                          sigma_phi=runif(1,0.25,0.30),sigma_p=runif(1,0.4,0.5)
)}


#Parameters
params<- c('mu_Hs','sd_Hs','mu_y2','sd_y2','maleEff_y2',
           'mu_k','maleEff_k','sd_k','sd_Lr','psi',
           'mean.p','mean.phi',
           'mu_phi','betaS.sex','betaS.svl',
           'mu_p','betaP.sex','betaP.svl',
           'sigma_phi','sigma_p'
           )

params_growth<- c('Hs', 'Lr_hat', 'y2', 'k','sex') #Individual Growth curve parameters



#####################################################################################################################################
####################################################################################################################################
set.seed(123)


n.adapt <- 20000
n.burnin <- 12000
n.iter <- 50000
n.thin <- 10  
n.chains <- 2

start.time = Sys.time()
out <- jags(data = jags.data,
            inits=inits,
            params,
            "Best_model_code.txt",
            n.chains=n.chains,
            n.adapt=n.adapt,
            n.thin=n.thin,
            n.iter=n.iter,
            n.burnin=n.burnin)
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'),2)
cat(paste(paste('Posterior computed in ', elapsed.time, sep=' '),'minutes\n', sep=' '))


library(MCMCvis)
MCMCsummary(out)
MCMCplot(out, params=c("betaS.sex","betaS.svl","betaP.sex","betaP.svl"))
MCMCtrace(out)

save(out, file = 'Plilfordi_MCMCsamples.RData')

############################################################################################################################################