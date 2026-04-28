	model {    
	#####################
	####Growth model####
	##################### 
  	for(i in 1:nind){
    	#SVL at hatching size
    	Hs[i] <- mu_Hs + eps_Hs[i]
    	eps_Hs[i] ~ dnorm(0, tau_Hs)
    	#maximum size
    	y2[i] <- mu_y2 + maleEff_y2 * sex[i] + eps_y2[i]
    	eps_y2[i] ~ dnorm(0, tau_y2)
    	#characteristic grwth rate
    	log(k[i]) <- mu_k + maleEff_k * sex[i] + eps_k[i]
    	eps_k[i] ~ dnorm(0, tau_k)
    
    	# first capture
    	Lr[i,f[i]] ~ dnorm(Lr_hat[i,f[i]], tau_Lr)         # SVL - observed
    	Lr_hat[i,f[i]] ~ dunif(30,80)   # expected size of animal at unknown age at first capture
    
    	sex[i] ~ dbern(psi)
  	}
  
  	#Priors Hsize
  	mu_Hs ~ dunif(25, 35)
  	tau_Hs <- 1 / sd_Hs^2
  	sd_Hs ~ dunif(2, 3) 
  
  	#Priors for y2
  	mu_y2 ~ dunif(65, 70)
  	maleEff_y2 ~ dunif(0, 15)
  	tau_y2 <- 1 / sd_y2^2
  	sd_y2 ~ dunif(0, 5)
  
  	#Priors for K
  	mu_k ~ dunif(-8,-5)
  	maleEff_k ~ dunif(0,1.5)
  	tau_k <- 1 / sd_k^2
  	sd_k ~ dunif(0,2)
  
  	psi ~ dunif(0, 1)
  	tau_Lr <- 1 / sd_Lr^2
  	sd_Lr ~ dunif(0, 5)

	for(i in 1:nind){
	  for(t in (f[i]+1):T){
	    #Schnute growth equation
	    Lr[i, t] ~ dnorm(Lr_hat[i, t], tau_Lr)   
	    Lr_hat[i, t] <- Lr_hat[i, t - 1] * exp(-k[i] * deltaT[t - 1]) + (y2[i] - Hs[i] * exp(-k[i] * (5110 - 0))) * (1 - exp(-k[i] * deltaT[t - 1])) / (1 - exp(-k[i] * (5110 - 0)))
	  }
	  Linf[i] <- (y2[i] - Hs[i] * exp(-k[i] * (5110 - 0))) / 1 - exp(-k[i] * (5110 - 0))
	  t0[i] <- 1 / k[i] * log((y2[i] - Hs[i]) / (y2[i] - Hs[i] * exp(-k[i] * (5110-0))))
	  
	  # standardize SVL
	  for(t in (f[i]):T){
	    SVL_st[i, t] <- (Lr_hat[i, t] - mean.svl) / sd.svl
	  }
	}
  	
  	#####################
  	####CJS model####
  	##################### 
  	
  	#------------------------------------------------------------------#
  	#                            Survival                     #
  	#------------------------------------------------------------------#
  	
  	#priors
  	mu_phi <- log(mean.phi / (1-mean.phi))
  	mean.phi ~ dunif(0, 1)
  	
  	betaS.sex~dunif(-5, 5)
  	betaS.svl~dunif(-5, 5)
  	
  	#Survival
  	for (i in 1:nind){
  	  for (t in f[i]:(T - 1)){
  	    logit(phi[i, t]) <- mu_phi +  betaS.sex * sex[i] + betaS.svl * SVL_st[i, t] + eps_phi[t] 
  	  }
  	}
  	
  	#------------------------------------------------------------------#
  	#                            Recapture                             #
  	#------------------------------------------------------------------#
  	
  	#priors
  	mu_p <- log(mean.p / (1 - mean.p))
  	mean.p ~ dunif(0, 1)
  	
  	betaP.sex~dunif(-5, 5)
  	betaP.svl~dunif(-5, 5)
  	
  	#Detectability
  	for (i in 1:nind){
  	  for (t in f[i]:(T - 1)){
  	    logit(p[i, t]) <-  mu_p + betaP.sex * sex[i] + betaP.svl * SVL_st[i, t] + eps_p[t]
  	  }
  	}
  	
  	for (t in 1:(T - 1)){
  	  eps_phi[t] ~ dnorm(0, tau_phi)
  	  eps_p[t] ~ dnorm(0, tau_p)
  	}
  	
  	tau_phi <- 1 / sigma_phi^2
  	sigma_phi ~ dunif(0, 5)
  	
  	tau_p <- 1 / sigma_p^2
  	sigma_p ~ dunif(0, 5)
  	
  	#likelihood
  	for (i in 1:nind){
  	  # Define latent state at first capture
  	  z[i, f[i]] <- 1
  	  for (t in (f[i] + 1):T){
  	    # State process
  	    z[i, t] ~ dbern(mu1[i, t])
  	    mu1[i, t] <- phi[i, t - 1] * z[i, t - 1]
  	    # Observation process
  	    y[i, t] ~ dbern(mu2[i, t])
  	    mu2[i, t] <- p[i, t - 1] * z[i, t]
  	  } 
  	}
	}
	
	
	
	
