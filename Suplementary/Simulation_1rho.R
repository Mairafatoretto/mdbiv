##################################################################
# Fatoretto et al. A bivariate framework to jointly model count and continuous responses - Simulation for 1 correlation parameters#
#Submitted to Biometrics - Journal of the International Biometric Society in October 2020#
##################################################################

#packages
require("parallel")
require(optimx)
require(pbmcapply)
require(numDeriv)
require(Matrix)
require(tidyverse)
require(bbmle)

#creating the scenarios
n <- 50 #100, 300 and 600
treat <- gl(2, n / 2)
plot <- gl(n / 5, 5)
dataset <- data.frame(treat,plot)

X1 <- model.matrix(~ 0 + treat)
X2 <- model.matrix(~ 0 + treat)
Xphi1 <-  model.matrix(~ 0 + treat)
Xphi2 <-  model.matrix(~ 0 + treat)

Xrho <- model.matrix(~1,treat)

link1 <- "log"
link2 <- "identity"
rho1 <- -0.8 #0.2, 0.5 and 0.8


#core processors using 
ncores <- 8



results <- pbmclapply(1:1000, function(j){   


   generic_function <- function(link){
      link_function <- switch (link,
                               "log" = log,
                               "sqrt" = sqrt,
                               "identity" = function(x) x,
                               "inverse" = function(x) 1/x
      )
      
      inverse_link_function <- switch(link,
                                      "log" = exp,
                                      "sqrt" = function(x) x^2,
                                      "inverse" = function(x) 1/x,
                                      "identity" = function(x) x,
                                      "logit" = function(x) plogis(x)
      )
      
      
      
      dfuncao <- switch(link,
                        "log"=expression({
                           .value <- log(x)
                           .grad <- array(0, c(length(.value), 1L), list(NULL,                           c("x")))
                           .grad[, "x"] <- 1/x
                           attr(.value, "gradient") <- .grad
                           .value
                        }),
                        "identity" = expression({
                           .value <- identity(x)
                           .grad <- array(0, c(length(.value), 1L), list(NULL,                         c("x")))
                           .grad[, "x"] <- 1
                           attr(.value, "gradient") <- .grad
                           .value
                        })
                        
      )         
      
      result <- list("link_fun"=link_function,
                        "inv_link_fun"=inverse_link_function,
                        "deriv" =dfuncao )
      return(result)
}
   

mdbiv <- function(X1,X2,Xphi1,Xphi2,Xrho,link1,link2,data){
      
      beta1 <- c(2, 3)
      phi1 <- c(.5, 2)
      beta2 <- c(10, 300)
      phi2 <- c(1, 10)

      n_rhos <- ncol(Xrho)
      

      link_fun1 <- generic_function(link1)$link_fun
      inv_link_fun1 <- generic_function(link1)$inv_link_fun
      deriv1 <- generic_function(link1)$deriv
      
      link_fun2 <- generic_function(link2)$link_fun
      inv_link_fun2 <- generic_function(link2)$inv_link_fun
      deriv2 <- generic_function(link2)$deriv
   
      

 
   
##################################################################
##############STEP 1 - First outcome######################
##################################################################

   
   #simulating
   mu1 <- inv_link_fun1(X1%*%beta1)
   Dmu1 <- diag(as.numeric(mu1))
   Dphi1 <- diag(as.numeric(Xphi1%*%phi1))
         
         
   y1 <- mvtnorm::rmvnorm(1,mu1,as.matrix(Dmu1%*%Dphi1))
         
         
   while (any(y1 < 0)) {
         y1 <- mvtnorm::rmvnorm(1, mu1, as.matrix(Dmu1 %*% Dphi1))
   }
         
   y1 <- c(y1)  
   
   
#fitting   
   beta01<- solve(crossprod(X1,X1))%*%crossprod(X1,link_fun1(y1))
   eta01 <- X1%*%beta01
   mu.hat01 <- inv_link_fun1(eta01)
   Diag.mu01 <- diag(as.numeric(mu.hat01))
   res01 <- y1 - mu.hat01
   phi.init01 <- solve(crossprod(Xphi1,Xphi1))%*%t(Xphi1)%*%((solve(Diag.mu01)%*%res01%*%t(res01)*diag(1,n))%*%as.matrix(rep(1,n)))
   Diag.phi01 <-  diag(as.numeric(Xphi1%*%phi.init01)) 
   all_parms02 <- list(beta01, phi.init01)
      
   criterion <- 0.2
   eps <- 1e-8 #epslon 
   iter <- 2
      

   while(criterion > eps) {
         all_parms01 <- all_parms02
         
         eta01 <- X1%*%all_parms01[[1]]
         mu.hat01 <- inv_link_fun1(eta01)
         Diag.mu01 <- diag(as.numeric(mu.hat01))
         Diag.phi01 <-  diag(as.numeric(Xphi1%*%all_parms01[[2]])) 
         
         
         x <- mu.hat01
         deriv_mu1 <- attr(eval(deriv1),"gradient")
         
         z01 <- eta01+(y1-mu.hat01)*deriv_mu1
         V <- Diag.phi01%*%Diag.mu01
         Diag_deriv_mu01 <- diag(as.numeric(deriv_mu1))
         w01 <-solve(V%*%Diag_deriv_mu01%*%Diag_deriv_mu01)
         
         beta01 <- solve(t(X1)%*%w01%*%X1)%*%t(X1)%*%w01%*%z01
         eta01 <- X1%*%beta01
         mu.hat01 <- inv_link_fun1(eta01)
         res01 <- y1 - mu.hat01
         phi.init01<- solve(crossprod(Xphi1,Xphi1))%*%t(Xphi1)%*%(diag(t(solve(Diag.mu01))%*%res01%*%t(res01)))
         
         all_parms02 <- list(beta01, phi.init01)
         criterion <- sum((unlist(all_parms02)-unlist(all_parms01))/unlist(all_parms01))^2
         iter <- iter + 1
         
   }
   
   betasML01 <-all_parms02[[1]]
   PhiML01 <- all_parms02[[2]]
   
   Sigma_1 <- diag(as.numeric(Xphi1%*%PhiML01))%*%diag(as.numeric(inv_link_fun1(X1%*%betasML01)))
   
##################################################################
##############STEP 2 - Cecond outcome#####################
##################################################################
 
#simulation  
   mu2 <- inv_link_fun2(X2%*%beta2)
   Dmu2 <- diag(as.numeric(mu2))
   Dphi2 <- diag(as.numeric(Xphi2%*%phi2))
   
   
   y2 <- mvtnorm::rmvnorm(1,mu2,as.matrix(Dmu2%*%Dphi2))
   
   y2 <- c(y2)
   
#fitting
   beta02 <- solve(crossprod(X2,X2))%*%crossprod(X2,link_fun2((y2)))
   eta02 <- X2%*%beta02
   mu.hat02 <- inv_link_fun2(eta02)
   Diag.mu02 <- diag(as.numeric(mu.hat02))
   res02 <- y2 - mu.hat02
   phi.init02 <- solve(crossprod(Xphi2,Xphi2))%*%t(Xphi2)%*%(diag(t(solve(Diag.mu02))%*%res02%*%t(res02)))
   Diag.phi02 <-  diag(as.numeric(Xphi2%*%phi.init02)) 
   all_parms12 <- list(beta02, phi.init02)
   
   
   criterion <- 0.2
   eps <- 1e-8
   iter <- 2
   

   while(criterion > eps) {
      all_parms11 <- all_parms12

      eta02 <- X2%*%all_parms11[[1]]
      mu.hat02 <- inv_link_fun2(eta02)
      Diag.mu02 <- diag(as.numeric(mu.hat02))
      Diag.phi02 <-  diag(as.numeric(Xphi2%*%all_parms11[[2]])) 
      
      x <- mu.hat02
      deriv_mu2 <- attr(eval(deriv2),"gradient")
      
      z02 <- eta02+(y2-mu.hat02)*deriv_mu2
      V2 <- Diag.phi02%*%Diag.mu02
      Diag_deriv_mu02 <- diag(as.numeric(deriv_mu2))
      w02 <-solve(V2%*%Diag_deriv_mu02%*%Diag_deriv_mu02)
      
      
      beta02 <- solve(t(X2)%*%w02%*%X2)%*%t(X2)%*%w02%*%z02
      eta02 <- X2%*%beta02
      mu.hat02 <- inv_link_fun2(eta02)
      Diag.mu02 <- diag(as.numeric(mu.hat02))
      res02 <- y2 - mu.hat02
      phi.init02<- solve(crossprod(Xphi2,Xphi2))%*%t(Xphi2)%*%(diag(t(solve(Diag.mu02))%*%res02%*%t(res02)))

      
      all_parms12 <- list(beta02, phi.init02)
      criterion <- sum((unlist(all_parms12)-unlist(all_parms11))/unlist(all_parms11))^2
      
      
      iter <- iter + 1
      
   }
   
   betasML02 <-all_parms12[[1]]
   PhiML02 <- all_parms12[[2]]
   
   Sigma_2 <- diag(as.numeric(Xphi2%*%PhiML02))%*%diag(as.numeric(inv_link_fun2(X2%*%betasML02)))
   
   
##################################################################
##############STEP 3 - Bivariate model####################
##################################################################

   #simulation
   a <-  diag(c(Xrho%*%as.numeric(c(rho1))))
   b <- diag(1,n)
   Sigma_b <- rbind(cbind(b,a),cbind(a,b))
   
   Sigma_init1 <- diag(as.numeric(Xphi1%*%phi1))%*%diag(as.numeric(inv_link_fun1(X1%*%beta1)))
   Sigma_init2 <- diag(as.numeric(Xphi2%*%phi2))%*%diag(as.numeric(inv_link_fun2(X2%*%beta2)))

   l1 <- list(t(chol(Sigma_init1)),t(chol(Sigma_init2)))
   l2 <- list(chol(Sigma_init1),chol(Sigma_init2))
   Sigma <- bdiag(l1)%*%Sigma_b%*%bdiag(l2)
   
   mu_init1 <- inv_link_fun1(X1%*%beta1)
   mu_init2 <- inv_link_fun2(X2%*%beta2)
   
   
   y <- mvtnorm::rmvnorm(1,c(mu_init1,mu_init2),sigma=as.matrix(Sigma))
   while (any(y[1:n] < 0)) {
   y <- mvtnorm::rmvnorm(1,c(mu_init1,mu_init2),sigma=as.matrix(Sigma))
   }
   
   y <- c(y)
   
   
   #fitting
   get_sigma_inv <- function(Sigma_1, Sigma_2, Sigma_b, Id,Diag_deriv_mu01,Diag_deriv_mu02) {
      chol_Sigma_1 <- chol(Sigma_1)
      chol_Sigma_2 <- chol(Sigma_2)
      l1 <- list(t(chol_Sigma_1),t(chol_Sigma_2))
      l2 <- list(chol_Sigma_1,chol_Sigma_2)
      Sigma <- bdiag(l1)%*%Sigma_b%*%bdiag(l2)
      Sigma.inv <- solve(Sigma)
      W0 <- Sigma.inv%*%bdiag(solve(Diag_deriv_mu01),solve(Diag_deriv_mu02))
      return(list("inverted_sigma" = Sigma.inv,
                  "sigma" = Sigma,
                  "W" = W0))
   }
   
   get_all_parms <- function(X, W0, n, y, Xphi, X1, X2) {
      beta <- as.matrix(solve(t(X)%*%W0%*%X)%*%t(X)%*%W0%*%c(link_fun1(y[1:n]),link_fun2(y[(n+1):(2*n)])))
      mu.hat <- rbind(inv_link_fun1(X1%*%beta[1:(length(beta)/2),]),inv_link_fun2(X2%*%beta[(length(beta)/2+1):(length(beta)),]))
      eta <- c(link_fun1(mu.hat[1:n]),link_fun2(mu.hat[(n+1):(2*n)]))
      Diag.mu <- diag(as.numeric(mu.hat))
      res <- y - mu.hat
      phi.init <- solve(crossprod(Xphi,Xphi))%*%t(Xphi)%*%(diag(t(solve(Diag.mu))%*%res%*%t(res)))
      all_parms2 <- list(beta, phi.init)
      
   }
   
   ##################################################################
   #Function to estimate correlation parameters by profile likelihood 
   #and update regression and dispersion parameters#
   ##################################################################  
   profloglik <- function(rho,y,X1,X2,Xphi1,Xphi2,Sigma_1,Sigma_2,Diag_deriv_mu01,Diag_deriv_mu02,Xrho,REML = FALSE) {
      
      Id <- (diag(length(y)/2))
      n <- dim(Id)[1]
      a <-  diag(c(Xrho%*%rho[1:n_rhos]))
      b <- diag(1,n)
      Sigma_b <- rbind(cbind(b,a),cbind(a,b))
      X <- as.matrix(bdiag(X1,X2))
      Xphi <- as.matrix(bdiag(Xphi1,Xphi2))
      
      # Get initial values
      W0 <- get_sigma_inv(Sigma_1 = Sigma_1,
                          Sigma_2 = Sigma_2,
                          Sigma_b = Sigma_b,
                          Id = Id,
                          Diag_deriv_mu01,
                          Diag_deriv_mu02)$W
      
      
      all_parms2 <- get_all_parms(X = X,
                                  W0 = W0,
                                  n = n,
                                  y = y,
                                  Xphi = Xphi,
                                  X1 = X1,
                                  X2 = X2)
      
      
      
      
      criterion <- 0.2
      eps <- 1e-8
      iter <- 2
      
      #Update values
      while(criterion > eps) {  
         all_parms <- all_parms2
         
         beta <- all_parms[[1]]
         mu.hat <- rbind(inv_link_fun1(X1%*%beta[1:(length(beta)/2),]),inv_link_fun2(X2%*%beta[(length(beta)/2+1):(length(beta)),]))
         
         Sigma01 <- diag(as.numeric(Xphi%*%all_parms[[2]]))%*%diag(as.numeric(mu.hat))
         
         x <- mu.hat[1:n]
         eta <- c(link_fun1(mu.hat[1:n]),link_fun2(mu.hat[(n+1):(2*n)]))
         deriv_mu1 <- attr(eval(deriv1),"gradient")
         z01 <- eta[1:n]+(y[1:n]-mu.hat[1:n])*deriv_mu1
         Diag_deriv_mu1 <- diag(as.numeric(deriv_mu1))
         
         
         x <- mu.hat[(n+1):(2*n)]
         deriv_mu2 <- attr(eval(deriv2),"gradient")
         z02 <- eta[(n+1):(2*n)]+(y[(n+1):(2*n)]-mu.hat[(n+1):(2*n)])*deriv_mu2
         Diag_deriv_mu2 <- diag(as.numeric(deriv_mu2))
         
         
         Sigma.inv2 <- get_sigma_inv(Sigma_1 = Sigma01[1:n,1:n],
                                     Sigma_2 = Sigma01[(n+1):(2*n),(n+1):(2*n)],
                                     Sigma_b = Sigma_b,
                                     Id = Id,
                                     Diag_deriv_mu01=Diag_deriv_mu1,
                                     Diag_deriv_mu02=Diag_deriv_mu2
         )$inverted_sigma
         
         
         W <- get_sigma_inv(Sigma_1 = Sigma01[1:n,1:n],
                            Sigma_2 = Sigma01[(n+1):(2*n),(n+1):(2*n)],
                            Sigma_b = Sigma_b,
                            Id = Id,
                            Diag_deriv_mu01=Diag_deriv_mu1,
                            Diag_deriv_mu02=Diag_deriv_mu2
         )$W
         
         
         
         beta <- as.matrix(solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%c(z01,z02))
         mu.hat <- rbind(inv_link_fun1(X1%*%beta[1:(length(beta)/2),]),inv_link_fun2(X2%*%beta[(length(beta)/2+1):(length(beta)),]))
         Diag.mu <- diag(as.numeric(mu.hat))
         res <- y - mu.hat
         phi.init <- solve(crossprod(Xphi,Xphi))%*%t(Xphi)%*%(diag(t(solve(Diag.mu))%*%res%*%t(res)))
         all_parms2 <- list(beta, phi.init)
         
         #Update criterion
         criterion <- sum((unlist(all_parms2)-unlist(all_parms))/unlist(all_parms))^2
         iter <- iter + 1
      }
      
      
      #Upadte Sigma
      Sigma1 <- diag(as.numeric(Xphi%*%phi.init))%*%diag(as.numeric(mu.hat))
      
      Sigma_3 <- get_sigma_inv(Sigma_1 = Sigma1[1:n,1:n],
                               Sigma_2 = Sigma1[(n+1):(2*n),(n+1):(2*n)],
                               Sigma_b = Sigma_b,
                               Id = Id,
                               Diag_deriv_mu01=Diag_deriv_mu1,
                               Diag_deriv_mu02=Diag_deriv_mu2)
      
      Sigma3 <- Sigma_3$sigma
      Sigma.inv3 <- Sigma_3$inverted_sigma
      llik <- -.5*(determinant(Sigma3)$modulus + t(res)%*%Sigma.inv3%*%res)
      a <- print(all_parms2)
      if(REML) llik <- llik -.5*determinant(t(X)%*%Sigma.inv3%*%X)$modulus
      return(-as.numeric(llik))
      
   }  

   # initial values
   rho <- runif(n_rhos,-1:1)
   names(rho) <- paste("rho",1:n_rhos,sep="")
   parnames(profloglik) <- names(rho) 


   # maximising profile loglikelihood 
   fitML <- mle2(profloglik, start  = rho,
                 method="L-BFGS-B",lower=c(rho=-0.999),
                 upper=c(rho=0.999),
                 data = list(y = y,X1=X1,X2=X2,Xphi1=Xphi1,Xphi2=Xphi2,Sigma_1=Sigma_1,Sigma_2=Sigma_2,Diag_deriv_mu01=Diag_deriv_mu01,Diag_deriv_mu02=Diag_deriv_mu02,Xrho=Xrho, REML = FALSE))

   #To extract correlation, regression and dispersion parameters
   param <- capture.output(profloglik(fitML@coef[1:n_rhos],y,X1,X2,Xphi1,Xphi2,Sigma_1,Sigma_2,Diag_deriv_mu01,Diag_deriv_mu02,Xrho))
   options(warn=-1)
   fvals <-na.omit(as.numeric(unlist(str_split(param, "\\s"))))
   
   X <- bdiag(X1,X2)
   Xphi <- bdiag(Xphi1,Xphi2)
   Beta_hat <- fvals[1:ncol(X)]
   Phi_hat <- fvals[(ncol(X)+1):(ncol(X)+ncol(Xphi))]
   
   
   
   #To extract the Covariance Matrix
   mu.hat <- rbind(inv_link_fun1(X1%*%Beta_hat[1:ncol(X1)]),inv_link_fun2(X2%*%Beta_hat[(ncol(X1)+1):(ncol(X))]))
   Sigma1 <-  (diag(as.numeric(Xphi%*%Phi_hat))%*%diag(as.numeric(mu.hat)))
   Id <- (diag(length(y)/2))
   n <- dim(Id)[1]
   a <-  diag(c(Xrho%*%fitML@coef))
   b <- diag(1,n)
   Sigma_b <- rbind(cbind(b,a),cbind(a,b))
   l11 <- list(t(chol(Sigma1[1:n,1:n])),t(chol(Sigma1[(n+1):(2*n),(n+1):(2*n)]))) 
   l12 <- list(chol(Sigma1[1:n,1:n]),chol(Sigma1[(n+1):(2*n),(n+1):(2*n)])) 
   Sigma3 <- bdiag(l11)%*%Sigma_b%*%bdiag(l12)
   
   
   #Numerically finding the Hessian matrix for the regression parameter using the estimated correlation parameter
   Hessian1 <- function(x){
      Id <- (diag(length(y)/2))
      n <- dim(Id)[1]
      mu.hat <- rbind(inv_link_fun1(X1%*%x[1:ncol(X1)]),inv_link_fun2(X2%*%x[(ncol(X1)+1):(ncol(X))]))
      res <- y - mu.hat
      Diag.mu <- diag(as.numeric(mu.hat))
      phi <- solve(crossprod(Xphi,Xphi))%*%t(Xphi)%*%(diag(t(solve(Diag.mu))%*%res%*%t(res)))
      Sigma1 <- diag(as.numeric(Xphi%*%phi))%*%diag(as.numeric(mu.hat))
      
      l11 <- list(t(chol(Sigma1[1:n,1:n])),t(chol(Sigma1[(n+1):(2*n),(n+1):(2*n)]))) 
      l12 <- list(chol(Sigma1[1:n,1:n]),chol(Sigma1[(n+1):(2*n),(n+1):(2*n)])) 
      
      a <- diag(c(Xrho%*%fitML@coef))
      b <- diag(1,n)
      Sigma_b <- rbind(cbind(b,a),cbind(a,b))
      Sigma3 <- bdiag(l11)%*%Sigma_b%*%bdiag(l12)
      Sigma.inv3 <- solve(Sigma3)
      return(-.5*(determinant(Sigma3)$modulus + t(res)%*%Sigma.inv3%*%res)[1,1])
   }
   
   #Standard errors for Betas
   x0 <- c(Beta_hat)
   variance_beta <- diag(solve(-hessian(func=Hessian1,
                                        x=x0)))
   stder_beta <- sqrt(variance_beta)
   
   
   #Numerically finding the Hessian matrix for the dispersion parameter using the estimated correlation and regression parameters
   Hessian2 <- function(phi,y,X1,X2,Xphi){
      Id <- (diag(length(y)/2))
      n <- dim(Id)[1]
      X <- as.matrix(bdiag(X1,X2))
      mu.hat <- rbind(inv_link_fun1(X1%*%Beta_hat[1:ncol(X1)]),inv_link_fun2(X2%*%Beta_hat[(ncol(X1)+1):(ncol(X))]))
      res <- y - mu.hat
      Diag.mu <- diag(as.numeric(mu.hat))
      Sigma1 <- diag(as.numeric(Xphi%*%phi))%*%diag(as.numeric(mu.hat))
      
      l11 <- list(t(chol(Sigma1[1:n,1:n])),t(chol(Sigma1[(n+1):(2*n),(n+1):(2*n)])))
      l12 <- list(chol(Sigma1[1:n,1:n]),chol(Sigma1[(n+1):(2*n),(n+1):(2*n)])) 
      
      a <-diag(c(Xrho%*%fitML@coef))
      b <- diag(1,n)
      Sigma_b <- rbind(cbind(b,a),cbind(a,b))
      Sigma3 <- bdiag(l11)%*%Sigma_b%*%bdiag(l12)
      Sigma.inv3 <- solve(Sigma3)
      llik <- -.5*(determinant(Sigma3)$modulus + t(res)%*%Sigma.inv3%*%res)
      return(-as.numeric(llik))
   }
   
   #Standard errors for phi
   phi <- runif(ncol(Xphi),0:1)
   op <- optim(c(phi),lower=0.0001,
               Hessian2, y=y, X1=X1,X2=X2,
               Xphi=Xphi,hessian = TRUE)
   variance_phi <-diag(abs(solve(-op$hessian)))
   stder_phi <- sqrt(variance_phi) 
   
   
   
   stder_par1 <- c(stder_beta,stder_phi)
   
   
   coefs <- c(Beta_hat,Phi_hat,fitML@coef)
   
   names(coefs) <- c(paste('beta1',0:(ncol(X1)-1),sep = "."),
                     paste('beta2',0:(ncol(X2)-1),sep = "."),
                     paste('phi1',0:(ncol(Xphi1)-1),sep = "."),
                     paste('phi2',0:(ncol(Xphi2)-1),sep = "."),
                     paste('rho',1:n_rhos,sep = "."))
   
   se.coefs <- c(stder_par1,coef(summary(fitML))[, 2])
   
   names(se.coefs) <- names(coefs)
   
   list(
      
      "betasML1" =Beta_hat[1:ncol(X1)],
      "error_betasML1" = stder_par1[1:ncol(X1)],
      
      "PhiML1" = Phi_hat[1:ncol(Xphi1)],
      "error_PhiML1" =stder_par1[(ncol(X)+1):(ncol(X)+ncol(Xphi1))],
      
      "betasML2"= Beta_hat[(ncol(X1)+1):ncol(X)],
      "error_betasML2" =stder_par1[(ncol(X1)+1):ncol(X)],
      
      "PhiML2" = Phi_hat[(ncol(Xphi1)+1):ncol(Xphi)],
      "error_PhiML2" = stder_par1[(ncol(X)+ncol(Xphi1)+1):(ncol(X)+ncol(Xphi))],
      
      "rho" = fitML@coef,
      "error_rho" = coef(summary(fitML))[, 2],
      
      "coefs" = coefs, "covariance" = Sigma3,"n" = n,"mu1"=mu.hat[1:n],"mu2" = mu.hat[(n+1):(2*n)],
      "X" = X,"Xphi"=Xphi, "fitted" = fitted, "resid" = resid, "loglik" = logLik(fitML),"Y" = y,"se.coefs"=se.coefs
      
   )
   
   }
   
   fit0 <- mdbiv(X1,X2,Xphi1,Xphi2,Xrho,link1,link2,dataset)
   

 return(fit0)

},mc.cores=ncores)


###############################################################
############Extracting values #################################
###############################################################
Sim <- rbind(results[1:1000])

coefs2 <- as.data.frame(do.call("rbind", Sim))


phi11 <- do.call(cbind, coefs2$PhiML1)[1,1:1000]
phi12 <- do.call(cbind, coefs2$PhiML1)[2,1:1000]
phi21 <- do.call(cbind, coefs2$PhiML2)[1,1:1000]
phi22 <- do.call(cbind, coefs2$PhiML2)[2,1:1000]

error_phi11 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML1))[,1]
error_phi12 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML1))[,2]
error_phi21 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML2))[,1]
error_phi22 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML2))[,2]


Beta11 <- do.call(cbind, coefs2$betasML1)[1,1:1000]
Beta12 <- do.call(cbind, coefs2$betasML1)[2,1:1000]
Beta21 <- do.call(cbind, coefs2$betasML2)[1,1:1000]
Beta22 <- do.call(cbind, coefs2$betasML2)[2,1:1000]

error_betasML11 <-  as.data.frame(do.call("rbind", coefs2$error_betasML1))[,1]
error_betasML12 <-  as.data.frame(do.call("rbind", coefs2$error_betasML1))[,2]
error_betasML21 <-  as.data.frame(do.call("rbind", coefs2$error_betasML2))[,1]
error_betasML22 <-  as.data.frame(do.call("rbind", coefs2$error_betasML2))[,2]



rho<-  as.data.frame(do.call("rbind", coefs2$rho))[,1]
error_rho <-  as.data.frame(do.call("rbind", coefs2$error_rho))[,1]


#getting results to Figure 2
table <- cbind(c(
   #bias
   mean(phi11)-0.5,
   mean(phi12)-2,
   mean(phi21)-1,
   mean(phi22)-10,
   mean(Beta11)-2,
   mean(Beta12)-3,
   mean(Beta21)-10,
   mean(Beta22)-300,
   mean(rho)-rho1),
   
   c(#lower limit
      mean(phi11)-0.5-1.96*mean(error_phi11),
      mean(phi12)-2-1.96*mean(error_phi12),
      mean(phi21)-1-1.96*mean(error_phi21),
      mean(phi22)-10-1.96*mean(error_phi22),
      mean(Beta11)-2-1.96*mean(error_betasML11),
      mean(Beta12)-3-1.96*mean(error_betasML12),
      mean(Beta21)-10-1.96*mean(error_betasML21),
      mean(Beta22)-300-1.96*mean(error_betasML22),
      mean(rho)-rho1-1.96*mean(error_rho)),
   
   c(#upper limit
      mean(phi11)-0.5+1.96*mean(error_phi11),
      mean(phi12)-2+1.96*mean(error_phi12),
      mean(phi21)-1+1.96*mean(error_phi21),
      mean(phi22)-10+1.96*mean(error_phi22),
      mean(Beta11)-2+1.96*mean(error_betasML11),
      mean(Beta12)-3+1.96*mean(error_betasML12),
      mean(Beta21)-10+1.96*mean(error_betasML21),
      mean(Beta22)-300+1.96*mean(error_betasML22),
      mean(rho)-rho1+1.96*mean(error_rho)))



############################################################
########Coverage rate#######################################
############################################################
int_phi11 <- cbind(lower = phi11-1.96*error_phi11,
                   upper = phi11+1.96*error_phi11)

int_phi12 <- cbind(lower = phi12-1.96*error_phi12,
                   upper = phi12+1.96*error_phi12)


int_phi21 <- cbind(lower =phi21-1.96*error_phi21,
                   upper =phi21+1.96*error_phi21)

int_phi22 <- cbind(lower = phi22-1.96*error_phi22,
                   upper = phi22+1.96*error_phi22)

int_Beta11 <- cbind(lower = Beta11-1.96*error_betasML11,
                    upper =  Beta11+1.96*error_betasML11)


int_Beta12 <- cbind(lower = Beta12-1.96*error_betasML12,
                    upper =  Beta12+1.96*error_betasML12)


int_Beta21 <- cbind(lower =Beta21-1.96*error_betasML21,
                    upper = Beta21+1.96*error_betasML21)

int_Beta22 <- cbind(lower = Beta22-1.96*error_betasML22,
                    upper =  Beta22+1.96*error_betasML22)

int_rho <- cbind(lower = rho -1.96*error_rho,
                 upper =  rho+1.96*error_rho)

phi11_c <- c(1:1000)
phi12_c <- c(1:1000)
phi21_c <- c(1:1000)
phi22_c <- c(1:1000)
beta11_c <- c(1:1000)
beta12_c <- c(1:1000)
beta21_c <- c(1:1000)
beta22_c <- c(1:1000)
rho_c <- c(1:1000)

for(i in 1:1000){
   
   if (int_phi11[i,1]<=0.5 & int_phi11[i,2]>=0.5)
   {
      phi11_c[i]<- 1
   }
   else{ 
      phi11_c[i] <-  0
   }
   
   
   if (int_phi12[i,1]<=2 & int_phi12[i,2]>=2)
   {
      phi12_c[i]<- 1
   }
   else{ 
      phi12_c[i] <-  0
   }
   
   if (int_phi21[i,1]<=1 & int_phi21[i,2]>=1)
   {
      phi21_c[i]<- 1
   }
   else{ 
      phi21_c[i] <-  0
   }
   
   if (int_phi22[i,1]<=10 & int_phi22[i,2]>=10)
   {
      phi22_c[i]<- 1
   }
   else{ 
      phi22_c[i] <-  0
   }
   
   
   if (int_Beta11[i,1]<=2 & int_Beta11[i,2]>=2)
   {
      beta11_c[i]<- 1
   }
   else{ 
      beta11_c[i] <-  0
   }
   
   if (int_Beta12[i,1]<=3 & int_Beta12[i,2]>=3)
   {
      beta12_c[i]<- 1
   }
   else{ 
      beta12_c[i] <-  0
   }
   
   if (int_Beta21[i,1]<=10 & int_Beta21[i,2]>=10)
   {
      beta21_c[i]<- 1
   }
   else{ 
      beta21_c[i] <-  0
   }
   
   
   if (int_Beta22[i,1]<=300 & int_Beta22[i,2]>=300)
   {
      beta22_c[i]<- 1
   }
   else{ 
      beta22_c[i] <-  0
   }
   
   
   if (int_rho[i,1]<=rho1 & int_rho[i,2]>=rho1)
   {
      rho_c[i]<- 1
   }
   else{ 
      rho_c[i] <-  0
   }
   
}

#getting results to Figure 3
coverage <-rbind(round(sum(phi11_c)/1000,2),
round(sum(phi12_c)/1000,2),
round(sum(phi21_c)/1000,2),
round(sum(phi22_c)/1000,2),
round(sum(beta11_c)/1000,2),
round(sum(beta12_c)/1000,2),
round(sum(beta21_c)/1000,2),
round(sum(beta22_c)/1000,2),
round(sum(rho_c)/1000,2))

#################################################################
#########saving results#########################################
#################################################################
save.image(paste(n,rho1,sep = "_",".RData"))

