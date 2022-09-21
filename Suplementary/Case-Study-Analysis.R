##################################################################
# Fatoretto et al. A mean-dispersion bivariate framework to model counts and continuous responses - Podisus data analysis #
#Submitted to Biometrics - Journal of the International Biometric Society in October 2020#
##################################################################

#PACKAGES
require(optimx)
require(bbmle)
require(numDeriv)
require(Matrix)
require(tidyverse)
require(bivrp)
require(ggplot2)



dataset <- read.table("podisus_weight_eggs.csv", 
                    header = TRUE,
                    sep = ",")

#Scaterplot
variable_names <- list(
   "anti" = "Anticarsia gemmatalis" ,
   "diat" = "Diatraea saccharalis"
)

variable_labeller <- function(variable,value){
   return(variable_names[value])
}

ggplot(dataset,aes(x=weight_day_18,y=eggs_laid))  +
   geom_point(size=3,shape=21,fill='white') +theme_bw() + ylab("Number of eggs\n") + xlab("\nFemale weight")+
   facet_grid(~prey,labeller=variable_labeller)+
   theme(strip.text.x = element_text(size = 15),panel.grid.major = element_line(colour = 'white'),legend.position=c(0.2,0.4), legend.background = element_rect(linetype='solid',color=1), legend.title = element_text(size=10), legend.text = element_text(size=10), axis.title.x=element_text(size=15), axis.text.x=element_text(size=10), axis.title.y=element_text(size=15), axis.text.y=element_text(size=10),legend.key.height=unit(1.7,"line"),legend.key.width=unit(2,"line")) 




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


## Function to fit the model
mdbiv <- function(y,X1,X2,Xphi1,Xphi2,Xrho,link1,link2,data) {

   n <- nrow(data)
   y1 <- y[1:n]
   y2 <- y[(n+1):(2*n)]
   X <- as.matrix(bdiag(X1,X2))
   Xphi <- as.matrix(bdiag(Xphi1,Xphi2))
   n_rhos <- ncol(Xrho)

   link_fun1 <- generic_function(link1)$link_fun
   inv_link_fun1 <- generic_function(link1)$inv_link_fun
   deriv1 <- generic_function(link1)$deriv
   
   link_fun2 <- generic_function(link2)$link_fun
   inv_link_fun2<- generic_function(link2)$inv_link_fun
   deriv2 <- generic_function(link2)$deriv
   
##################################################################
##############STEP 1 - Fitting first outcome######################
##################################################################
   beta01<- solve(crossprod(X1,X1))%*%crossprod(X1,link_fun1(y1))
   eta01 <- X1%*%beta01
   mu.hat01 <- inv_link_fun1(eta01)
   Diag.mu01 <- diag(as.numeric(mu.hat01))
   res01 <- y1 - mu.hat01
   phi.init01 <- solve(crossprod(Xphi1,Xphi1))%*%t(Xphi1)%*%(diag(t(solve(Diag.mu01))%*%res01%*%t(res01)))
   Diag.phi01 <-  diag(as.numeric(Xphi1%*%phi.init01)) 
   all_parms02 <- list(beta01, phi.init01)

   criterion <- 0.2
   eps <- 1e-15
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
##############STEP 2 - Fitting second outcome#####################
##################################################################
   beta02 <- solve(crossprod(X2,X2))%*%crossprod(X2,link_fun2((y2)))
   eta02 <- X2%*%beta02
   mu.hat02 <- inv_link_fun2(eta02)
   Diag.mu02 <- diag(as.numeric(mu.hat02))
   res02 <- y2 - mu.hat02
   phi.init02 <- solve(crossprod(Xphi2,Xphi2))%*%t(Xphi2)%*%(diag(t(solve(Diag.mu02))%*%res02%*%t(res02)))
   Diag.phi02 <-  diag(as.numeric(Xphi2%*%phi.init02)) 
   all_parms12 <- list(beta02, phi.init02)


   criterion <- 0.2
   eps <- 1e-15
   iter <-2

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
##############STEP 3 - Fitting Bivariate model####################
##################################################################
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
      mu.hat <- rbind(inv_link_fun1(X1%*%beta[1:ncol(X1),]),inv_link_fun2(X2%*%beta[(ncol(X1)+1):(ncol(X)),]))
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
      eps <- 1e-15
      iter <- 2
      
      #Update values
      while(criterion > eps) {  
         all_parms <- all_parms2
         
         beta <- all_parms[[1]]
         mu.hat <- rbind(inv_link_fun1(X1%*%beta[1:ncol(X1),]),inv_link_fun2(X2%*%beta[(ncol(X1)+1):ncol(X),]))
         
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
         mu.hat <- rbind(inv_link_fun1(X1%*%beta[1:ncol(X1),]),inv_link_fun2(X2%*%beta[(ncol(X1)+1):ncol(X),]))
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
   
   coefs <- c(Beta_hat,Phi_hat,coef(summary(fitML))[,1])
   
   names(coefs) <- c(paste('beta1',0:(ncol(X1)-1),sep = "."),
                     paste('beta2',0:(ncol(X2)-1),sep = "."),
                     paste('phi1',0:(ncol(Xphi1)-1),sep = "."),
                     paste('phi2',0:(ncol(Xphi2)-1),sep = "."),
                     paste('rho',1:n_rhos,sep = "."))
   
   se.coefs <- c(stder_par1,coef(summary(fitML))[, 2])
   
   names(se.coefs) <- names(coefs)
      
   fitted <- c(mu.hat)
   resid <- y - fitted
   
   
      ret <- list("coefs" = coefs, "covariance" = Sigma3,"n" = n,"mu1"=mu.hat[1:n],"mu2" = mu.hat[(n+1):(2*n)],
                  "X" = X,"Xphi"=Xphi, "fitted" = fitted, "resid" = resid, "loglik" = logLik(fitML),"Y" = y,"se.coefs"=se.coefs)
      class(ret) <- "bivnormfit"
      return(ret)
}


#Defining the order of the response variables
y <- c(dataset$eggs_laid,dataset$weight_day_18)

#Defining model matrix
X1 <- model.matrix(~ 0+prey,dataset)
X2 <- model.matrix(~0+prey,dataset)
Xphi1 <-  model.matrix(~ 0+prey,dataset)
Xphi2 <-  model.matrix(~ 0+prey,dataset)
Xrho <- model.matrix(~0+prey,dataset)


#Fitting the model
fit0 <- mdbiv(y,X1,X2,Xphi1,Xphi2,Xrho,"log","identity",dataset)

#Extrating coefficients
round(fit0$coefs,4)
round(fit0$se.coefs,4)



###########################################################
## Function for extracting diagnostics (raw residuals)#########
##############################################################
dfun <- function(obj) {
   r <- obj$resid
   n <- obj$n
   return(list(r[1:n], r[(n+1):(2*n)]))
}


## Function for simulating new response variables
sfun <- function(obj) {
   n <- obj$n
   fitted <- obj$fitted
      y <- mvtnorm::rmvnorm(1,c(obj$mu1,obj$mu2),sigma=as.matrix(obj$covariance))
   y <- c(y)
   return(list(y[1:n], y[(n+1):(2*n)]))
}


ffun <- function(new.obj) {
   Ynew <- c(abs(new.obj[[1]]), new.obj[[2]])
   mdbiv(Ynew,X1,X2,Xphi1,Xphi2,Xrho,"log","identity",dataset)
}


## Bivariate residual plot for model fit0
   plot1 <- bivrp(fit0, diagfun=dfun, simfun=sfun, fitfun=ffun,verb=TRUE)
   


########################################################
############Likelihood ratio test
#######################################################

#testing the treatment effect on the mean of the first response variable -number of eggs laid  
 X1 <- model.matrix(~1,dataset)
 X2 <- model.matrix(~0+prey,dataset)
 fit1 <- mdbiv(y,X1,X2,Xphi1,Xphi2,Xrho,"log","identity",dataset)

 test <- 2*(fit0$loglik - fit1$loglik);test
 pchisq(test, 1, lower=FALSE) 
 
 
 #testing the treatment effect on the mean of the second response variable - weight
 X1 <- model.matrix(~0+prey,dataset)
 X2 <- model.matrix(~1,dataset)
 fit1 <- mdbiv(y,X1,X2,Xphi1,Xphi2,Xrho,"log","identity",dataset)
 
 test <- 2*(fit0$loglik - fit1$loglik);test
 pchisq(test, 1, lower=FALSE) 
   
   
#testing the treatment effect on the dispersion of the first response variable -number of eggs laid  
 X1 <- model.matrix(~0+prey,dataset)
 X2 <- model.matrix(~0+prey,dataset)
 
 Xphi1 <- model.matrix(~1,dataset)
 Xphi2 <- model.matrix(~0+prey,dataset)
 
 fit1 <- mdbiv(y,X1,X2,Xphi1,Xphi2,Xrho,"log","identity",dataset)

 test <- 2*(fit0$loglik - fit1$loglik);test
 pchisq(test, 2, lower=FALSE) 
 
 
 #testing the treatment effect on the dispersion of the second response variable - weight 
 Xphi1 <- model.matrix(~0+prey,dataset)
 Xphi2 <- model.matrix(~1,dataset)
 
 fit1 <- mdbiv(y,X1,X2,Xphi1,Xphi2,Xrho,"log","identity",dataset)
 
 test <- 2*(fit0$loglik - fit1$loglik);test
 pchisq(test, 2, lower=FALSE) 
 
 



