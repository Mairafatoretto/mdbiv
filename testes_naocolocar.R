rm(list = ls()) 

load("600_-08_1.RData")
Sim1 <- results[1:500]

load("600_-08_2.RData")
Sim2 <- results[1:500]

results <- c(Sim1,Sim2)
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


coverage <-rbind(round(sum(phi11_c)/1000,2),
                 round(sum(phi12_c)/1000,2),
                 round(sum(phi21_c)/1000,2),
                 round(sum(phi22_c)/1000,2),
                 round(sum(beta11_c)/1000,2),
                 round(sum(beta12_c)/1000,2),
                 round(sum(beta21_c)/1000,2),
                 round(sum(beta22_c)/1000,2),
                 round(sum(rho_c)/1000,2))



save.image(paste(600,0.2,sep="_",".RData"),compress="xz")

           


table
coverage






save.image("test.RData")


load("test.RData")
