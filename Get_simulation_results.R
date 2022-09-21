n <- 50
rho1 <- -0.8

load(paste(n,rho1,sep="_",".RData"))

Sim <- rbind(results[1:1000])

coefs2 <- as.data.frame(do.call("rbind", Sim))


phi11 <- do.call(cbind, coefs2$PhiML1)[1,1:1000]
phi12 <- do.call(cbind, coefs2$PhiML1)[2,1:1000]
phi21 <- do.call(cbind, coefs2$PhiML2)[1,1:1000]
phi22 <- do.call(cbind, coefs2$PhiML2)[2,1:1000]

error_phi1 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML1))
error_phi2 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML2))

Beta11 <- do.call(cbind, coefs2$betasML1)[1,1:1000]
Beta12 <- do.call(cbind, coefs2$betasML1)[2,1:1000]
Beta21 <- do.call(cbind, coefs2$betasML2)[1,1:1000]
Beta22 <- do.call(cbind, coefs2$betasML2)[2,1:1000]

error_betasML2 <-  as.data.frame(do.call("rbind", coefs2$error_betasML2))
error_betasML1 <-  as.data.frame(do.call("rbind", coefs2$error_betasML1))


rho<-  as.data.frame(do.call("rbind", coefs2$rho))
error_rho <-  as.data.frame(do.call("rbind", coefs2$error_rho))



table1 <- cbind(c(
#bias
mean(phi11)-0.5,
mean(phi12)-2,
mean(phi21)-1,
mean(phi22)-10,
mean(Beta11)-2,
mean(Beta12)-3,
mean(Beta21)-10,
mean(Beta22)-300,
mean(rho$V1)-rho1),

c(#lower limit
mean(phi11)-0.5-1.96*mean(error_phi1$phi1),
mean(phi12)-2-1.96*mean(error_phi1$phi2),
mean(phi21)-1-1.96*mean(error_phi2$phi3),
mean(phi22)-10-1.96*mean(error_phi2$phi4),
mean(Beta11)-2-1.96*mean(error_betasML1$V1),
mean(Beta12)-3-1.96*mean(error_betasML1$V2),
mean(Beta21)-10-1.96*mean(error_betasML2$V1),
mean(Beta22)-300-1.96*mean(error_betasML2$V2),
mean(rho$V1)-rho1-1.96*mean(error_rho$V1)),

c(#upper limit
mean(phi11)-0.5+1.96*mean(error_phi1$phi1),
mean(phi12)-2+1.96*mean(error_phi1$phi2),
mean(phi21)-1+1.96*mean(error_phi2$phi3),
mean(phi22)-10+1.96*mean(error_phi2$phi4),
mean(Beta11)-2+1.96*mean(error_betasML1$V1),
mean(Beta12)-3+1.96*mean(error_betasML1$V2),
mean(Beta21)-10+1.96*mean(error_betasML2$V1),
mean(Beta22)-300+1.96*mean(error_betasML2$V2),
mean(rho$V1)-rho1+1.96*mean(error_rho$V1)))


n <- 50
rho1 <- 0.2

load(paste(n,rho1,sep="_",".RData"))

Sim <- rbind(results[1:1000])

coefs2 <- as.data.frame(do.call("rbind", Sim))


phi11 <- do.call(cbind, coefs2$PhiML1)[1,1:1000]
phi12 <- do.call(cbind, coefs2$PhiML1)[2,1:1000]
phi21 <- do.call(cbind, coefs2$PhiML2)[1,1:1000]
phi22 <- do.call(cbind, coefs2$PhiML2)[2,1:1000]

error_phi1 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML1))
error_phi2 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML2))

Beta11 <- do.call(cbind, coefs2$betasML1)[1,1:1000]
Beta12 <- do.call(cbind, coefs2$betasML1)[2,1:1000]
Beta21 <- do.call(cbind, coefs2$betasML2)[1,1:1000]
Beta22 <- do.call(cbind, coefs2$betasML2)[2,1:1000]

error_betasML2 <-  as.data.frame(do.call("rbind", coefs2$error_betasML2))
error_betasML1 <-  as.data.frame(do.call("rbind", coefs2$error_betasML1))


rho<-  as.data.frame(do.call("rbind", coefs2$rho))
error_rho <-  as.data.frame(do.call("rbind", coefs2$error_rho))



table2 <- cbind(c(
   #bias
   mean(phi11)-0.5,
   mean(phi12)-2,
   mean(phi21)-1,
   mean(phi22)-10,
   mean(Beta11)-2,
   mean(Beta12)-3,
   mean(Beta21)-10,
   mean(Beta22)-300,
   mean(rho$V1)-rho1),
   
   c(#lower limit
      mean(phi11)-0.5-1.96*mean(error_phi1$phi1),
      mean(phi12)-2-1.96*mean(error_phi1$phi2),
      mean(phi21)-1-1.96*mean(error_phi2$phi3),
      mean(phi22)-10-1.96*mean(error_phi2$phi4),
      mean(Beta11)-2-1.96*mean(error_betasML1$V1),
      mean(Beta12)-3-1.96*mean(error_betasML1$V2),
      mean(Beta21)-10-1.96*mean(error_betasML2$V1),
      mean(Beta22)-300-1.96*mean(error_betasML2$V2),
      mean(rho$V1)-rho1-1.96*mean(error_rho$V1)),
   
   c(#upper limit
      mean(phi11)-0.5+1.96*mean(error_phi1$phi1),
      mean(phi12)-2+1.96*mean(error_phi1$phi2),
      mean(phi21)-1+1.96*mean(error_phi2$phi3),
      mean(phi22)-10+1.96*mean(error_phi2$phi4),
      mean(Beta11)-2+1.96*mean(error_betasML1$V1),
      mean(Beta12)-3+1.96*mean(error_betasML1$V2),
      mean(Beta21)-10+1.96*mean(error_betasML2$V1),
      mean(Beta22)-300+1.96*mean(error_betasML2$V2),
      mean(rho$V1)-rho1+1.96*mean(error_rho$V1)))


n <- 50
rho1 <- 0.5

load(paste(n,rho1,sep="_",".RData"))

Sim <- rbind(results[1:1000])

coefs2 <- as.data.frame(do.call("rbind", Sim))


phi11 <- do.call(cbind, coefs2$PhiML1)[1,1:1000]
phi12 <- do.call(cbind, coefs2$PhiML1)[2,1:1000]
phi21 <- do.call(cbind, coefs2$PhiML2)[1,1:1000]
phi22 <- do.call(cbind, coefs2$PhiML2)[2,1:1000]

error_phi1 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML1))
error_phi2 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML2))

Beta11 <- do.call(cbind, coefs2$betasML1)[1,1:1000]
Beta12 <- do.call(cbind, coefs2$betasML1)[2,1:1000]
Beta21 <- do.call(cbind, coefs2$betasML2)[1,1:1000]
Beta22 <- do.call(cbind, coefs2$betasML2)[2,1:1000]

error_betasML2 <-  as.data.frame(do.call("rbind", coefs2$error_betasML2))
error_betasML1 <-  as.data.frame(do.call("rbind", coefs2$error_betasML1))


rho<-  as.data.frame(do.call("rbind", coefs2$rho))
error_rho <-  as.data.frame(do.call("rbind", coefs2$error_rho))



table3 <- cbind(c(
   #bias
   mean(phi11)-0.5,
   mean(phi12)-2,
   mean(phi21)-1,
   mean(phi22)-10,
   mean(Beta11)-2,
   mean(Beta12)-3,
   mean(Beta21)-10,
   mean(Beta22)-300,
   mean(rho$V1)-rho1),
   
   c(#lower limit
      mean(phi11)-0.5-1.96*mean(error_phi1$phi1),
      mean(phi12)-2-1.96*mean(error_phi1$phi2),
      mean(phi21)-1-1.96*mean(error_phi2$phi3),
      mean(phi22)-10-1.96*mean(error_phi2$phi4),
      mean(Beta11)-2-1.96*mean(error_betasML1$V1),
      mean(Beta12)-3-1.96*mean(error_betasML1$V2),
      mean(Beta21)-10-1.96*mean(error_betasML2$V1),
      mean(Beta22)-300-1.96*mean(error_betasML2$V2),
      mean(rho$V1)-rho1-1.96*mean(error_rho$V1)),
   
   c(#upper limit
      mean(phi11)-0.5+1.96*mean(error_phi1$phi1),
      mean(phi12)-2+1.96*mean(error_phi1$phi2),
      mean(phi21)-1+1.96*mean(error_phi2$phi3),
      mean(phi22)-10+1.96*mean(error_phi2$phi4),
      mean(Beta11)-2+1.96*mean(error_betasML1$V1),
      mean(Beta12)-3+1.96*mean(error_betasML1$V2),
      mean(Beta21)-10+1.96*mean(error_betasML2$V1),
      mean(Beta22)-300+1.96*mean(error_betasML2$V2),
      mean(rho$V1)-rho1+1.96*mean(error_rho$V1)))


n <- 50
rho1 <- 0.8

load(paste(n,rho1,sep="_",".RData"))

Sim <- rbind(results[1:1000])

coefs2 <- as.data.frame(do.call("rbind", Sim))


phi11 <- do.call(cbind, coefs2$PhiML1)[1,1:1000]
phi12 <- do.call(cbind, coefs2$PhiML1)[2,1:1000]
phi21 <- do.call(cbind, coefs2$PhiML2)[1,1:1000]
phi22 <- do.call(cbind, coefs2$PhiML2)[2,1:1000]

error_phi1 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML1))
error_phi2 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML2))

Beta11 <- do.call(cbind, coefs2$betasML1)[1,1:1000]
Beta12 <- do.call(cbind, coefs2$betasML1)[2,1:1000]
Beta21 <- do.call(cbind, coefs2$betasML2)[1,1:1000]
Beta22 <- do.call(cbind, coefs2$betasML2)[2,1:1000]

error_betasML2 <-  as.data.frame(do.call("rbind", coefs2$error_betasML2))
error_betasML1 <-  as.data.frame(do.call("rbind", coefs2$error_betasML1))


rho<-  as.data.frame(do.call("rbind", coefs2$rho))
error_rho <-  as.data.frame(do.call("rbind", coefs2$error_rho))



table4 <- cbind(c(
   #bias
   mean(phi11)-0.5,
   mean(phi12)-2,
   mean(phi21)-1,
   mean(phi22)-10,
   mean(Beta11)-2,
   mean(Beta12)-3,
   mean(Beta21)-10,
   mean(Beta22)-300,
   mean(rho$V1)-rho1),
   
   c(#lower limit
      mean(phi11)-0.5-1.96*mean(error_phi1$phi1),
      mean(phi12)-2-1.96*mean(error_phi1$phi2),
      mean(phi21)-1-1.96*mean(error_phi2$phi3),
      mean(phi22)-10-1.96*mean(error_phi2$phi4),
      mean(Beta11)-2-1.96*mean(error_betasML1$V1),
      mean(Beta12)-3-1.96*mean(error_betasML1$V2),
      mean(Beta21)-10-1.96*mean(error_betasML2$V1),
      mean(Beta22)-300-1.96*mean(error_betasML2$V2),
      mean(rho$V1)-rho1-1.96*mean(error_rho$V1)),
   
   c(#upper limit
      mean(phi11)-0.5+1.96*mean(error_phi1$phi1),
      mean(phi12)-2+1.96*mean(error_phi1$phi2),
      mean(phi21)-1+1.96*mean(error_phi2$phi3),
      mean(phi22)-10+1.96*mean(error_phi2$phi4),
      mean(Beta11)-2+1.96*mean(error_betasML1$V1),
      mean(Beta12)-3+1.96*mean(error_betasML1$V2),
      mean(Beta21)-10+1.96*mean(error_betasML2$V1),
      mean(Beta22)-300+1.96*mean(error_betasML2$V2),
      mean(rho$V1)-rho1+1.96*mean(error_rho$V1)))




n <- 100
rho1 <- -0.8

load(paste(n,rho1,sep="_",".RData"))

Sim <- rbind(results[1:1000])

coefs2 <- as.data.frame(do.call("rbind", Sim))


phi11 <- do.call(cbind, coefs2$PhiML1)[1,1:1000]
phi12 <- do.call(cbind, coefs2$PhiML1)[2,1:1000]
phi21 <- do.call(cbind, coefs2$PhiML2)[1,1:1000]
phi22 <- do.call(cbind, coefs2$PhiML2)[2,1:1000]

error_phi1 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML1))
error_phi2 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML2))

Beta11 <- do.call(cbind, coefs2$betasML1)[1,1:1000]
Beta12 <- do.call(cbind, coefs2$betasML1)[2,1:1000]
Beta21 <- do.call(cbind, coefs2$betasML2)[1,1:1000]
Beta22 <- do.call(cbind, coefs2$betasML2)[2,1:1000]

error_betasML2 <-  as.data.frame(do.call("rbind", coefs2$error_betasML2))
error_betasML1 <-  as.data.frame(do.call("rbind", coefs2$error_betasML1))


rho<-  as.data.frame(do.call("rbind", coefs2$rho))
error_rho <-  as.data.frame(do.call("rbind", coefs2$error_rho))



table5 <- cbind(c(
   #bias
   mean(phi11)-0.5,
   mean(phi12)-2,
   mean(phi21)-1,
   mean(phi22)-10,
   mean(Beta11)-2,
   mean(Beta12)-3,
   mean(Beta21)-10,
   mean(Beta22)-300,
   mean(rho$V1)-rho1),
   
   c(#lower limit
      mean(phi11)-0.5-1.96*mean(error_phi1$V1),
      mean(phi12)-2-1.96*mean(error_phi1$V2),
      mean(phi21)-1-1.96*mean(error_phi2$V1),
      mean(phi22)-10-1.96*mean(error_phi2$V2),
      mean(Beta11)-2-1.96*mean(error_betasML1$V1),
      mean(Beta12)-3-1.96*mean(error_betasML1$V2),
      mean(Beta21)-10-1.96*mean(error_betasML2$V1),
      mean(Beta22)-300-1.96*mean(error_betasML2$V2),
      mean(rho$V1)-rho1-1.96*mean(error_rho$V1)),
   
   c(#upper limit
      mean(phi11)-0.5+1.96*mean(error_phi1$V1),
      mean(phi12)-2+1.96*mean(error_phi1$V2),
      mean(phi21)-1+1.96*mean(error_phi2$V1),
      mean(phi22)-10+1.96*mean(error_phi2$V2),
      mean(Beta11)-2+1.96*mean(error_betasML1$V1),
      mean(Beta12)-3+1.96*mean(error_betasML1$V2),
      mean(Beta21)-10+1.96*mean(error_betasML2$V1),
      mean(Beta22)-300+1.96*mean(error_betasML2$V2),
      mean(rho$V1)-rho1+1.96*mean(error_rho$V1)))


n <- 100
rho1 <- 0.2

load(paste(n,rho1,sep="_",".RData"))

Sim <- rbind(results[1:1000])

coefs2 <- as.data.frame(do.call("rbind", Sim))


phi11 <- do.call(cbind, coefs2$PhiML1)[1,1:1000]
phi12 <- do.call(cbind, coefs2$PhiML1)[2,1:1000]
phi21 <- do.call(cbind, coefs2$PhiML2)[1,1:1000]
phi22 <- do.call(cbind, coefs2$PhiML2)[2,1:1000]

error_phi1 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML1))
error_phi2 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML2))

Beta11 <- do.call(cbind, coefs2$betasML1)[1,1:1000]
Beta12 <- do.call(cbind, coefs2$betasML1)[2,1:1000]
Beta21 <- do.call(cbind, coefs2$betasML2)[1,1:1000]
Beta22 <- do.call(cbind, coefs2$betasML2)[2,1:1000]

error_betasML2 <-  as.data.frame(do.call("rbind", coefs2$error_betasML2))
error_betasML1 <-  as.data.frame(do.call("rbind", coefs2$error_betasML1))


rho<-  as.data.frame(do.call("rbind", coefs2$rho))
error_rho <-  as.data.frame(do.call("rbind", coefs2$error_rho))



table6 <- cbind(c(
   #bias
   mean(phi11)-0.5,
   mean(phi12)-2,
   mean(phi21)-1,
   mean(phi22)-10,
   mean(Beta11)-2,
   mean(Beta12)-3,
   mean(Beta21)-10,
   mean(Beta22)-300,
   mean(rho$V1)-rho1),
   
   c(#lower limit
      mean(phi11)-0.5-1.96*mean(error_phi1$phi1),
      mean(phi12)-2-1.96*mean(error_phi1$phi2),
      mean(phi21)-1-1.96*mean(error_phi2$phi3),
      mean(phi22)-10-1.96*mean(error_phi2$phi4),
      mean(Beta11)-2-1.96*mean(error_betasML1$V1),
      mean(Beta12)-3-1.96*mean(error_betasML1$V2),
      mean(Beta21)-10-1.96*mean(error_betasML2$V1),
      mean(Beta22)-300-1.96*mean(error_betasML2$V2),
      mean(rho$V1)-rho1-1.96*mean(error_rho$V1)),
   
   c(#upper limit
      mean(phi11)-0.5+1.96*mean(error_phi1$phi1),
      mean(phi12)-2+1.96*mean(error_phi1$phi2),
      mean(phi21)-1+1.96*mean(error_phi2$phi3),
      mean(phi22)-10+1.96*mean(error_phi2$phi4),
      mean(Beta11)-2+1.96*mean(error_betasML1$V1),
      mean(Beta12)-3+1.96*mean(error_betasML1$V2),
      mean(Beta21)-10+1.96*mean(error_betasML2$V1),
      mean(Beta22)-300+1.96*mean(error_betasML2$V2),
      mean(rho$V1)-rho1+1.96*mean(error_rho$V1)))


n <- 100
rho1 <- 0.5

load(paste(n,rho1,sep="_",".RData"))

Sim <- rbind(results[1:1000])

coefs2 <- as.data.frame(do.call("rbind", Sim))


phi11 <- do.call(cbind, coefs2$PhiML1)[1,1:1000]
phi12 <- do.call(cbind, coefs2$PhiML1)[2,1:1000]
phi21 <- do.call(cbind, coefs2$PhiML2)[1,1:1000]
phi22 <- do.call(cbind, coefs2$PhiML2)[2,1:1000]

error_phi1 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML1))
error_phi2 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML2))

Beta11 <- do.call(cbind, coefs2$betasML1)[1,1:1000]
Beta12 <- do.call(cbind, coefs2$betasML1)[2,1:1000]
Beta21 <- do.call(cbind, coefs2$betasML2)[1,1:1000]
Beta22 <- do.call(cbind, coefs2$betasML2)[2,1:1000]

error_betasML2 <-  as.data.frame(do.call("rbind", coefs2$error_betasML2))
error_betasML1 <-  as.data.frame(do.call("rbind", coefs2$error_betasML1))


rho<-  as.data.frame(do.call("rbind", coefs2$rho))
error_rho <-  as.data.frame(do.call("rbind", coefs2$error_rho))



table7 <- cbind(c(
   #bias
   mean(phi11)-0.5,
   mean(phi12)-2,
   mean(phi21)-1,
   mean(phi22)-10,
   mean(Beta11)-2,
   mean(Beta12)-3,
   mean(Beta21)-10,
   mean(Beta22)-300,
   mean(rho$V1)-rho1),
   
   c(#lower limit
      mean(phi11)-0.5-1.96*mean(error_phi1$phi1),
      mean(phi12)-2-1.96*mean(error_phi1$phi2),
      mean(phi21)-1-1.96*mean(error_phi2$phi3),
      mean(phi22)-10-1.96*mean(error_phi2$phi4),
      mean(Beta11)-2-1.96*mean(error_betasML1$V1),
      mean(Beta12)-3-1.96*mean(error_betasML1$V2),
      mean(Beta21)-10-1.96*mean(error_betasML2$V1),
      mean(Beta22)-300-1.96*mean(error_betasML2$V2),
      mean(rho$V1)-rho1-1.96*mean(error_rho$V1)),
   
   c(#upper limit
      mean(phi11)-0.5+1.96*mean(error_phi1$phi1),
      mean(phi12)-2+1.96*mean(error_phi1$phi2),
      mean(phi21)-1+1.96*mean(error_phi2$phi3),
      mean(phi22)-10+1.96*mean(error_phi2$phi4),
      mean(Beta11)-2+1.96*mean(error_betasML1$V1),
      mean(Beta12)-3+1.96*mean(error_betasML1$V2),
      mean(Beta21)-10+1.96*mean(error_betasML2$V1),
      mean(Beta22)-300+1.96*mean(error_betasML2$V2),
      mean(rho$V1)-rho1+1.96*mean(error_rho$V1)))


n <- 100
rho1 <- 0.8

load(paste(n,rho1,sep="_",".RData"))

Sim <- rbind(results[1:1000])

coefs2 <- as.data.frame(do.call("rbind", Sim))


phi11 <- do.call(cbind, coefs2$PhiML1)[1,1:1000]
phi12 <- do.call(cbind, coefs2$PhiML1)[2,1:1000]
phi21 <- do.call(cbind, coefs2$PhiML2)[1,1:1000]
phi22 <- do.call(cbind, coefs2$PhiML2)[2,1:1000]

error_phi1 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML1))
error_phi2 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML2))

Beta11 <- do.call(cbind, coefs2$betasML1)[1,1:1000]
Beta12 <- do.call(cbind, coefs2$betasML1)[2,1:1000]
Beta21 <- do.call(cbind, coefs2$betasML2)[1,1:1000]
Beta22 <- do.call(cbind, coefs2$betasML2)[2,1:1000]

error_betasML2 <-  as.data.frame(do.call("rbind", coefs2$error_betasML2))
error_betasML1 <-  as.data.frame(do.call("rbind", coefs2$error_betasML1))


rho<-  as.data.frame(do.call("rbind", coefs2$rho))
error_rho <-  as.data.frame(do.call("rbind", coefs2$error_rho))



table8 <- cbind(c(
   #bias
   mean(phi11)-0.5,
   mean(phi12)-2,
   mean(phi21)-1,
   mean(phi22)-10,
   mean(Beta11)-2,
   mean(Beta12)-3,
   mean(Beta21)-10,
   mean(Beta22)-300,
   mean(rho$V1)-rho1),
   
   c(#lower limit
      mean(phi11)-0.5-1.96*mean(error_phi1$phi1),
      mean(phi12)-2-1.96*mean(error_phi1$phi2),
      mean(phi21)-1-1.96*mean(error_phi2$phi3),
      mean(phi22)-10-1.96*mean(error_phi2$phi4),
      mean(Beta11)-2-1.96*mean(error_betasML1$V1),
      mean(Beta12)-3-1.96*mean(error_betasML1$V2),
      mean(Beta21)-10-1.96*mean(error_betasML2$V1),
      mean(Beta22)-300-1.96*mean(error_betasML2$V2),
      mean(rho$V1)-rho1-1.96*mean(error_rho$V1)),
   
   c(#upper limit
      mean(phi11)-0.5+1.96*mean(error_phi1$phi1),
      mean(phi12)-2+1.96*mean(error_phi1$phi2),
      mean(phi21)-1+1.96*mean(error_phi2$phi3),
      mean(phi22)-10+1.96*mean(error_phi2$phi4),
      mean(Beta11)-2+1.96*mean(error_betasML1$V1),
      mean(Beta12)-3+1.96*mean(error_betasML1$V2),
      mean(Beta21)-10+1.96*mean(error_betasML2$V1),
      mean(Beta22)-300+1.96*mean(error_betasML2$V2),
      mean(rho$V1)-rho1+1.96*mean(error_rho$V1)))


n <- 300
rho1 <- -0.8

load(paste(n,rho1,sep="_",".RData"))

Sim <- rbind(results[1:1000])

coefs2 <- as.data.frame(do.call("rbind", Sim))


phi11 <- do.call(cbind, coefs2$PhiML1)[1,1:1000]
phi12 <- do.call(cbind, coefs2$PhiML1)[2,1:1000]
phi21 <- do.call(cbind, coefs2$PhiML2)[1,1:1000]
phi22 <- do.call(cbind, coefs2$PhiML2)[2,1:1000]

error_phi1 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML1))
error_phi2 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML2))

Beta11 <- do.call(cbind, coefs2$betasML1)[1,1:1000]
Beta12 <- do.call(cbind, coefs2$betasML1)[2,1:1000]
Beta21 <- do.call(cbind, coefs2$betasML2)[1,1:1000]
Beta22 <- do.call(cbind, coefs2$betasML2)[2,1:1000]

error_betasML2 <-  as.data.frame(do.call("rbind", coefs2$error_betasML2))
error_betasML1 <-  as.data.frame(do.call("rbind", coefs2$error_betasML1))


rho<-  as.data.frame(do.call("rbind", coefs2$rho))
error_rho <-  as.data.frame(do.call("rbind", coefs2$error_rho))



table9 <- cbind(c(
   #bias
   mean(phi11)-0.5,
   mean(phi12)-2,
   mean(phi21)-1,
   mean(phi22)-10,
   mean(Beta11)-2,
   mean(Beta12)-3,
   mean(Beta21)-10,
   mean(Beta22)-300,
   mean(rho$rho.rho1)-rho1),
   
   c(#lower limit
      mean(phi11)-0.5-1.96*mean(error_phi1$V1),
      mean(phi12)-2-1.96*mean(error_phi1$V2),
      mean(phi21)-1-1.96*mean(error_phi2$V1),
      mean(phi22)-10-1.96*mean(error_phi2$V2),
      mean(Beta11)-2-1.96*mean(error_betasML1$V1),
      mean(Beta12)-3-1.96*mean(error_betasML1$V2),
      mean(Beta21)-10-1.96*mean(error_betasML2$V1),
      mean(Beta22)-300-1.96*mean(error_betasML2$V2),
      mean(rho$rho.rho1)-rho1-1.96*mean(error_rho$rho.rho1)),
   
   c(#upper limit
      mean(phi11)-0.5+1.96*mean(error_phi1$V1),
      mean(phi12)-2+1.96*mean(error_phi1$V2),
      mean(phi21)-1+1.96*mean(error_phi2$V1),
      mean(phi22)-10+1.96*mean(error_phi2$V2),
      mean(Beta11)-2+1.96*mean(error_betasML1$V1),
      mean(Beta12)-3+1.96*mean(error_betasML1$V2),
      mean(Beta21)-10+1.96*mean(error_betasML2$V1),
      mean(Beta22)-300+1.96*mean(error_betasML2$V2),
      mean(rho$rho.rho1)-rho1+1.96*mean(error_rho$rho.rho1)))


n <- 300
rho1 <- 0.2

load(paste(n,rho1,sep="_",".RData"))

Sim <- rbind(results[1:1000])

coefs2 <- as.data.frame(do.call("rbind", Sim))


phi11 <- do.call(cbind, coefs2$PhiML1)[1,1:1000]
phi12 <- do.call(cbind, coefs2$PhiML1)[2,1:1000]
phi21 <- do.call(cbind, coefs2$PhiML2)[1,1:1000]
phi22 <- do.call(cbind, coefs2$PhiML2)[2,1:1000]

error_phi1 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML1))
error_phi2 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML2))

Beta11 <- do.call(cbind, coefs2$betasML1)[1,1:1000]
Beta12 <- do.call(cbind, coefs2$betasML1)[2,1:1000]
Beta21 <- do.call(cbind, coefs2$betasML2)[1,1:1000]
Beta22 <- do.call(cbind, coefs2$betasML2)[2,1:1000]

error_betasML2 <-  as.data.frame(do.call("rbind", coefs2$error_betasML2))
error_betasML1 <-  as.data.frame(do.call("rbind", coefs2$error_betasML1))


rho<-  as.data.frame(do.call("rbind", coefs2$rho))
error_rho <-  as.data.frame(do.call("rbind", coefs2$error_rho))



table10 <- cbind(c(
   #bias
   mean(phi11)-0.5,
   mean(phi12)-2,
   mean(phi21)-1,
   mean(phi22)-10,
   mean(Beta11)-2,
   mean(Beta12)-3,
   mean(Beta21)-10,
   mean(Beta22)-300,
   mean(rho$rho.rho1)-rho1),
   
   c(#lower limit
      mean(phi11)-0.5-1.96*mean(error_phi1$V1),
      mean(phi12)-2-1.96*mean(error_phi1$V2),
      mean(phi21)-1-1.96*mean(error_phi2$V1),
      mean(phi22)-10-1.96*mean(error_phi2$V2),
      mean(Beta11)-2-1.96*mean(error_betasML1$V1),
      mean(Beta12)-3-1.96*mean(error_betasML1$V2),
      mean(Beta21)-10-1.96*mean(error_betasML2$V1),
      mean(Beta22)-300-1.96*mean(error_betasML2$V2),
      mean(rho$rho.rho1)-rho1-1.96*mean(error_rho$rho.rho1)),
   
   c(#upper limit
      mean(phi11)-0.5+1.96*mean(error_phi1$V1),
      mean(phi12)-2+1.96*mean(error_phi1$V2),
      mean(phi21)-1+1.96*mean(error_phi2$V1),
      mean(phi22)-10+1.96*mean(error_phi2$V2),
      mean(Beta11)-2+1.96*mean(error_betasML1$V1),
      mean(Beta12)-3+1.96*mean(error_betasML1$V2),
      mean(Beta21)-10+1.96*mean(error_betasML2$V1),
      mean(Beta22)-300+1.96*mean(error_betasML2$V2),
      mean(rho$rho.rho1)-rho1+1.96*mean(error_rho$rho.rho1)))


n <- 300
rho1 <- 0.5

load(paste(n,rho1,sep="_",".RData"))

Sim <- rbind(results[1:1000])

coefs2 <- as.data.frame(do.call("rbind", Sim))


phi11 <- do.call(cbind, coefs2$PhiML1)[1,1:1000]
phi12 <- do.call(cbind, coefs2$PhiML1)[2,1:1000]
phi21 <- do.call(cbind, coefs2$PhiML2)[1,1:1000]
phi22 <- do.call(cbind, coefs2$PhiML2)[2,1:1000]

error_phi1 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML1))
error_phi2 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML2))

Beta11 <- do.call(cbind, coefs2$betasML1)[1,1:1000]
Beta12 <- do.call(cbind, coefs2$betasML1)[2,1:1000]
Beta21 <- do.call(cbind, coefs2$betasML2)[1,1:1000]
Beta22 <- do.call(cbind, coefs2$betasML2)[2,1:1000]

error_betasML2 <-  as.data.frame(do.call("rbind", coefs2$error_betasML2))
error_betasML1 <-  as.data.frame(do.call("rbind", coefs2$error_betasML1))


rho<-  as.data.frame(do.call("rbind", coefs2$rho))
error_rho <-  as.data.frame(do.call("rbind", coefs2$error_rho))



table11 <- cbind(c(
   #bias
   mean(phi11)-0.5,
   mean(phi12)-2,
   mean(phi21)-1,
   mean(phi22)-10,
   mean(Beta11)-2,
   mean(Beta12)-3,
   mean(Beta21)-10,
   mean(Beta22)-300,
   mean(rho$V1)-rho1),
   
   c(#lower limit
      mean(phi11)-0.5-1.96*mean(error_phi1$phi1),
      mean(phi12)-2-1.96*mean(error_phi1$phi2),
      mean(phi21)-1-1.96*mean(error_phi2$phi3),
      mean(phi22)-10-1.96*mean(error_phi2$phi4),
      mean(Beta11)-2-1.96*mean(error_betasML1$V1),
      mean(Beta12)-3-1.96*mean(error_betasML1$V2),
      mean(Beta21)-10-1.96*mean(error_betasML2$V1),
      mean(Beta22)-300-1.96*mean(error_betasML2$V2),
      mean(rho$V1)-rho1-1.96*mean(error_rho$V1)),
   
   c(#upper limit
      mean(phi11)-0.5+1.96*mean(error_phi1$phi1),
      mean(phi12)-2+1.96*mean(error_phi1$phi2),
      mean(phi21)-1+1.96*mean(error_phi2$phi3),
      mean(phi22)-10+1.96*mean(error_phi2$phi4),
      mean(Beta11)-2+1.96*mean(error_betasML1$V1),
      mean(Beta12)-3+1.96*mean(error_betasML1$V2),
      mean(Beta21)-10+1.96*mean(error_betasML2$V1),
      mean(Beta22)-300+1.96*mean(error_betasML2$V2),
      mean(rho$V1)-rho1+1.96*mean(error_rho$V1)))


n <- 300
rho1 <- 0.8

load(paste(n,rho1,sep="_",".RData"))

Sim <- rbind(results[1:1000])

coefs2 <- as.data.frame(do.call("rbind", Sim))


phi11 <- do.call(cbind, coefs2$PhiML1)[1,1:1000]
phi12 <- do.call(cbind, coefs2$PhiML1)[2,1:1000]
phi21 <- do.call(cbind, coefs2$PhiML2)[1,1:1000]
phi22 <- do.call(cbind, coefs2$PhiML2)[2,1:1000]

error_phi1 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML1))
error_phi2 <-  as.data.frame(do.call("rbind", coefs2$error_PhiML2))

Beta11 <- do.call(cbind, coefs2$betasML1)[1,1:1000]
Beta12 <- do.call(cbind, coefs2$betasML1)[2,1:1000]
Beta21 <- do.call(cbind, coefs2$betasML2)[1,1:1000]
Beta22 <- do.call(cbind, coefs2$betasML2)[2,1:1000]

error_betasML2 <-  as.data.frame(do.call("rbind", coefs2$error_betasML2))
error_betasML1 <-  as.data.frame(do.call("rbind", coefs2$error_betasML1))


rho<-  as.data.frame(do.call("rbind", coefs2$rho))
error_rho <-  as.data.frame(do.call("rbind", coefs2$error_rho))



table12 <- cbind(c(
   #bias
   mean(phi11)-0.5,
   mean(phi12)-2,
   mean(phi21)-1,
   mean(phi22)-10,
   mean(Beta11)-2,
   mean(Beta12)-3,
   mean(Beta21)-10,
   mean(Beta22)-300,
   mean(rho$rho.rho1)-rho1),
   
   c(#lower limit
      mean(phi11)-0.5-1.96*mean(error_phi1$V1),
      mean(phi12)-2-1.96*mean(error_phi1$V2),
      mean(phi21)-1-1.96*mean(error_phi2$V1),
      mean(phi22)-10-1.96*mean(error_phi2$V2),
      mean(Beta11)-2-1.96*mean(error_betasML1$V1),
      mean(Beta12)-3-1.96*mean(error_betasML1$V2),
      mean(Beta21)-10-1.96*mean(error_betasML2$V1),
      mean(Beta22)-300-1.96*mean(error_betasML2$V2),
      mean(rho$rho.rho1)-rho1-1.96*mean(error_rho$rho.rho1)),
   
   c(#upper limit
      mean(phi11)-0.5+1.96*mean(error_phi1$V1),
      mean(phi12)-2+1.96*mean(error_phi1$V2),
      mean(phi21)-1+1.96*mean(error_phi2$V1),
      mean(phi22)-10+1.96*mean(error_phi2$V2),
      mean(Beta11)-2+1.96*mean(error_betasML1$V1),
      mean(Beta12)-3+1.96*mean(error_betasML1$V2),
      mean(Beta21)-10+1.96*mean(error_betasML2$V1),
      mean(Beta22)-300+1.96*mean(error_betasML2$V2),
      mean(rho$rho.rho1)-rho1+1.96*mean(error_rho$rho.rho1)))





