require(HDInterval)

interval1 <- vector("list")
interval_bci1 <- vector("list")
p <- vector("list")
q <- as.numeric()

for (j in 1:100){
   coefs <- results[[j]][[6]]
   coefs2 <- cbind(
      phi11 = do.call(cbind, coefs$PhiML1_b)[1,],
      phi12 = do.call(cbind, coefs$PhiML1_b)[2,],
      phi21 = do.call(cbind, coefs$PhiML2_b)[1,],
      phi22 = do.call(cbind, coefs$PhiML2_b)[2,],
      Beta11 = do.call(cbind, coefs$betasML1_b)[1,],
      Beta12 = do.call(cbind, coefs$betasML1_b)[2,],
      Beta21 = do.call(cbind, coefs$betasML2_b)[1,],
      Beta22 = do.call(cbind, coefs$betasML2_b)[2,],
      rho = do.call(cbind, coefs$rho_b)[1,]
   )
   
   
   sorted <- apply(coefs2,2,sort)
   
   #Interval HDI for Betas
   interval1[j] <- list(apply(coefs2,2,hdi)) 
   
   
   #Interval bias-corrected for phi and rho
   for (i in 1:9){
      p <- qnorm(sum(coefs2[,i] < apply(coefs2,2,mean)[i])/1000)
      q[i] <- list(rbind(lower = sorted[pnorm(2*p-1.96)*1000,i],
                         upper = sorted[pnorm(2*p+1.96)*1000,i]))
   }
   
   interval_bci1[j] <- list(cbind(q[[1]],q[[2]],q[[3]],q[[4]],q[[5]],q[[6]],q[[7]],q[[8]],q[[9]])
   )
}


phi11_c <- c(1:100)
phi12_c <- c(1:100)
phi21_c <- c(1:100)
phi22_c <- c(1:100)
beta11_c <- c(1:100)
beta12_c <- c(1:100)
beta21_c <- c(1:100)
beta22_c <- c(1:100)
rho_c <- c(1:100)


for(i in 1:100){
   
   if (interval_bci1[[i]][,1][1]<=0.5 & interval_bci1[[i]][,1][2]>=0.5)
   {
      phi11_c[i]<- 1
   }
   else{ 
      phi11_c[i] <-  0
   }
   
   
   if (interval_bci1[[i]][,2][1]<=2 & interval_bci1[[i]][,2][2]>=2)
   {
      phi12_c[i]<- 1
   }
   else{ 
      phi12_c[i] <-  0
   }
   
   if (interval_bci1[[i]][,3][1]<=1 & interval_bci1[[i]][,3][2]>=1)
   {
      phi21_c[i]<- 1
   }
   else{ 
      phi21_c[i] <-  0
   }
   
   if (interval_bci1[[i]][,4][1]<=10 & interval_bci1[[i]][,4][2]>=10)
   {
      phi22_c[i]<- 1
   }
   else{ 
      phi22_c[i] <-  0
   }
   
   
   if (interval_bci1[[i]][,5][1]<=2 & interval_bci1[[i]][,5][2]>=2)
   {
      beta11_c[i]<- 1
   }
   else{ 
      beta11_c[i] <-  0
   }
   
   if (interval_bci1[[i]][,6][1]<=3 & interval_bci1[[i]][,6][2]>=3)
   {
      beta12_c[i]<- 1
   }
   else{ 
      beta12_c[i] <-  0
   }
   
   if (interval_bci1[[i]][,7][1]<=10 & interval_bci1[[i]][,7][2]>=10)
   {
      beta21_c[i]<- 1
   }
   else{ 
      beta21_c[i] <-  0
   }
   
   
   if (interval_bci1[[i]][,8][1]<=300 & interval_bci1[[i]][,8][2]>=300)
   {
      beta22_c[i]<- 1
   }
   else{ 
      beta22_c[i] <-  0
   }
   
   
   if (interval_bci1[[i]][,9][1]<=rho1 & interval_bci1[[i]][,9][2]>=rho1)
   {
      rho_c[i]<- 1
   }
   else{ 
      rho_c[i] <-  0
   }
   
}

coverage <- rbind(round(sum(beta11_c)/100,2),
                  round(sum(beta12_c)/100,2),
                  round(sum(beta21_c)/100,2),
                  round(sum(beta22_c)/100,2),
                  round(sum(phi11_c)/100,2),
                  round(sum(phi12_c)/100,2),
                  round(sum(phi21_c)/100,2),
                  round(sum(phi22_c)/100,2),
                  round(sum(rho_c)/100,2))



save.image("50_-0.8boot_.RData")
