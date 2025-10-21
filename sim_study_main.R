###           Simulations studies for                 ### 
### Incomplete U-Statistics of Equireplicate Designs: ###
### Berry-Esseen Bound and Efficient Construction     ###

source("sim_study_functions.R")

load("C10_list.Rdata") #CIFAR-10 images

library(MASS)

library(ggplot2)

library(doParallel)

## Remark 6, Section 4 ##

t_even <- system.time({

 equirep_even(n = 10^6,r = 10^2)

})

t_odd <- system.time({
  
  equirep_odd(n = 10^6+1,r = 10^2)
  
})

## Simulation Study 1 - CLT validation (see Theorem 3.8 and Corollary 3.10) ##

n_vec <- c(100,200,400,800,1600) 

r_vec <- c(1, ceiling(log(n_vec[1])), ceiling(log(n_vec[1])^2),ceiling(log(n_vec[1])^3),n_vec[1]/2,n_vec[1]-1)   

p <- 1    #Dimensionality (number of variables)

#rbf_sigma <- exact_medheur(d = p)

nsim <- 500 

mu_x <- rep(0, p)

eps <- 0

sigma_xy <- diag(p)

parallel::detectCores() 

ncores <- 100

cl <- makeCluster(ncores) 

registerDoParallel(cl)

getDoParWorkers() 

itpar <- ncores

results1 <- foreach(s = 1:itpar, .combine = 'c', .inorder = TRUE, 
                   .packages = c('ggplot2', 'MASS'), .errorhandling = 'pass') %dopar% { 
                     
                     tab_mat1 <- matrix(data = rep(NA_integer_,length(n_vec)*length(r_vec)),nrow = length(n_vec),ncol = length(r_vec))
                     
                     tab_mat2 <- matrix(data = rep(NA_integer_,length(n_vec)*length(r_vec)),nrow = length(n_vec),ncol = length(r_vec))
                     
                     tab_mat3 <- matrix(data = rep(NA_integer_,length(n_vec)*length(r_vec)),nrow = length(n_vec),ncol = length(r_vec))
                     
                     for (j in 1:length(n_vec)) {
                       
                       n <- n_vec[j]
                       
                       r_vec <- c(1, ceiling(log(n)), ceiling(log(n)^2),ceiling(log(n)^3), ceiling(n/2), n-1)   
                       
                       for (q in 1: length(r_vec)){
                         
                         r <- r_vec[q]
                        
                         D_size <- n*r/2
                         
                         D_equi <- equirep_even(n = n,r = r)
                         
                         z_vec <- rep(NA_integer_,nsim)
                         
                         z_vec_iid <- rep(NA_integer_,nsim)
                         
                         m_vec <- rep(NA_integer_,nsim)
                         
                         if (r==1) {
                           
                           sig_vec_iid <- rep(NA_integer_,nsim)
                           
                         }
                         
                         for (i in 1:nsim) {
                           
                           #timela <- system.time({
                             
                             set.seed(i+ s*nsim)
                             
                             # Sample from multivariate normal distributions
                             X <- mvrnorm(n, mu = mu_x, Sigma = sigma_xy)
                             Y <- mvrnorm(n, mu = eps, Sigma = sigma_xy)
                             
                             # Gaussian kernel
                             
                             #rbf_sigma <- median_heur_MMD(X = X,Y = Y)
                             
                             #hS_vec_equi <- mmd_iustat(X = X,Y = Y,D = D_equi[[2]],kernel_function = rbf_kernelpp, kernel_params = list(sigma = rbf_sigma))
                             
                             # Linear kernel
                             
                             hS_vec_equi <- mmd_iustat(X = X,Y = Y,D = D_equi[[2]],kernel_function = linear_kernelpp)
                             
                             m_vec[i] <- mean(hS_vec_equi)
                             
                             z_vec[i] <- (sqrt(D_size)*mean(hS_vec_equi))/sqrt(var(hS_vec_equi))
                             
                             if(r == 1) {
                               
                               sig_vec_iid[i] <- sqrt(var(hS_vec_equi))
                               
                               z_vec_iid[i] <- z_vec[i]
                               
                             } else {
                               
                               z_vec_iid[i] <- (sqrt(D_size)*mean(hS_vec_equi))/sig_vec_iid[i]
                               
                             }
                             
                         #  }
                        #   )
                           
                         }
                         
                         z_vec_mc <- m_vec/sqrt(var(m_vec))
 
                         tab_mat1[j,q] <- ks.test(z_vec, "pnorm", mean = 0, sd = 1)$statistic
                         
                         tab_mat2[j,q] <- ks.test(z_vec_iid, "pnorm", mean = 0, sd = 1)$statistic
                         
                         tab_mat3[j,q] <- ks.test(z_vec_mc, "pnorm", mean = 0, sd = 1)$statistic

                       }
                       
                     } 
                     
                     output <- list(
                       tab_s = tab_mat1,
                       tab_iid = tab_mat2,
                       tab_mc = tab_mat3
                     )
                     
                     list(output)
                     
                   }

stopCluster(cl)

save.image(file = "sim_study_one.Rdata")

## Simulation Study 2 - Minimum Variance (see Section 4: Efficient Construction of Equireplicate Designs) ##

n_vec <- c(100,200,400,800,1200,1600,2000)  

r_vec <- c(1, ceiling(log(n_vec[1])), ceiling(log(n_vec[1])^2),ceiling(log(n_vec[1])^3))   
  
p <- 1    #Dimensionality (number of variables)

ncomp <- 2

nsim <- 500 

mu_x <- rep(0, p)

eps <- 2 #to avoid sigma_1^2 equal to 0

sigma_xy <- diag(p)

parallel::detectCores()

ncores <- 100

cl <- makeCluster(ncores) 

registerDoParallel(cl)

getDoParWorkers() 

itpar <- ncores

results2 <- foreach(s = 1:itpar, .combine = 'c', .inorder = TRUE, 
                   .packages = c('ggplot2', 'MASS'), .errorhandling = 'pass') %dopar% { 
                     
                     tab_mat1 <- matrix(data = rep(NA_integer_,length(n_vec)*length(r_vec)),nrow = length(n_vec),ncol = length(r_vec))
                     
                     tab_mat2 <- matrix(data = rep(NA_integer_,length(n_vec)*length(r_vec)),nrow = length(n_vec),ncol = length(r_vec))
                     
                     for (j in 1:length(n_vec)) {
                       
                       n <- n_vec[j]
                       
                       r_vec <- c(1, ceiling(log(n)), ceiling(log(n)^2),ceiling(log(n)^3))   
                       
                       for (q in 1: length(r_vec)){
                         
                         r <- r_vec[q]
                         
                         D_size <- n*r/2
                         
                         D_equi <- equirep_even(n = n,r = r)
                        
                         set.seed(s)
                         
                         D_rand <- randes_full_k2(n = n,B = D_size)
                         
                         m_mat <- matrix(data = rep(NA_integer_,(ncomp)*nsim),nrow = (ncomp),ncol = nsim)
                         
                         for (i in 1:nsim) {
       
                           set.seed(i+ s*nsim)
                           
                           # Sample from multivariate normal distributions
                           X <- mvrnorm(n, mu = mu_x, Sigma = sigma_xy)
                           Y <- mvrnorm(n, mu = eps, Sigma = sigma_xy)
                           
                           # Gaussian kernel
                           
                           #rbf_sigma <- median_heur_MMD(X = X,Y = Y)
                           
                           #hS_vec_equi <- mmd_iustat(X = X,Y = Y,D = D_equi[[2]],kernel_function = rbf_kernelpp, kernel_params = list(sigma = rbf_sigma))
                           
                           #hS_vec_rand <- mmd_iustat(X = X,Y = Y,D = D_rand,kernel_function = rbf_kernelpp, kernel_params = list(sigma = rbf_sigma))
                           
                           # Linear case
                           
                           hS_vec_equi <- mmd_iustat(X = X,Y = Y,D = D_equi[[2]],kernel_function = linear_kernelpp)
                           
                           hS_vec_rand <- mmd_iustat(X = X,Y = Y,D = D_rand,kernel_function = linear_kernelpp)
                           
                           m_mat[1,i] <- mean(hS_vec_equi)
                           
                           m_mat[2,i] <- mean(hS_vec_rand)
                           
                         }
                         
                         tab_mat1[j,q] <- var(m_mat[1,])/var(m_mat[2,])
                         
                         tab_mat2[j,q] <- var(m_mat[2,])/var(m_mat[1,])
                         
                       }
                       
                     } 
                     
                     output <- list(
                       var_ratio_mat1 = tab_mat1,
                       var_ratio_mat2 = tab_mat2
                     )
                     
                     list(output)
                     
                   }

stopCluster(cl)

save.image(file = "sim_study_two.Rdata")

## Simulation Study 3 - CLT validation for k>2 (see Theorem 3.8) ##

eta <- function(x) 2^x

k <- 4

n_vec <- c(200,400,800,1600)  #choose n to be divisible by 4

r_vec <- c(1, 
           k*round(log(n_vec[1])/k), 
           k*round(log(n_vec[1])^2/k),
           k*round(log(n_vec[1])^3/k)
           )
          #must be multiple of k  

n_vec > 3*eta(k-1)*(eta(k-1)-eta(0))

ncomp <- 1 #number of competitors

p <- 1    # Dimensionality (number of variables)

nsim <- 500

rho <- 0

mu_x <- rep(0, p)

sigma_xy <- diag(p)

#rbf_sigmaX <- rbf_sigmaY <- exact_medheur(d = p,tau = sigma_xy)

parallel::detectCores()

ncores <- 100 
  
cl <- makeCluster(ncores) 

registerDoParallel(cl)

getDoParWorkers() 

itpar <- ncores

results3 <- foreach(s = 1:itpar, .combine = 'c', .inorder = TRUE, 
                   .packages = c('ggplot2', 'MASS'), .errorhandling = 'pass') %dopar% { 
                     
                     tab_mat1 <- matrix(data = rep(NA_integer_,length(n_vec)*length(r_vec)),nrow = length(n_vec),ncol = length(r_vec))
                     
                     tab_mat2 <- matrix(data = rep(NA_integer_,length(n_vec)*length(r_vec)),nrow = length(n_vec),ncol = length(r_vec))
                     
                     tab_mat3 <- matrix(data = rep(NA_integer_,length(n_vec)*length(r_vec)),nrow = length(n_vec),ncol = length(r_vec))
                     
                     for (j in 1:length(n_vec)) {
                       
                       n <- n_vec[j]
                       
                       r_vec <- c(1, 
                                  k*round(log(n)/k), 
                                  k*round((log(n)^2)/k),
                                  k*round((log(n)^3)/k)
                                 )
                       
                       for (q in 1: length(r_vec)){
                         
                         r <- r_vec[q]
                         
                         D_size <- n*r/k
                         
                         if (r==1) {
                           
                           D_equi <- matrix(1:n, ncol = k, byrow = TRUE)
                           
                           sig_mat_iid <- matrix(data = rep(NA_integer_,(ncomp)*nsim),nrow = (ncomp),ncol = nsim)
                           
                         } else {
                           
                           D <- cyclic_design(n = n,k = k,r = r,eta = eta)
                           
                           D_equi <- D[[2]]
                           
                         }
                         
                         z_mat <- matrix(data = rep(NA_integer_,(ncomp)*nsim),nrow = (ncomp),ncol = nsim)
                         
                         z_mat_iid <- matrix(data = rep(NA_integer_,(ncomp)*nsim),nrow = (ncomp),ncol = nsim)
                         
                         m_mat <- matrix(data = rep(NA_integer_,(ncomp)*nsim),nrow = (ncomp),ncol = nsim)
                         
                         for (i in 1:nsim) {
                           
                           set.seed(i+s*nsim)
                           
                           # Sample from multivariate normal distributions
                           X <- mvrnorm(n, mu = mu_x, Sigma = sigma_xy)
                           E <- mvrnorm(n, mu = mu_x, Sigma = sigma_xy)
                           
                           Y <- rho * sin(X) + sqrt(1 - rho^2) * E
                           
                           # Y <- rho * (X^2 - 1) + sqrt(1 - rho^2) * E #non-linear case
                           
                           # Gaussian kernel
                           
                           #rbf_sigmaX <- median_heur_HSIC(Z = X)
                           
                           #rbf_sigmaY <- median_heur_HSIC(Z = Y)
                           
                           #hS_vec_equi <- hsic_iustat(X = X,Y = Y,D = D_equi,kernel_function = rbf_kernelpp, kernel_paramsX = list(sigma = rbf_sigmaX),kernel_paramsY = list(sigma = rbf_sigmaY))
                           
                           # Linear kernel
                           
                           hS_vec_equi <- hsic_iustat(X = X,Y = Y,D = D_equi,kernel_function = linear_kernelpp)
                           
                           m_mat[1,i] <- mean(hS_vec_equi)
                           
                           z_mat[1,i] <- (sqrt(D_size)*mean(hS_vec_equi))/sqrt(var(hS_vec_equi))
                           
                           if(r == 1) {
                             
                             sig_mat_iid[1, i] <- sqrt(var(hS_vec_equi))
                             
                             z_mat_iid[1,i] <- z_mat[1,i]
                             
                           } else {
                             
                             z_mat_iid[1,i] <- (sqrt(D_size)*mean(hS_vec_equi))/sig_mat_iid[1,i]
                             
                           }
                           
                         }
                         
                         z_mat_mc <- m_mat[1,]/sqrt(var(m_mat[1,]))
                         
                         tab_mat1[j,q] <- ks.test(z_mat[1,], "pnorm", mean = 0, sd = 1)$statistic
                         
                         tab_mat2[j,q] <- ks.test(z_mat_iid[1,], "pnorm", mean = 0, sd = 1)$statistic
                         
                         tab_mat3[j,q] <- ks.test(z_mat_mc, "pnorm", mean = 0, sd = 1)$statistic
                         
                       }
                       
                     } 
                     
                     output <- list(
                       tab_s = tab_mat1,
                       tab_iid = tab_mat2,
                       tab_mc = tab_mat3
                     )
                     
                     list(output)
                     
                   }

stopCluster(cl)

save.image(file = "sim_study_three.Rdata")

## Simulation Study 4 - Minimum Variance (see Section 4: Efficient Construction of Equireplicate Designs) ##

eta <- function(x) 2^x

k <- 4

n_vec <- c(200,400,800,1600)  #choose n to be divisible by 4

n_vec > 3*eta(k-1)*(eta(k-1)-eta(0))

r_vec <- c(1, 
           k*round(log(n_vec[1])/k), 
           k*round(log(n_vec[1])^2/k),
           k*round(log(n_vec[1])^3/k)
          )#must be multiple of k  

p <- 1    #Dimensionality (number of variables)

ncomp <- 2

nsim <- 500 

mu_x <- rep(0, p)

rho <- 0.5

sigma_xy <- diag(p)

parallel::detectCores()

ncores <- 100

cl <- makeCluster(ncores) 

registerDoParallel(cl)

getDoParWorkers() 

itpar <- ncores

results4 <- foreach(s = 1:itpar, .combine = 'c', .inorder = TRUE, 
                    .packages = c('ggplot2', 'MASS'), .errorhandling = 'pass') %dopar% { 
                      
                      tab_mat1 <- matrix(data = rep(NA_integer_,length(n_vec)*length(r_vec)),nrow = length(n_vec),ncol = length(r_vec))
                      
                      tab_mat2 <- matrix(data = rep(NA_integer_,length(n_vec)*length(r_vec)),nrow = length(n_vec),ncol = length(r_vec))
                      
                      for (j in 1:length(n_vec)) {
                        
                        n <- n_vec[j]
                        
                        r_vec <- c(1, 
                                   k*round((log(n))/k), 
                                   k*round((log(n)^2)/k),
                                   k*round((log(n)^3)/k)
                        )#must be multiple of k  
                        
                        for (q in 1: length(r_vec)){
                          
                          r <- r_vec[q]
                          
                          D_size <- n*r/k
                          
                          if (r==1) {
                            
                            D_equi <- matrix(1:n, ncol = k, byrow = TRUE)
                            
                          } else {
                            
                            D <- cyclic_design(n = n,k = k,r = r,eta = eta)
                            
                            D_equi <- D[[2]]
                            
                          }
                     
                          set.seed(s)
                          
                          #D_rand <- randes_full_k4(n = n,B = D_size) #comp.intensive
                          
                          D_rand <- t(replicate(D_size, sort(sample.int(n, k, replace = FALSE))))
                          
                          m_mat <- matrix(data = rep(NA_integer_,(ncomp)*nsim),nrow = (ncomp),ncol = nsim)
                          
                          for (i in 1:nsim) {
                            
                            set.seed(i+s*nsim)
                            
                            # Sample from multivariate normal distributions
                            X <- mvrnorm(n, mu = mu_x, Sigma = sigma_xy)
                            E <- mvrnorm(n, mu = mu_x, Sigma = sigma_xy)
                            
                            Y <- rho * sin(X) + sqrt(1 - rho^2) * E
                            
                            # Y <- rho * (X^2 - 1) + sqrt(1 - rho^2) * E #non-linear case
                            
                            # Gaussian kernel
                            
                            rbf_sigmaX <- median_heur_HSIC(Z = X)
                            
                            rbf_sigmaY <- median_heur_HSIC(Z = Y)
                            
                            hS_vec_equi <- hsic_iustat(X = X,Y = Y,D = D_equi,kernel_function = rbf_kernelpp, kernel_paramsX = list(sigma = rbf_sigmaX), kernel_paramsY = list(sigma = rbf_sigmaY))
                            
                            hS_vec_rand <- hsic_iustat(X = X,Y = Y,D = D_rand,kernel_function = rbf_kernelpp, kernel_paramsX = list(sigma = rbf_sigmaX), kernel_paramsY = list(sigma = rbf_sigmaY))
                            
                            # Linear case
                            
                            #hS_vec_equi <- hsic_iustat(X = X,Y = Y,D = D_equi,kernel_function = linear_kernelpp)
                            
                            #hS_vec_rand <- hsic_iustat(X = X,Y = Y,D = D_rand,kernel_function = linear_kernelpp)
                            
                            m_mat[1,i] <- mean(hS_vec_equi)
                            
                            m_mat[2,i] <- mean(hS_vec_rand)
                            
                          }
                          
                          tab_mat1[j,q] <- var(m_mat[1,])/var(m_mat[2,])
                          
                          tab_mat2[j,q] <- var(m_mat[2,])/var(m_mat[1,])
                          
                        }
                        
                      } 
                      
                      output <- list(
                        var_ratio_mat1 = tab_mat1,
                        var_ratio_mat2 = tab_mat2
                      )
                      
                      list(output)
                      
                    }

stopCluster(cl)

save.image(file = "sim_study_four.Rdata")

## Simulation Study 5 - Real Data Example (Power, Type I error, Computational burden) ##

lab_CIFAR10 <- c("airplane", "automobile", "bird","cat","deer","dog","frog","horse","ship","truck")

n_vec <- c(100,200,400,800,1600) 

r_vec <- c(ceiling(log(n_vec[1])), ceiling(log(n_vec[1])^2))  
    
kx <- c(3,4,8, 9) 

ky <- c(0,1,5,7) 

c <- 1 #0  

alpha <- 0.05
  
n_sim <- 500 
  
n_perm <- 1000

parallel::detectCores()

ncores <- 125

cl <- makeCluster(ncores) 

registerDoParallel(cl)

getDoParWorkers() 

resultsA <- foreach(s = 251:n_sim, .combine = 'c', .inorder = TRUE, 
                    .packages = c('ggplot2', 'MASS'), .errorhandling = 'pass') %dopar% { 
                      
                      set.seed(s)
                      
                      p_mat_equi <- matrix(data = rep(NA_integer_,length(n_vec)*length(r_vec)),nrow = length(n_vec),ncol = length(r_vec))
                      
                      t_mat_equi <- matrix(data = rep(NA_integer_,length(n_vec)*length(r_vec)),nrow = length(n_vec),ncol = length(r_vec))
                      
                      p_mat_rand <- matrix(data = rep(NA_integer_,length(n_vec)*length(r_vec)),nrow = length(n_vec),ncol = length(r_vec))
                      
                      t_mat_rand <- matrix(data = rep(NA_integer_,length(n_vec)*length(r_vec)),nrow = length(n_vec),ncol = length(r_vec))
                      
                      for (j in 1:length(n_vec)) {
                        
                        n <- n_vec[j]
                        
                        df <- build_XY_multi(P = P,kx = kx,ky = ky,n = n,c = c)
                        
                        sig <- median_heur_MMD(X = df$X,Y = df$Y)
                        
                        r_vec <- c(ceiling(log(n)), ceiling(log(n)^2))   
                        
                        for (q in 1: length(r_vec)){
                          
                          r <- r_vec[q]
                          
                          t_equi <- system.time({
                            
                            out_equi <- MMD_equirep_tst(X = df$X,Y = df$Y,r = r,kernel_function = rbf_kernelpp,
                                                    kernel_params = list(sigma = sig),alpha = alpha)
                            
                          })
                          
                          t_rand <- system.time({
                            
                            out_rand <- MMD_rand_perm_rbf_tst(X = df$X ,Y = df$Y,r = r,sig = sig, n_perm = n_perm,alpha = alpha)
                            
                          })
                          
                          p_mat_equi[j,q] <- out_equi$z_test
                          
                          p_mat_rand[j,q] <- out_rand$z_perm_test
                          
                          t_mat_equi[j,q] <- as.numeric(t_equi[3])
                            
                          t_mat_rand[j,q] <- as.numeric(t_rand[3])
                          
                        }
                        
                      } 
                      
                      output <- list(
                        p_mat_equi = p_mat_equi,
                        p_mat_rand = p_mat_rand,
                        t_mat_equi = t_mat_equi,
                        t_mat_rand = t_mat_rand
                      )
                      
                      list(output)
                      
                    }

stopCluster(cl)

rm(P)

save.image(file = "real_data_partA.Rdata")
