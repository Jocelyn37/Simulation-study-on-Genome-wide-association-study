#### Final version before power cox model 

library(parallel)
library(survival)
library(coxphf)
#library(miranda) #generate SNP matrix fast 
#library(microbenchmark)
library(Rcpp)
library(Rfast)
library(SPACox)
Rcpp::sourceCpp('simulation.cpp')

onestep <- function(u) {
  
  # sig <- data.frame(matrix(0, nrow = 1, ncol = 5,
  #                         dimnames = list(NULL, c("Score", "Wald", "LikeRatio" ,"SPACox", "Firth"))))
  
  
  sig_iteration <- rep(0, 10)
  #sig_iteration <- data.frame(matrix(0, nrow = 1, ncol = 5, dimnames = list(NULL, c("Score", "Wald", "LikeRatio" ,"SPACox", "Firth"))))
  
  #Generate survival data
  Surv_data <- Surv_gen(n, maf, scale, shape, cens_prob)
  
  #Generate survival analysis tests outputs
  # 1. Score test
  tryCatch({
    p_score <- score_notie(Surv_data$time, Surv_data$status, X)
    sig_iteration[1] <- length(which(p_score < alpha))
    
    reduced_snp <- which(p_score < threshold)
    reduced_X <- X[, c(reduced_snp)]
    
  },error = function(e) {sig_iteration[6] <- sig_iteration[6] + 1})
  
  
  # 2. Wald test & Likelihood-Ratio test & Firth correction based test
  # Test each SNP
  tryCatch({
    
    for(j in reduced_snp){
      #for(j in 1:p){
      
      # j <- 12
      #print(j)
      x <- X[,j]
      
      cox <- summary(coxph(Surv(time, status) ~ x, data = Surv_data))
      wald <- cox$waldtest[3]
      #print(wald)
      loglike <- cox$logtest[3]
      # print(".....")
      # print(loglike)
      if (!is.na(wald))
        sig_iteration[2] <- sig_iteration[2] + 1*(wald < alpha)
      else
        sig_iteration[7] <- sig_iteration[7] +1
      if (!is.na(loglike))
        sig_iteration[3] <- sig_iteration[3] + 1*(loglike < alpha)
      else
        sig_iteration[8] <- sig_iteration[8] +1
      
      #
      # #Problem: when x are all 0, a convergence check and issues a warning in case of non-convergence.
      # Error: Loglik converged before variable  1 ; coefficient may be infinite.
      
      if (1 %in% x | 2 %in% x){
        #cox_Firth <- summary(coxphf(Surv(time, status) ~ x,  data = Surv_data, alpha = alpha))
        Surv_data$x <- x
        cox_Firth <- coxphftest(as.formula("Surv(time, status) ~ x"), test = as.formula("~x"), data = Surv_data)
        Firth <- cox_Firth$prob
        if (!is.na(Firth))
          sig_iteration[5] <- sig_iteration[5] + 1*(Firth < alpha)
        else
          sig_iteration[10] <- sig_iteration[10] +1 
      } #Firth test
    } # test each reduced_snp
  }) #trycatch
  # 3. SPACox
  
  tryCatch({
    mull <- SPACox_Null_Model(Surv(time,status) ~ 1, data = Surv_data, pIDs = 1:n, gIDs = 1:n);
    
    colnames(reduced_X) <- reduced_snp;
    #colnames(X) <- 1:p
    rownames(reduced_X) <- 1:n;
    p_spacox <- SPACox(mull, reduced_X)
    sig_iteration[4] <- length(which(p_spacox[,3] < alpha))
  },
  
  error = function(e) {sig_iteration[9] <- sig_iteration[9] +1},
  warning = function(w) {sig_iteration[9] <- sig_iteration[9] +1})
  
  #print(sig_iteration)
  #return(list("sig" = sig))
  return(list(sig_iteration))
}


# Fast version of Score test ----------------------------------------------
score_notie <- function(time,status,X,covariates=NULL) {
  #  X <- X+0.1
  n<-length(time)
  q <- ncol(covariates)
  
  IO <- Rfast::Order(time) #order time from small to large
  time <- time[rev(IO)] # time are now ordered from biggest to smallest(50)
  status <- status[rev(IO)] # order corresponding status(50)
  X <- X[rev(IO),] # order X 
  
  sn <- status/(1:n) 
  
  X2 <- X^2
  sX2 <- Rfast::colCumSums(X2)
  term1 <- Rfast::colsums(sX2*sn)
  
  sX <- Rfast::colCumSums(X)
  mX <- sX/(1:n)
  mX2 <- mX^2
  term2 <- Rfast::colsums(mX2*status)
  vr <- term1-term2  #calculate variance
  
  fit <- coxph(Surv(time,status)~1)
  r <- fit$residuals
  cv <- cov(X,r) #test statistics
  
  stats <- cv^2*(n-1)^2/vr  #nomalized statistics
  return(pchisq(stats,df=1,lower.tail=FALSE)) #p value
}


# Choose upper bound of Uniform censoring ------------------------------------
calc_cens_weibull <- function(c, lambda, k) {
  integrand <- function(x) {
    pweibull(x, k, 1/lambda,lower.tail=FALSE)
  }
  
  integrate(integrand, 0, c)$value * (1/c)
}

calc_up_weibull <- function(p, lambda, k, max = 100) {
  # scale parameter: lambda
  # shape parameter: k
  # p: probability of censoring 
  c0 <- 0
  c1 <- max
  pnold <- calc_cens_weibull(c1, lambda, k)
  diff <- 1-calc_cens_weibull(c1, lambda, k)
  
  while (diff > 1e-6) {
    cnew <- (c0+c1)/2
    pn <- calc_cens_weibull(cnew, lambda, k)
    if (pn <= p) {
      c1 <- (c0+c1)/2
    } else {
      c0 <- (c0+c1)/2 
    }
    diff <- abs(pn - p)
    pnold <- pn
  }
  return(cnew)
}


# Simulate Survival data -------------------------------------------------------
Surv_gen <- function(n, maf, scale, shape, cens_prob){
  
  #survival time(Weibull)
  coef = 0 # null hypothesis
  surv_time = (-log(runif(n)) / (scale * exp(coef)))^(1/shape)
  
  #censoring time(uniform)
  cens_up <- calc_up_weibull(p = cens_prob, lambda = scale, k = shape)
  cens_time <- runif(n = n, min = 0, max = cens_up)
  
  # time
  time <- pmin(surv_time, cens_time)
  
  # event status
  status <- as.numeric(surv_time <= cens_time)
  
  data.frame(time, status)
}

# Set Initial Value ----------------------------------------------------------

alpha <- 5e-8
threshold <- alpha*100

#MAF <- c(0.02, 0.05, 0.1, 0.2, 0.5)
MAF <- c(0.5, 0.2, 0.1, 0.05, 0.02)
N <- c(350) #sample size
p_total <- 1e6 #total SNP 
K <- 100 #simulation steps
p <- p_total/K
scale <- 1
SHAPE <- c(1, 1.5)
CENS_PROB <- c(0.25, 0.5)

n <- 10
p <- 20
maf <- 0.25
X <- cpp_gensnp(n, p, maf)
scale <- 1
shape <- 1.5
cens_prob<- 0.25
Surv_data <- Surv_gen(n, maf, scale, shape, cens_prob)


x <- X[,1]

cox <- summary(coxph(Surv(time, status) ~ x, data = Surv_data))
wald <- cox$waldtest[3]


scenario <- length(N)* length(MAF)* length(SHAPE)* length(CENS_PROB)


# Main part and output results -------------------------------------------------
dir.create("/home/y014g/Simulation_sig")
path <- "/home/y014g/Simulation_sig/"
# dir.create("/Users/jocelyn/Desktop/Biolab/Project/Simulation code/Simulation_sig")
# path <- "/Users/jocelyn/Desktop/Biolab/Project/Simulation code/Simulation_sig/"

type1error <- data.frame(matrix(nrow = scenario, ncol = 9,
                                dimnames = list(NULL, c("N", "MAF", "Shape_para","Cens_Prob", "Score", "Wald", "LikeRatio" ,"SPACox", "Firth"))))

Error_num_Score <- 0 # record error
Error_num_wald <- 0 
Error_num_likelihood <- 0
Error_num_SPA <- 0
Error_num_Firth <- 0

iteration <- 1

#set seed
set.seed(1)

start <- proc.time()

for(n in N){
  
  for(maf in MAF){
    
    #Generate SNP matrix(n*p matrix)
    X <- cpp_gensnp(n, p, maf)
    
    for(shape in SHAPE){
      
      for(cens_prob in CENS_PROB){
        
        # Simulation steps
        senario_result <- mclapply(1:K,function(u) {
          
          # sig <- data.frame(matrix(0, nrow = 1, ncol = 5,
          #                         dimnames = list(NULL, c("Score", "Wald", "LikeRatio" ,"SPACox", "Firth"))))
          
          tryCatch({
            sig_iteration <- rep(0, 10)
            #sig_iteration <- data.frame(matrix(0, nrow = 1, ncol = 5, dimnames = list(NULL, c("Score", "Wald", "LikeRatio" ,"SPACox", "Firth"))))
            
            #Generate survival data
            Surv_data <- Surv_gen(n, maf, scale, shape, cens_prob)})
          
          #Generate survival analysis tests outputs
          # 1. Score test
          tryCatch({
            p_score <- score_notie(Surv_data$time, Surv_data$status, X)
            sig_iteration[1] <- length(which(p_score < alpha))
            
            reduced_snp <- which(p_score < threshold)
            reduced_X <- X[, c(reduced_snp)]
            
          },error = function(e) {sig_iteration[6] <- sig_iteration[6] + 1})
          
          
          # 2. Wald test & Likelihood-Ratio test & Firth correction based test
          # Test each SNP
          tryCatch({
            
            for(j in reduced_snp){
              #for(j in 1:p){
              
              # j <- 12
              #print(j)
              x <- X[,j]
              
              cox <- summary(coxph(Surv(time, status) ~ x, data = Surv_data))
              wald <- cox$waldtest[3]
              #print(wald)
              loglike <- cox$logtest[3]
              # print(".....")
              # print(loglike)
              if (!is.na(wald))
                sig_iteration[2] <- sig_iteration[2] + 1*(wald < alpha)
              else
                sig_iteration[7] <- sig_iteration[7] +1
              if (!is.na(loglike))
                sig_iteration[3] <- sig_iteration[3] + 1*(loglike < alpha)
              else
                sig_iteration[8] <- sig_iteration[8] +1
              
              #
              # #Problem: when x are all 0, a convergence check and issues a warning in case of non-convergence.
              # Error: Loglik converged before variable  1 ; coefficient may be infinite.
              
              if (1 %in% x | 2 %in% x){
                #cox_Firth <- summary(coxphf(Surv(time, status) ~ x,  data = Surv_data, alpha = alpha))
                Surv_data$x <- x
                cox_Firth <- coxphftest(as.formula("Surv(time, status) ~ x"), test = as.formula("~x"), data = Surv_data)
                Firth <- cox_Firth$prob
                if (!is.na(Firth))
                  sig_iteration[5] <- sig_iteration[5] + 1*(Firth < alpha)
                else
                  sig_iteration[10] <- sig_iteration[10] + 1 
              } #Firth test
            } # test each reduced_snp
          }) #trycatch
          
          # 3. SPACox
          tryCatch({
            mull <- SPACox_Null_Model(Surv(time,status) ~ 1, data = Surv_data, pIDs = 1:n, gIDs = 1:n);
            
            colnames(reduced_X) <- reduced_snp;
            #colnames(X) <- 1:p
            rownames(reduced_X) <- 1:n;
            p_spacox <- SPACox(mull, reduced_X)
            sig_iteration[4] <- length(which(p_spacox[,3] < alpha))
          },
          error = function(e) {sig_iteration[9] <- sig_iteration[9] +1},
          warning = function(w) {sig_iteration[9] <- sig_iteration[9] +1})
          
          #print(sig_iteration)
          #return(list("sig" = sig))
          return(list(sig_iteration))
        }
        , mc.cores = 50) # loop K steps mclapply
        
        sig <- data.frame(matrix(unlist(senario_result), nrow=K, byrow=T, 
                                 dimnames = list(NULL, c("Score", "Wald", "LikeRatio" ,"SPACox", "Firth", "Err_S", "Err_W", "Err_L", "Err_SPA","Err_F"))))
        
        # simu_info <- data.frame(
        #   MAF = maf,
        #   shape_para = shape,
        #   censoring_rate = cens_prob)
        
        simu_info <- paste(paste("Sam", n, sep =""),
                           paste("Sh", shape, sep =""),
                           paste("Cens", cens_prob, sep =""),
                           paste("MAF", maf, sep =""),
                           sep = "-")
        
        # type 1 error
        type1error[iteration, "N"] <- n
        type1error[iteration, "MAF"] <- maf
        type1error[iteration, "Shape_para"] <- shape
        type1error[iteration, "Cens_Prob"] <- cens_prob
        type1error[iteration, 5:9] <- apply(sig[,c(1:5)], 2, sum) / p_total
        
        iteration <- iteration + 1
        
        # Record number of test errors
        
        Error_num_Score <- Error_num_Score + apply(sig["Err_S"], 2, sum)
        Error_num_wald <- Error_num_wald + apply(sig["Err_W"], 2, sum) 
        Error_num_likelihood <- Error_num_likelihood + apply(sig["Err_L"], 2, sum)
        Error_num_SPA <- Error_num_SPA + apply(sig["Err_SPA"], 2, sum)
        Error_num_Firth <- Error_num_Firth + apply(sig["Err_F"], 2, sum)
        
        # Save all p values
        FileName <- paste(simu_info, '.csv',sep = "")
        write.csv(sig[,c(1:5)], paste(path, FileName), row.names = TRUE)
        
      } # censoring prob *2
    } # shape *2
  } # loop MAF *6
} # sample size *3

write.csv(type1error, paste(path, 'Type1_error.csv'), row.names = TRUE)

simu_time <- proc.time() - start

