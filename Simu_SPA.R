
library(parallel)
library(survival)
library(coxphf)
#library(miranda) #generate SNP matrix fast 
#library(microbenchmark)
library(Rcpp)
library(Rfast)
library(SPACox)
Rcpp::sourceCpp('simulation.cpp')

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
Surv_gen <- function(n, scale, shape, cens_prob){
  
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

alpha <- 1e-3

#MAF <- c(0.02, 0.05, 0.1, 0.2, 0.5)
MAF <- c(0.5, 0.2, 0.1, 0.05, 0.02)
N <- c(350) #sample size
p_total <- 1e9 #total SNP 
K <- 1e4 #simulation steps
p <- p_total/K
scale <- 1
SHAPE <- c(1, 1.5)
CENS_PROB <- c(0.25, 0.5)


scenario <- length(N)* length(MAF)* length(SHAPE)* length(CENS_PROB)


# Main part and output results -------------------------------------------------
dir.create("/home/y014g/Simulation_sig_SPA_1billion_1e-3")
path <- "/home/y014g/Simulation_sig_SPA_1billion_1e-3/"
# dir.create("/Users/jocelyn/Desktop/Biolab/Project/Simulation code/Simulation_SPA")
# path <- "/Users/jocelyn/Desktop/Biolab/Project/Simulation code/Simulation_SPA/"

type1error <- data.frame(matrix(nrow = scenario, ncol = 5,
                                dimnames = list(NULL, c("N", "MAF", "Shape_para","Cens_Prob","SPACox"))))

Error_num_SPA <- 0

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
          
          tryCatch({
            sig_iteration <- rep(0, 2)
            #sig_iteration <- data.frame(matrix(0, nrow = 1, ncol = 5, dimnames = list(NULL, c("Score", "Wald", "LikeRatio" ,"SPACox", "Firth"))))
            
            #Generate survival data
            Surv_data <- Surv_gen(n, scale, shape, cens_prob)
          })
          
          # 3. SPACox
          tryCatch({
            mull <- SPACox_Null_Model(Surv(time,status) ~ 1, data = Surv_data, pIDs = 1:n, gIDs = 1:n);
            
            colnames(X) <- 1:p;
            #colnames(X) <- 1:p
            rownames(X) <- 1:n;
            p_spacox <- SPACox(mull, X)
            sig_iteration[1] <- length(which(p_spacox[,3] < alpha))
          },
          error = function(e) {sig_iteration[2] <- sig_iteration[2] +1},
          warning = function(w) {sig_iteration[2] <- sig_iteration[2] +1})
          
          #print(sig_iteration)
          #return(list("sig" = sig))
          return(list(sig_iteration))
        }
        , mc.cores = 50) # loop K steps mclapply
        
        sig <- data.frame(matrix(unlist(senario_result), nrow=K, byrow=T, 
                                 dimnames = list(NULL, c("SPACox", "Err"))))
        
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
        type1error[iteration, 5] <- sum(sig[,1]) / p_total
        
        iteration <- iteration + 1
        
        # Record number of test errors
        
        Error_num_SPA <- Error_num_SPA + sum(sig[,2])
        
        # Save all p values
        FileName <- paste(simu_info, '.csv',sep = "")
        write.csv(sig[,1], paste(path, FileName), row.names = TRUE)
        
      } # censoring prob *2
    } # shape *2
  } # loop MAF *6
} # sample size *3

write.csv(type1error, paste(path, 'Type1_error.csv'), row.names = TRUE)

simu_time <- proc.time() - start

