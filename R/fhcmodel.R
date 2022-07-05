#' Fit Flexible-Hazards Cure Model
#' @description Fits a flexible-hazards cure model. Tied failure times are dealt with Breslow approximation.
#' @param data a data.frame with no missing values contains survival time, event status, any continuous variable(s) and/or dummy variable(s) of categorical variable(s).
#' The categorical variables with more than 2 levels need to be converted into dummy variables before fitting the function.
#' @param event_time the name of observed survival time in the data.
#' @param event_status the name of event status in the data.
#' @param id the name of patient id in the data.
#' @param beta_variable the name(s) of variable(s) defined as long-term covariate(s) in the model, which cannot be NULL.
#' @param gamma_variable the name(s) of variable(s) defined as short-term covariate(s) in the model. By default gamma_variable = NULL.
#' If there is no short-term covariate that needs to be defined, gamma should be set as NULL and the model becomes a proportional hazards cure model.
#' @param coef a logical option for obtaining point estimates of regression parameters. By default coef = TRUE.
#' @param se a logical option for obtaining standard errors of regression parameters. By default se = TRUE.
#'        If it is TRUE, the program returns standard error by perturbation method. If it is set to FALSE, the program only returns estimated regression parameters.
#' @param max.int maximum number of iterations. The default is 200.
#' @param se.pert the number of runs for perturbation. The default is 100.
#' @param no_cores the number of cores used while computing the standard errors by applying the perturbation. The default is 7.
#' @param absconverge.par absolute difference \eqn{(current-previous)} for regression parameters as the convergence criteria. The default is 1e-6.
#' @param relconverge.par relative difference \eqn{[(current-previous)/previous]} for regression parameters as the convergence criteria. The default is 1e-4.
#' @param absconverge.F0t the average of absolute differences for \eqn{F_0(t)} as the convergence criteria. The default is 1e-6.
#' @param relconverge.F0t the average of relative differences for \eqn{F_0(t)} as the convergence criteria. The default is 1e-3.
#'
#' @return a list containing results of the fit. The following items coef, iter, and data are returned when the argument coef=TRUE.
#' The item std.err and perturb are returned when the argument se=TRUE.
#'     \item{coef}{estimated regression parameters (beta and gamma) if coef=TRUE}
#'     \item{iter}{the number of iterations used to complete the point estimation if coef=TRUE}
#'     \item{data}{the final data including the estimated \eqn{F_0(t)} and \eqn{f_0(t)} if coef=TRUE}
#'     \item{std.err}{standard errors of regression paramters if se=TRUE}
#'     \item{perturb}{a matrix containing the perturbation estimates, whose number of rows is se.pert, if se=TRUE}
#' @export
#'
#' @examples data(fhc_dat)
#' result <- fhcmodel(data=fhc_dat,event_status = "event",event_time="time",id="id",beta_variable = c("age","sex"),gamma_variable = c("age","sex"),se=F)
#' result <- fhcmodel(data=fhc_dat,event_status = "event",event_time="time",id="id",beta_variable = c("age","sex"),se=T)
#' result$coef
#' @import foreach
#' @import parallel
#' @import survival
#' @importFrom maxLik maxNR

#' @importFrom zoo na.locf

fhcmodel <- function(data, event_time, event_status, id, beta_variable, gamma_variable=NULL, coef=T, se=T, max.int=200,
                     se.pert=100, no_cores=7,
                     absconverge.par=1e-6, relconverge.par=1e-4, absconverge.F0t=1e-6, relconverge.F0t=1e-3){

  cat("Program is running. Please be patient.","\n")

  n <- nrow(data)

  point_est <- function(data,pert){

    cat("Point estimation starts.","\n")

    colnames(data)[colnames(data)==event_status] <- "event"
    colnames(data)[colnames(data)==event_time] <- "T_new"

    Z_var <- beta_variable
    if (!is.null(gamma_variable)){
      X_var <- gamma_variable
    } else{
      X_var <- beta_variable
    }

    ### Functions for estimation
    #### Step 2: Update beta, gamma, and beta0 ###
    # Function 1
    beta_gamma_comp <- function(par,beta0){

      beta_gamma_loglik <- function(par){

        beta <- par[1:length(beta_variable)]
        if (!is.null(gamma_variable)){
          gamma <- par[(length(beta_variable)+1):(length(beta_variable)+length(gamma_variable))]
        } else{
          gamma <- rep(0,length(X_var))
        }

        data_temp <- as.matrix(data)
        Z_i <- data_temp[,Z_var]
        X_i <- data_temp[,X_var]
        betaZ <- Z_i%*%as.matrix(beta,nrow=length(beta))
        gammaX <- X_i%*%as.matrix(gamma,nrow=length(gamma))
        pert_weight_i <- data_temp[,"p_weight"]
        event_i <- data_temp[,"event"]

        failure_part <- event_i*(beta0+betaZ+gammaX+(exp(gammaX)-1)*log(data_temp[,"base_cdf"])+log(data_temp[,"base_f"]))
        failure_part[which(is.na(failure_part))] <- 0
        LL <- sum(pert_weight_i*(failure_part-exp(beta0)*exp(betaZ)*(data_temp[,"base_cdf"]^exp(gammaX))))

        return(LL)
      }

      beta_gamma_est <-optim(par=par,fn=beta_gamma_loglik,method = "Nelder-Mead",control=list("fnscale"=-1,maxit=10000), hessian=F)

      if (beta_gamma_est$convergence==0){
        beta_gamma_k <- beta_gamma_est$par
        beta <- beta_gamma_k[1:length(beta_variable)]
        if (!is.null(gamma_variable)){
          gamma <- beta_gamma_k[(length(beta_variable)+1):(length(beta_variable)+length(gamma_variable))]
        } else {
          gamma <- rep(0,length(Z_var))
        }
      } else{
        if(pert==F){
          stop(paste("beta and gamma does not converge at iteration",iteration))
        }else{
          beta <- "beta does not converge"
          gamma <- "gamma does not converge"
        }
      }

      return(list(beta,gamma))
    }

    # Function 2
    beta0_comp <- function(beta,gamma){
      pert_weight <- data$p_weight
      data_temp <- as.matrix(data)
      Z_i <- data_temp[,Z_var]
      X_i <- data_temp[,X_var]
      betaZ <- Z_i%*%as.matrix(beta,nrow=length(beta))
      gammaX <- X_i%*%as.matrix(gamma,nrow=length(gamma))
      beta_0 <- log(pert_weight%*%data_temp[,"event"]/(t(pert_weight*exp(betaZ))%*%(data_temp[,"base_cdf"]^(exp(gammaX)))))
      return(as.vector(beta_0))
    }

    ######################### Step 3: Update F0(t) #########################
    # Function 3
    f_est <- function(gamma,beta,beta0){

      t <- data$T_new[data$event==1]

      data_temp <- as.matrix(data)
      # for patients who are at risk
      Z_i <- lapply(unique(t),function(x) data_temp[data_temp[,"T_new"]>=x,Z_var])
      X_i <- lapply(unique(t),function(x) data_temp[data_temp[,"T_new"]>=x,X_var])
      base_F_i <- lapply(unique(t), function(x) data_temp[data_temp[,"T_new"]>=x,"base_cdf"])

      # for patients who are at risk and fail in the future
      X_i_failure <- lapply(unique(t),function(x) data_temp[data_temp[,"T_new"]>=x & data_temp[,"event"]==1,X_var])
      base_F_i_failure <- lapply(unique(t), function(x) data_temp[data_temp[,"T_new"]>=x & data_temp[,"event"]==1,"base_cdf"])

      # weight
      q_w_i <- sapply(unique(t),function(x) sum(data$p_weight[data$event==1&data$T_new==x])) # with event
      q_w_j <- lapply(unique(t),function(x) data_temp[data_temp[,"T_new"]>=x,"p_weight"])# at risk
      q_w_failure <- lapply(unique(t),function(x) data_temp[data_temp[,"T_new"]>=x & data_temp[,"event"]==1, "p_weight"])#at risk and fail in the future

      beta <- matrix(beta,nrow=length(beta))
      gamma <- matrix(gamma,nrow=length(gamma))

      # denominator
      deno <- sapply(seq(unique(t)),function(x)
        t(q_w_j[[x]]*exp(beta0)*exp(Z_i[[x]]%*%beta)*exp(X_i[[x]]%*%gamma))%*%((base_F_i[[x]])^(exp(X_i[[x]]%*%gamma)-1))
        -t(q_w_failure[[x]]*(exp(X_i_failure[[x]]%*%gamma)-1))%*%(1/base_F_i_failure[[x]]))

      #obtain the value of alpha
      alpha_est <- function(alpha){
        sum(q_w_i/(deno+alpha))-1
      }

      alpha <- "error"
      times <- 2
      while(alpha=="error" & times < 5){
        range_root <- lapply(seq(-10^times,10^times),function(x) c(x,(x+1)))
        alpha_result <- sapply(seq(range_root),function(i)tryCatch(uniroot(alpha_est,range_root[[i]],extendInt = "yes")$root,error=function(x) return("error")))
        alpha <- ifelse(sum(alpha_result=="error")<length(range_root),"root","error")
        if (alpha=="root"){
          alpha_select <- unique(round(as.numeric(alpha_result[alpha_result!="error"]),4))

          for (i in seq(alpha_select)){
            check <- as.numeric(alpha_select[i])
            if (abs(sum(q_w_i/(deno+check))-1) <1e-4 & sum(q_w_i/(deno+check)<0)==0){
              alpha_select[i] <- check
            } else{
              alpha_select[i] <- "error"
            }
          }
          if (sum(alpha_select!="error")>0){
            alpha <- as.numeric(alpha_select[alpha_select!="error"])
            if (length(alpha)>1){
              alpha <- alpha[1]
            }
          }else{alpha <- "error"}
        }
        times <- times +1
      }

      new_base_f <- q_w_i/(deno+alpha)
      new_base_cdf <- cumsum(new_base_f)
      event_base_f <- unlist(sapply(seq(length(unique(t))),function(x) rep(new_base_f[x],sum(data$T_new[data$event==1]==unique(t)[x]))))
      event_base_cdf <- unlist(sapply(seq(length(unique(t))),function(x) rep(new_base_cdf[x],sum(data$T_new[data$event==1]==unique(t)[x]))))

      remove <- which(colnames(data)%in%c("base_cdf","base_f"))
      data <- data[,-remove]

      data$base_f <- NA
      data$base_cdf <- NA
      data$base_f[data$event==1] <- event_base_f
      data$base_f[data$event==0] <- 0
      data$base_cdf[data$event==1] <- event_base_cdf

      cen <- which(data$event==0)
      death <- which(data$event==1)
      #check beginning
      if (sum(cen<death[1])!=0){
        pos <- cen[cen<death[1]]
        data$base_cdf[pos] <- 0
      }
      data$base_cdf <- na.locf(data$base_cdf)

      return(list(new_base_cdf,data))
    }

    ### Start estimation ###
    #### Step 1: find initial values ###
    data <- data[order(data$T_new),]
    t <- data$T_new[data$event==1]
    formula.exp <- paste("Surv(T_new, event)", paste(beta_variable, collapse=" + "), sep=" ~ ")
    initial <- coxph(as.formula(formula.exp), data = data, weights = p_weight, ties="breslow")

    H0 <- basehaz(initial, centered=FALSE) # considered ties using breslow estimator
    H0 <- H0[H0[, 2] %in% unique(t), ]

    if (nrow(H0)!=length(unique(t))){
      compute.h.base <- function(k){
        Z_j <- data[data$T_new>=t[k],beta_variable]
        h_base <- 1/sum(exp(t(initial$coefficients)%*%t(Z_j)))
        return(h_base)
      }
      h_base <- rep(NA,length(t))
      for (i in 1:length(t)){
        h_base[i] <- compute.h.base(i)
      }
      H_base <- cumsum(h_base)
      H0 <- data.frame(hazard=H_base,time=t)
    }

    #initial values of beta
    beta <- initial$coefficients
    #initial value of beta_0
    beta_0 <- log(H0$hazard[H0$time==max(H0$time)])
    #initial values of baseline CDF
    base_cdf <- sapply(1:nrow(H0),function(x) H0$hazard[x]/H0$hazard[length(H0$hazard)])
    #initial values of baseline pdf
    base_f <- c(base_cdf[1],sapply(2:length(base_cdf), function(x) base_cdf[x]-base_cdf[x-1]))
    base <- data.frame(T_new=H0$time,base_cdf=base_cdf,base_f=base_f)

    remove <- which(colnames(data)%in%c("base_cdf","base_f"))
    if (length(remove)!=0){
      data <- data[,-remove]
    }
    data[data$event==1,c("base_cdf","base_f")] <- merge(data[data$event==1,],base,by="T_new")[,c("base_cdf","base_f")]
    data$base_f[data$event==0] <- 0

    cen <- which(data$event==0)
    death <- which(data$event==1)
    #check beginning
    if (sum(cen<death[1])!=0){
      pos <- cen[cen<death[1]]
      data$base_cdf[pos] <- 0
    }
    data$base_cdf <- na.locf(data$base_cdf)

    ## first iteration
    iteration <- 1
    # Step 2 gamma and beta_0
    beta_gamma_k1 <- beta_gamma_comp(par=c(beta,rep(0,length(Z_var))),beta0=beta_0)
    beta_k1 <- beta_gamma_k1[[1]]
    gamma_k1 <- beta_gamma_k1[[2]]
    beta_0_k1 <- beta0_comp(beta = beta_k1,gamma = gamma_k1)

    # Step 3 baseline CDF
    base_f_0 <- f_est(gamma=as.vector(gamma_k1),beta=beta_k1,beta0=beta_0_k1)
    base_cdf_k1 <- base_f_0[[1]]
    data <- base_f_0[[2]]

    cat("Iteration 1 is done.", "\n")

    # print(beta_0_k1)
    # print(beta_k1)
    # print(gamma_k1)

    # second iteration
    iteration <- 2
    # step 2
    beta_gamma_k2 <- beta_gamma_comp(par=c(beta_k1,gamma_k1),beta0=beta_0_k1)
    beta_k2 <- beta_gamma_k2[[1]]
    gamma_k2 <- beta_gamma_k2[[2]]
    beta_0_k2 <- beta0_comp(beta = beta_k2,gamma = gamma_k2)

    # Step 3 baseline CDF
    base_f_0 <- f_est(gamma=as.vector(gamma_k2),beta=beta_k2,beta0=beta_0_k2)
    base_cdf_k2 <- base_f_0[[1]]
    data <- base_f_0[[2]]

    cat("Iteration 2 is done.", "\n")

    # print(beta_0_k2)
    # print(beta_k2)
    # print(gamma_k2)

    ## Keep doing iterations until convergence
    beta_diff <- min(c(sum(abs((beta_k2-beta_k1)/beta_k1) > relconverge.par),sum((abs(beta_k2-beta_k1) > absconverge.par))))
    beta_0_diff <- min(c(abs((beta_0_k2-beta_0_k1)/beta_0_k1) > relconverge.par, abs(beta_0_k2-beta_0_k1) > absconverge.par))
    base_cdf_diff <- min(c(mean(abs((base_cdf_k2-base_cdf_k1)/base_cdf_k1)) > relconverge.F0t, mean(abs((base_cdf_k2-base_cdf_k1))) > absconverge.F0t))

    if(!is.null(gamma_variable)){
      gamma_diff <- min(c(sum(abs((gamma_k2-gamma_k1)/gamma_k1) > relconverge.par),sum((abs(gamma_k2-gamma_k1) > absconverge.par))))
      deter <- ( beta_diff !=0 | beta_0_diff !=0 | gamma_diff!=0 | base_cdf_diff !=0)
    }else{
      deter <- ( beta_diff !=0 | beta_0_diff !=0 | base_cdf_diff !=0)
    }

    while( deter & iteration <= max.int){

      iteration  <- iteration + 1
      beta_k1 <- beta_k2
      beta_0_k1 <- beta_0_k2
      gamma_k1 <- gamma_k2
      base_cdf_k1 <- base_cdf_k2

      # step 2
      beta_gamma_k2 <- beta_gamma_comp(par=c(beta_k1,gamma_k1),beta0=beta_0_k1)
      beta_k2 <- beta_gamma_k2[[1]]
      gamma_k2 <- beta_gamma_k2[[2]]
      beta_0_k2 <- beta0_comp(beta = beta_k2,gamma = gamma_k2)

      # Step 3 baseline CDF
      base_f_0 <- f_est(gamma=gamma_k2,beta=beta_k2,beta0=beta_0_k2)
      base_cdf_k2 <- base_f_0[[1]]
      data <- base_f_0[[2]]

      cat("Iteration", iteration, "is done.", "\n")

      # print(iteration)
      # print(beta_0_k2)
      # print(beta_k2)
      # print(gamma_k2)

      beta_diff <- min(c(sum(abs((beta_k2-beta_k1)/beta_k1) > relconverge.par),sum((abs(beta_k2-beta_k1) > absconverge.par))))
      beta_0_diff <- min(c(abs((beta_0_k2-beta_0_k1)/beta_0_k1) > relconverge.par, abs(beta_0_k2-beta_0_k1) > absconverge.par))
      base_cdf_diff <- min(c(mean(abs((base_cdf_k2-base_cdf_k1)/base_cdf_k1)) > relconverge.F0t, mean(abs((base_cdf_k2-base_cdf_k1))) > absconverge.F0t))

      if(!is.null(gamma_variable)){
        gamma_diff <- min(c(sum(abs((gamma_k2-gamma_k1)/gamma_k1) > relconverge.par),sum((abs(gamma_k2-gamma_k1) > absconverge.par))))
        deter <- ( beta_diff !=0 | beta_0_diff !=0 | gamma_diff!=0 | base_cdf_diff !=0)
      }else{
        deter <- ( beta_diff !=0 | beta_0_diff !=0 | base_cdf_diff !=0)
      }

    }

    if (!is.null(gamma_variable)){
      result.coef <- matrix(c(beta_0_k2,beta_k2,gamma_k2,iteration),nrow=1)
      colnames(result.coef) <- c("intercept",paste0("beta_",beta_variable),paste0("gamma_",gamma_variable),"num.iter")
    }else{
      result.coef <- matrix(c(beta_0_k2,beta_k2,iteration),nrow=1)
      colnames(result.coef) <- c("intercept",paste0("beta_",beta_variable),"num.iter")
    }

    data <- data[order(data[,id]),]
    data <- subset(data, select = -p_weight )
    colnames(data)[colnames(data)=="event"] <- event_status
    colnames(data)[colnames(data)=="T_new"] <- event_time

    return(list(coef=result.coef,data=data))

  }

  # point estimation
  if (coef){
    data$p_weight <- 1
    result.coef <- point_est(data,pert=FALSE)
    para_coef <- result.coef$coef[1:(1+length(beta_variable)+length(gamma_variable))]
    names(para_coef) <- colnames(result.coef$coef)[1:(1+length(beta_variable)+length(gamma_variable))]
    iter <- result.coef$coef[2+length(beta_variable)+length(gamma_variable)]
    data_final <- result.coef$data
  }

  # standard deviation
  if (se){

    set.seed(1)
    seeds_pert <- round(runif(se.pert)*10^3)
    k <- length(seeds_pert)

    no_cores <- no_cores
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)
    cat("Standard error estimation starts.Please be patient.","\n")
    result.pert <- foreach(i=icount(k),.packages =c("survival","maxLik","zoo"),.combine=rbind,.errorhandling = "remove") %dopar% {
      #write(paste("Starting perturbation run",i,"\n"),file="log.txt",append=TRUE)
      set.seed(seeds_pert[i])
      data$p_weight <- rexp(n,rate=1)
      one_run <- point_est(data,pert=TRUE)
      one_run$coef[1:(1+length(beta_variable)+length(gamma_variable))]}
    stopCluster(cl)

    para_se <- apply(result.pert,2,sd)
    if (!is.null(gamma_variable)){
      names(para_se) <- c("intercept",paste0("beta_",beta_variable),paste0("gamma_",gamma_variable))
    } else {
      names(para_se) <- c("intercept",paste0("beta_",beta_variable))
    }

    cat("Standard error estimation is done.")
  }

  if (coef==T & se==T){
    result <- list(coef=para_coef,iter=iter,std.err=para_se,perturb=result.pert,data=data_final)
  } else if(coef==T & se==F){
    result <- list(coef=para_coef,iter=iter,data=data_final)
  } else if (coef==F & se==T){
    result <- list(std.err=para_se,perturb=result.pert)
  }

  return(result)

}
