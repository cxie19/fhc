#' Compute the estimated survival function for a subgroup
#' @description Estimates the survival function for a subgroup of patients with specified characteristics.
#' @param object an object of fhcmodel.
#' @param z_value value(s) of the long-term covariate(s) in the order specified at the argument beta_variable of the fhcmodel().
#' @param x_value value(s) of the short-term covariate(s) in the order specified at the argument gamma_variable of the fhcmodel().
#' @param event_time the name of observed survival time in the data.
#'
#' @return a data.frame with sorted survival times and estimated survival probabilities.
#' @export
#'
#' @examples
#' result <- fhcmodel(data=fhc_dat,event_status = "event",event_time="time",id="id",
#'                    beta_variable = c("age","sex"),gamma_variable = c("age","sex"),se=F)
#' survival_female <- estsurvival(object=result,z_value=c(mean(fhc_dat$age),1),
#'                                x_value=c(mean(fhc_dat$age),1),event_time="time")
#' survival_male <- estsurvival(object=result,z_value=c(mean(fhc_dat$age),0),
#'                              x_value=c(mean(fhc_dat$age),0),event_time="time")
#' plot(survival_male$time,survival_male$survival,ylim=c(0,1),xlim=c(0,40),xlab="Time",
#'      ylab="Survival probability",type="l",col="red")
#' lines(survival_female$time,survival_female$survival,col="blue")
#' legend(25,1,c("Male","Female"),cex=1,lty=c(1,1),col=c("red","blue"),horiz=F,bty="n")

estsurvival <- function(object,z_value,x_value,event_time){
  coef <- object$coef
  data <- object$data
  F_0t <- data$base_cdf[order(data[,event_time])]
  est <- exp(-as.vector(exp(coef[1]+coef[2:(1+length(z_value))]%*%z_value))*F_0t^as.vector(exp(coef[(length(coef)-length(x_value)+1):length(coef)]%*%x_value)))
  return(data.frame(time=sort(data[,event_time]),survival=est))
}
