#' Estimate cure rates for subgroups
#' @description Estimates a cure rate for a subgroup of patients with specified characteristics.
#' @param object an object of fhcmodel.
#' @param z_value value(s) of the long-term covariate(s) in the order specified at the argument beta_variable of the fhcmodel().

#'
#' @return an estimated cure rate
#' @export
#'
#' @examples
#' result <- fhcmodel(data=fhc_dat,event_status = "event",event_time="time",,id="id",beta_variable = c("age","sex"),gamma_variable = c("age","sex"),se=F)
#' estcure(object=result,estvalue=c(0,1)) # cure rate for a subgroup population with age as the mean and female
#' estcure(object=result,estvalue=c(-20,0)) # cure rate for a subgroup population with 20 years younger than the mean age and male

estcure <- function(object, z_value){
  coef <- object$coef
  exp(-exp(coef[1]+coef[2:(1+length(z_value))]%*%z_value))
}
