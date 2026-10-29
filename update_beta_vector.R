update_beta_vector <- function(
  beta_vec,
  beta_new,
  update_time
){
  stopifnot(is.numeric(beta_vec))
  stopifnot(length(beta_new) == 1)
  stopifnot(update_time >= 1)
  stopifnot(update_time <= length(beta_vec))

  beta_vec[update_time:length(beta_vec)] <- beta_new

  beta_vec
}