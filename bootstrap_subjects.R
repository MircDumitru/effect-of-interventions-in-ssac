bootstrap_subjects <- function(
  dat
){
  ids <- unique(dat$id)
  sampled_ids <- sample(
    ids,
    size = length(ids),
    replace = TRUE
  )
  boot_dat <- map_dfr(
    seq_along(sampled_ids),
    function(i){
      tmp <- dat[dat$id == sampled_ids[i], ]
      tmp$id <- i
      tmp
    }
  )
  boot_dat
}