#'
#'
#'

midpoint_resample = function(x){
  round(runif(n = 1, min = x-0.5, max = x+0.5),3)
}
midpoint_resample_vec = Vectorize(midpoint_resample)