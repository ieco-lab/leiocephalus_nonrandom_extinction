#' Sample Ties in Extinction Order
#'
#' The function \code{sample_extinction} handles ties in extinction order by sampling from the ties randomly to create a clean sequential order
#'
#'@param cluster_order A vector of numbers that indicate the order of species going extinct (1 is the first to go extinct) that includes duplicates that indicate clustered or unknown order of extinction
#'
#'@return A vector of the sampled extinction order that is in clean sequential order and preserves the order of the original input (can be \code{cbind} to the input)
#'
#'@export

sample_extinction <- function(cluster_order){
  #initial counter
  counter <- 0

  #dataframe the original input and the fixed
  cluster_object <- data.frame(
    cluster_input = cluster_order,
    cluster_fixed = rep(NA, times = length(cluster_order))
  )

  #the loop to populate
  for(a in sort(unique(cluster_order))){
    #check if there is a cluster at this stage
    if(length(cluster_order[cluster_order == a]) > 1){
      #sample the order
      recode_order <- counter + sample(x = 1:length(cluster_order[cluster_order == a]))
      #add the new order to the new column
      cluster_object$cluster_fixed[cluster_order == a] <- recode_order
      #update counter
      counter <- max(cluster_object$cluster_fixed, na.rm = T)

    }
    #case for lacking clustering
    else{
      #add 1 to counter
      counter <- counter + 1
      #recode the new order  for this
      cluster_object$cluster_fixed[cluster_object$cluster_input == a] <- counter
    }
  }
  return(cluster_object$cluster_fixed)
}
