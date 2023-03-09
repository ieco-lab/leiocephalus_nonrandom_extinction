#' Calculate the mean and variance as species are removed
#'
#' The function \code{stepwise_extinction} uses a set extinction order and repeatedly calculates the mean and variance as species are lost
#'
#'@param x A dataframe or tibble containing all columns of data
#'@param ordering Name of column in x that orders extinctions from 1 to n that sets the order of species to remove sequentially
#'@param trait Name of column in x that contains the data for the trait of interest
#'@param spp Name of column in x that contains the names of species for comparison
#'@param rescale Logical that determines if all variance and mean values through time are rescaled by the starting value
#'@param sliding_eval Logical that determines if sliding evaluation should be used, whereby each loss step calculates based on the \code{eval_step} loss step values before it by obtaining the mean across those steps and the current step
#'@param eval_step Integer value that determines the number of previous loss steps to average across for each loss step
#'
#'@return Dataframe (tibble) that contains columns that include the number of species, the mean at each loss step, the variance at each loss step, and the names of the species that are lost at that step
#'
#'@export

stepwise_extinction <- function(x, ordering, trait, spp, rescale = TRUE, sliding_eval = FALSE, eval_step = 1){

  #tripwire for sliding window evaluation if the size of eval is greater than the number of data points
  if(sliding_eval == TRUE & (eval_step > length(x[[spp]]) | eval_step < 1 | eval_step%%1 != 0)){
    stop("'eval_step' must be a nonzero integer smaller than 'n_ex' to evaluate with a sliding window")
  }

  #check that ordering is a column in data
  if(isFALSE(as.character(c(ordering, trait, spp)) %in% colnames(x))){
    stop("'ordering' must be a column within 'data'")
  }

  #build a simplified version of x
  coord <- tibble(id = as.integer(x[[ordering]]),
                  spp = x[[spp]],
                  trait = as.double(x[[trait]])
  )

  #initialize the results objects
  extinction_results <- as_tibble(as.data.frame(matrix(nrow = length(coord$id)+1, ncol = 4)))
  colnames(extinction_results) <- c("n_spp", "avg", "var", "spp_rm")


  #get the means and vars to start and we will modify depending on the arguments after
  extinction_results$avg <- coord %>%
    arrange(id) %>%
    slide(.x = .$trait, .f = ~mean(.x, na.rm = T), .before = 0, .after = Inf) %>%
    unlist(.) %>%
    c(., NA)

  extinction_results$var <- coord %>%
    arrange(id) %>%
    slide(.x = .$trait, .f = ~var(.x, na.rm = T), .before = 0, .after = Inf) %>%
    unlist(.) %>%
    c(., NA)

  #sliding eval == FALSE but still using slider package
  if(sliding_eval == FALSE && rescale == FALSE){
    #empty because we do not need to make any changes (we can remove this for brevity later)
  }
  else if(sliding_eval == FALSE && rescale == TRUE){
    #extra rescale step!
    extinction_results <- extinction_results %>%
      mutate(avg = avg - avg[1], var = var - var[1])

  }
  else if(sliding_eval == TRUE && rescale == TRUE){
    #now take the difference between a step and previous n steps based on eval_step and obtain the mean (only really matters for steps > 1)
    extinction_results <- extinction_results %>%
      mutate(
        avg = slide(.x = .$avg, .f = ~diff(.x), .before = eval_step) %>% #gets the difference between a loss step and the n ones before it based on eval_step size
          ifelse((sapply(., length) == 0), NA, .) %>% #changes things that are numeric(0) to NA
          lapply(X = ., FUN = function(x) mean(x, na.rm = T)) %>% #calculates the means of the differences returned, this generates NaN values for NA only elements
          rapply(., f = function(x) ifelse(is.nan(x), 0, x), how = "replace") %>% #replaces NaN values with zero (0)
          unlist(.), #unlists the results, as slider::slide() results are returned as a list. we will match the results to their original column in the dataframe!
        var = slide(.x = .$var, .f = ~diff(.x), .before = eval_step) %>%
          ifelse((sapply(., length) == 0), NA, .) %>%
          lapply(X = ., FUN = function(x) mean(x, na.rm = T)) %>%
          rapply(., f = function(x) ifelse(is.nan(x), 0, x), how = "replace") %>%
          unlist(.)
      )

  }
  else if(sliding_eval == TRUE && rescale == FALSE){
    #TESTING
    #extinction_results$avg %>%
    #  slide(.x = ., .f = ~.x, .before = eval_step) %>% #gets the difference between a loss step and the n ones before it based on eval_step size
    #  print(.)
    #
    #extinction_results$avg %>%
    #  slide(.x = ., .f = ~mean(.x, na.rm = T), .before = eval_step) %>% #gets the difference between a loss step and the n ones before it based on eval_step size
    #      #ifelse((sapply(., length) == 0), NA, .) %>% #changes things that are numeric(0) to NA
    #      #lapply(X = ., FUN = function(x) mean(x, na.rm = T)) %>% #calculates the means of the differences returned, this generates NaN values for NA only elements
    #      #rapply(., f = function(x) ifelse(is.nan(x), 0, x), how = "replace") %>% #replaces NaN values with zero (0)
    #      #unlist(.)
    #  print(.)

    #implement
    extinction_results <- extinction_results %>%
      mutate(
        avg = slide(.x = .$avg, .f = ~mean(.x, na.rm = T), .before = eval_step) %>%
          unlist(.),
        var = slide(.x = .$var, .f = ~mean(.x, na.rm = T), .before = eval_step) %>%
          unlist(.)
      )

  }


  extinction_results$n_spp <- length(unique(coord$spp)):0

  #get the species removed at each step
  extinction_results$spp_rm <- coord %>%
    arrange(id) %>%
    slide(.x = .$spp, .f = ~.x, .before = 0, .after = 0) %>%
    unlist(.) %>%
    c(NA, .)

  #output
  return(extinction_results)
}
