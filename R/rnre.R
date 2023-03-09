#' Simulates extinctions and calculates mean and variance for predicted loss
#'
#' The function \code{rnre} simulates predicted patterns of extinction for a set of species with a trait of interest, a number of which go extinct. For each function call, a single sequence of loss is calculated for random, directional (may be two), disruptive, and stabilizing extinction.
#'
#'@param trait_x A sequence of trait values to use for each species. Each value is considered to be the mean trait value for a species.
#'@param n_ex Sets the number of species to go extinct. If not specified, the number of trait values in \code{trait_x} is the default.
#'@param sides Determines if loss is simulated for a specific tail or both tails for directional extinction
#'@param rescale Logical that determines if all variance and mean values through time are rescaled by the starting value (default set to TRUE)
#'@param sliding_eval Logical that determines if sliding evaluation should be used, whereby each loss step calculates based on the \code{eval_step} loss step values before it by obtaining the mean across those steps and the current step (default set to FALSE)
#'@param eval_step Integer value that determines the number of previous loss steps to average across for each loss step (default set to 1)
#'
#'@return Dataframe (tibble) that contains columns that include the number of species, the mean at each loss step, the variance at each loss step, and the type of type of simulated extinction. One instance of each type of extinction is run per function call and the resulting simulations are stacked in the same dataframe.
#'
#'@export


rnre <- function(trait_x, n_ex = length(trait_x), sides = c("BOTH", "LEFT", "RIGHT"), rescale = TRUE, sliding_eval = FALSE, eval_step = 1){

  #match sides argument against valid input
  sides <- match.arg(sides)

  #tripwire for sliding window evaluation if the size of eval is greater than the number of data points
  if(sliding_eval == TRUE & (eval_step > n_ex | eval_step < 1 | eval_step%%1 != 0)){
    stop("'eval_step' must be a nonzero integer smaller than 'n_ex' to evaluate with a sliding window")
  }

  #######################################################
  #Initialize starting objects
  #######################################################

  #make coord object with IDs to make picking easier
  coord <- tibble(id = as.integer(1:length(trait_x)),
                  x = as.double(trait_x)
                  )

  #calculate the starting mean and variance
  start_1trait <- as_tibble(as.data.frame(matrix(nrow = 1, ncol = 3)))
  colnames(start_1trait) <- c("n_spp", "avg", "var")

  #calculate values to fill in
  start_1trait$n_spp <- nrow(coord)
  start_1trait$avg <- mean(coord$x, na.rm = T)
  start_1trait$var <- var(x = coord$x, na.rm = T)

  #set species to go extinct for n_ex loss events for each of the different "patterns"

  ########################################
  #1. random

  #initialize storage obj for each rep
  sim_1trait_ran <- as_tibble(as.data.frame(matrix(nrow = n_ex + 1, ncol = 3)))
  colnames(sim_1trait_ran) <- c("n_spp", "avg", "var")

  sim_1trait_ran <- sim_1trait_ran %>%
    mutate(n_spp = as.integer(n_spp), avg = as.double(avg), var = as.double(var))


  #get the lost species and generate an order
  #sample of points to remove
  rm_ran <- sample(x = coord$id, size = n_ex)


  ########################################
  #2. directional

  if(sides == "LEFT" || sides == "RIGHT"){
    #directional holder
    sim_1trait_dir <- as_tibble(as.data.frame(matrix(nrow = n_ex + 1, ncol = 3)))
    colnames(sim_1trait_dir) <- c("n_spp", "avg", "var")

    sim_1trait_dir <- sim_1trait_dir %>%
      mutate(n_spp = as.integer(n_spp), avg = as.double(avg), var = as.double(var))

    #get the ID's
    #1D version: need to choose which side of the tail based on user input
    #tail_choice <- sample(c(TRUE, FALSE), size = 1)

    #ifelse to get the points we need to pull TRUE is remove the left tail, FALSE is remove the right tail
    if(sides == "LEFT"){
      #get the ID's to remove
      rm_dir <-  coord %>%
        arrange(x) %>%
        dplyr::select(id) %>%
        head(id, n = n_ex) %>%
        unlist(.)

      #get rid of old names
      names(rm_dir) <- NULL

    } else if(sides == "RIGHT"){
      #get the ID's to remove
      rm_dir <-  coord %>%
        arrange(desc(x)) %>%  #note that there was a syntax change in jumping to R v.4.0 for tail(), or so it seems... if we keep the arrange call here as for left, then we can just use tail(.$id, n = n_ex) in the line with the head() call
        dplyr::select(id) %>%
        head(id, n = n_ex) %>%
        unlist(.)

      #get rid of old names
      names(rm_dir) <- NULL
    }

  } else if(sides == "BOTH"){
    #left tail holder
    sim_1trait_dirl <- as_tibble(as.data.frame(matrix(nrow = n_ex + 1, ncol = 3)))
    colnames(sim_1trait_dirl) <- c("n_spp", "avg", "var")

    sim_1trait_dirl <- sim_1trait_dirl %>%
      mutate(n_spp = as.integer(n_spp), avg = as.double(avg), var = as.double(var))

    #right tail holder
    sim_1trait_dirr <- as_tibble(as.data.frame(matrix(nrow = n_ex + 1, ncol = 3)))
    colnames(sim_1trait_dirr) <- c("n_spp", "avg", "var")

    sim_1trait_dirr <- sim_1trait_dirr %>%
      mutate(n_spp = as.integer(n_spp), avg = as.double(avg), var = as.double(var))

    #get the ID's to remove
    rm_dirl <-  coord %>%
      arrange(x) %>%
      dplyr::select(id) %>%
      head(id, n = n_ex) %>%
      unlist(.)


    #get rid of old names
    names(rm_dirl) <- NULL

    #get the ID's to remove
    rm_dirr <-  coord %>%
      arrange(desc(x)) %>%
      dplyr::select(id) %>%
      head(id, n = n_ex) %>%
      unlist(.)

    #get rid of old names
    names(rm_dirr) <- NULL
  }


  ####################################
  #3. disruptive

  #initialize storage obj for each rep
  sim_1trait_dis <- as_tibble(as.data.frame(matrix(nrow = n_ex + 1, ncol = 3)))
  colnames(sim_1trait_dis) <- c("n_spp", "avg", "var")

  sim_1trait_dis <- sim_1trait_dis %>%
    mutate(n_spp = as.integer(n_spp), avg = as.double(avg), var = as.double(var))


  #get the ID's
  #find the middle n_ex IDs to have go extinct
  #use the same tactic as the stabilizing method, but pull the other half that equals n_ex
  #get number of extant species at the end
  n_extant <- nrow(coord) - n_ex

  #now do the splits so that the number of species lost in each tail is stored
  dis_splits <- c(floor(n_extant/2), ceiling(n_extant/2))
  #randomize left and right tail assignments: a bit of randomization if a small bias in numbers on each side exists (extra +/-1)
  dis_sides <- sample(x = dis_splits, size = 2, replace = F)


  #get the LEFT tail to keep
  keep_dis_left <- coord %>%
    arrange(x) %>%
    dplyr::select(id) %>%
    head(id, n = dis_sides[1]) %>%
    unlist(.)


  #get the RIGHT tail to keep
  keep_dis_right <- coord %>%
    arrange(desc(x)) %>%
    dplyr::select(id) %>%
    head(id, n = dis_sides[2]) %>%
    unlist(.)


  #merge tail IDs to same object and remove names
  keep_dis <- c(keep_dis_left, keep_dis_right)
  names(keep_dis) <- NULL


  #now create an object that contains the IDs not in keep_dis
  rm_dis <- coord$id[!coord$id %in% keep_dis]

  #fix the rm_dis to be in the order from lowest abs values to greatest for removal
  rm_dis <- coord %>%
    filter(id %in% rm_dis) %>%
    #mutate(x = abs(x)) %>%
    #arrange(desc(x)) %>%
    #arrange(x) %>%
    mutate(xdist = abs(x - start_1trait$avg)) %>%
    arrange(xdist) %>%
    dplyr::select(id) %>%
    unlist(.)

  names(rm_dis) <- NULL


  ########################################
  #4. stabilizing

  #initialize storage obj for each rep
  sim_1trait_sta <- as_tibble(as.data.frame(matrix(nrow = n_ex + 1, ncol = 3)))
  colnames(sim_1trait_sta) <- c("n_spp", "avg", "var")

  sim_1trait_sta <- sim_1trait_sta %>%
    mutate(n_spp = as.integer(n_spp), avg = as.double(avg), var = as.double(var))

  #get the ID's
  #split the number of extinctions between the two sides and get the extremes
  #we need to have a way to split the two halves
  sta_splits <- c(floor(n_ex/2), ceiling(n_ex/2))

  #randomize the left and right tail assignments (use sta_splits and just sample them twice without replacement)
  sta_sides <- sample(x = sta_splits, size = 2, replace = F)


  #get the ID's to remove: LEFT TAIL
  rm_sta_left <-  coord %>%
    arrange(x) %>%
    dplyr::select(id) %>%
    head(id, n = sta_sides[1]) %>%
    unlist(.)

  #get the ID's to remove: RIGHT TAIL
  rm_sta_right <-  coord %>%
    arrange(desc(x)) %>%
    dplyr::select(id) %>%
    head(id, n = sta_sides[2]) %>%
    unlist(.)

  #merge the tail IDs into one object and remove old names
  rm_sta <- c(rm_sta_left, rm_sta_right)
  names(rm_sta) <- NULL

  #need to randomize the order of the IDs to have go extinct to avoid wonky half/half plot
  #fix the rm_sta to be in the order from greatest distance abs values to lowest for removal
  rm_sta <- coord %>%
    filter(id %in% rm_sta) %>%
    #mutate(x = abs(x)) %>%
    #arrange(desc(x)) %>%
    mutate(xdist = abs(x - start_1trait$avg)) %>%
    arrange(desc(xdist)) %>%
    dplyr::select(id) %>%
    unlist(.)

  names(rm_sta) <- NULL
  #rm_sta <- sample(x = rm_sta, size = length(rm_sta), replace = FALSE)


  ########################################


  #now to get the results populated in the holding objects
  #sliding eval version of adding the values
  #mean
  #random
  sim_1trait_ran$avg <- coord %>%
    arrange(match(id, c(rm_ran, coord$id[!coord$id %in% rm_ran]))) %>%
    slide(.x = .$x, .f = ~mean(.x, na.rm = T), .before = 0, .after = Inf) %>%
    unlist(.) %>%
    .[1:(n_ex + 1)]
  #disruptive
  sim_1trait_dis$avg <- coord %>%
    arrange(match(id, c(rm_dis, coord$id[!coord$id %in% rm_dis]))) %>%
    slide(.x = .$x, .f = ~mean(.x, na.rm = T), .before = 0, .after = Inf) %>%
    unlist(.) %>%
    .[1:(n_ex + 1)]
  #stabilizing
  sim_1trait_sta$avg <- coord %>%
    arrange(match(id, c(rm_sta, coord$id[!coord$id %in% rm_sta]))) %>%
    slide(.x = .$x, .f = ~mean(.x, na.rm = T), .before = 0, .after = Inf) %>%
    unlist(.) %>%
    .[1:(n_ex + 1)]
  #directional
  if(sides == "LEFT" || sides == "RIGHT"){

    #object is the same for either tail
    sim_1trait_dir$avg <- coord %>%
      arrange(match(id, c(rm_dir, coord$id[!coord$id %in% rm_dir]))) %>%
      slide(.x = .$x, .f = ~mean(.x, na.rm = T), .before = 0, .after = Inf) %>%
      unlist(.) %>%
      .[1:(n_ex + 1)]

  } else if(sides == "BOTH"){

    #left tail
    sim_1trait_dirl$avg <- coord %>%
      arrange(match(id, c(rm_dirl, coord$id[!coord$id %in% rm_dirl]))) %>%
      slide(.x = .$x, .f = ~mean(.x, na.rm = T), .before = 0, .after = Inf) %>%
      unlist(.) %>%
      .[1:(n_ex + 1)]

    #right tail
    sim_1trait_dirr$avg <- coord %>%
      arrange(match(id, c(rm_dirr, coord$id[!coord$id %in% rm_dirr]))) %>%
      slide(.x = .$x, .f = ~mean(.x, na.rm = T), .before = 0, .after = Inf) %>%
      unlist(.) %>%
      .[1:(n_ex + 1)]

  }

  #variance
  #random
  sim_1trait_ran$var <- coord %>%
    arrange(match(id, c(rm_ran, coord$id[!coord$id %in% rm_ran]))) %>%
    slide(.x = .$x, .f = ~var(.x), .before = 0, .after = Inf) %>%
    unlist(.) %>%
    .[1:(n_ex + 1)]
  #disruptive
  sim_1trait_dis$var <- coord %>%
    arrange(match(id, c(rm_dis, coord$id[!coord$id %in% rm_dis]))) %>%
    slide(.x = .$x, .f = ~var(.x), .before = 0, .after = Inf) %>%
    unlist(.) %>%
    .[1:(n_ex + 1)]
  #stabilizing
  sim_1trait_sta$var <- coord %>%
    arrange(match(id, c(rm_sta, coord$id[!coord$id %in% rm_sta]))) %>%
    slide(.x = .$x, .f = ~var(.x), .before = 0, .after = Inf) %>%
    unlist(.) %>%
    .[1:(n_ex + 1)]
  #directional
  if(sides == "LEFT" || sides == "RIGHT"){

    #object is the same for either tail
    sim_1trait_dir$var <- coord %>%
      arrange(match(id, c(rm_dir, coord$id[!coord$id %in% rm_dir]))) %>%
      slide(.x = .$x, .f = ~var(.x), .before = 0, .after = Inf) %>%
      unlist(.) %>%
      .[1:(n_ex + 1)]

  } else if(sides == "BOTH"){

    #left tail
    sim_1trait_dirl$var <- coord %>%
      arrange(match(id, c(rm_dirl, coord$id[!coord$id %in% rm_dirl]))) %>%
      slide(.x = .$x, .f = ~var(.x), .before = 0, .after = Inf) %>%
      unlist(.) %>%
      .[1:(n_ex + 1)]

    #right tail
    sim_1trait_dirr$var <- coord %>%
      arrange(match(id, c(rm_dirr, coord$id[!coord$id %in% rm_dirr]))) %>%
      slide(.x = .$x, .f = ~var(.x), .before = 0, .after = Inf) %>%
      unlist(.) %>%
      .[1:(n_ex + 1)]

  }

  #number of species
  sim_1trait_ran$n_spp <- (nrow(coord)):(nrow(coord) - n_ex)
  sim_1trait_dis$n_spp <- (nrow(coord)):(nrow(coord) - n_ex)
  sim_1trait_sta$n_spp <- (nrow(coord)):(nrow(coord) - n_ex)
  if(sides == "LEFT" || sides == "RIGHT"){
    sim_1trait_dir$n_spp <- (nrow(coord)):(nrow(coord) - n_ex)
  } else if(sides == "BOTH"){
    sim_1trait_dirl$n_spp <- (nrow(coord)):(nrow(coord) - n_ex)
    sim_1trait_dirr$n_spp <- (nrow(coord)):(nrow(coord) - n_ex)
  }


  ###########################################

  ###########################################

  #put together object for output: sim_1trait_ran, _dir-l-r, _dis, _sta, _bid-lr-rl

  #add starting data and "type" column so we can just rbind() them all and they can be broken out later by users
  sim_1trait_ran <- sim_1trait_ran %>%
    #rbind(start_1trait, .) %>%
    mutate(type = ("random"))
  sim_1trait_dis <- sim_1trait_dis %>%
    #rbind(start_1trait, .) %>%
    mutate(type = ("disruptive"))
  sim_1trait_sta <- sim_1trait_sta %>%
    #rbind(start_1trait, .) %>%
    mutate(type = ("stabilizing"))

  if(sides == "LEFT"){

    sim_1trait_dir <- sim_1trait_dir %>%
      #rbind(start_1trait, .) %>%
      mutate(type = ("directional_l"))
    #sim_1trait_bid <- sim_1trait_bid %>%
    #  rbind(start_1trait, .) %>%
    #  mutate(type = ("bidirectional_lr"))


  } else if(sides == "RIGHT"){

    sim_1trait_dir <- sim_1trait_dir %>%
      #rbind(start_1trait, .) %>%
      mutate(type = ("directional_r"))
    #sim_1trait_bid <- sim_1trait_bid %>%
    #  rbind(start_1trait, .) %>%
    #  mutate(type = ("bidirectional_rl"))


  } else if(sides == "BOTH"){

    #LEFT TAIL
    sim_1trait_dirl <- sim_1trait_dirl %>%
      #rbind(start_1trait, .) %>%
      mutate(type = ("directional_l"))
    # sim_1trait_bidlr <- sim_1trait_bidlr %>%
    #   rbind(start_1trait, .) %>%
    #   mutate(type = ("bidirectional_lr"))

    #RIGHT TAIL
    sim_1trait_dirr <- sim_1trait_dirr %>%
      #rbind(start_1trait, .) %>%
      mutate(type = ("directional_r"))
    #sim_1trait_bidrl <- sim_1trait_bidrl %>%
    #  rbind(start_1trait, .) %>%
    #  mutate(type = ("bidirectional_rl"))


  }


  #NOTE: WE DO NOT NEED A conditional for RESCALE == FALSE && SLIDING_EVAL == FALSE, this is already what we have at this point!
  #rescale for plotting by the original mean and variance such that the start is at 0
  if(sliding_eval == FALSE && rescale == TRUE){

    sim_1trait_ran <- sim_1trait_ran %>%
      mutate(avg = avg - avg[n_spp == nrow(coord)],
             var = var - var[n_spp == nrow(coord)])
    sim_1trait_dis <- sim_1trait_dis %>%
      mutate(avg = avg - avg[n_spp == nrow(coord)],
             var = var - var[n_spp == nrow(coord)])
    sim_1trait_sta <- sim_1trait_sta %>%
      mutate(avg = avg - avg[n_spp == nrow(coord)],
             var = var - var[n_spp == nrow(coord)])

    if(sides == "LEFT" || sides == "RIGHT"){

      sim_1trait_dir <- sim_1trait_dir %>%
        mutate(avg = avg - avg[n_spp == nrow(coord)],
               var = var - var[n_spp == nrow(coord)])
      #sim_1trait_bid <- sim_1trait_bid %>%
      #  mutate(avg = avg - avg[n_spp == nrow(coord)], var = var - var[n_spp == nrow(coord)])

    } else if(sides == "BOTH"){

      sim_1trait_dirl <- sim_1trait_dirl %>%
        mutate(avg = avg - avg[n_spp == nrow(coord)],
               var = var - var[n_spp == nrow(coord)])
      sim_1trait_dirr <- sim_1trait_dirr %>%
        mutate(avg = avg - avg[n_spp == nrow(coord)],
               var = var - var[n_spp == nrow(coord)])
      #sim_1trait_bidlr <- sim_1trait_bidlr %>%
      #  mutate(avg = avg - avg[n_spp == nrow(coord)], var = var - var[n_spp == nrow(coord)])
      #sim_1trait_bidrl <- sim_1trait_bidrl %>%
      #  mutate(avg = avg - avg[n_spp == nrow(coord)], var = var - var[n_spp == nrow(coord)])

    }

  }
  else if(sliding_eval == TRUE && rescale == TRUE){

    #test region
    #print(sim_1trait_sta %>% slide(.x = .$avg, .f = ~diff(.x), .before = eval_step) %>% ifelse((sapply(., length) == 0), NA, .))
    #print(sim_1trait_sta %>% slide(.x = .$avg, .f = ~diff(.x), .before = eval_step) %>% ifelse((sapply(., length) == 0), NA, .) %>% lapply(X = ., FUN = function(x) mean(x, na.rm = T)) %>% rapply(., f = function(x) ifelse(is.nan(x), 0, x), how = "replace") %>% unlist(.))
    #print(sim_1trait_sta %>% slide(.x = .$avg, .f = ~diff(.x), .before = eval_step) %>% ifelse((sapply(., length) == 0), 0, .) %>% unlist(.))
    #print(sim_1trait_sta %>% slide(.x = .$avg, .f = ~diff(.x), .before = eval_step) %>% ifelse((sapply(., length) == 0), NA, .) %>% lapply(X = ., FUN = function(x) mean(x, na.rm = T)) %>% rapply(., f = function(x) ifelse(is.nan(x), 0, x), how = "replace") %>% unlist(.))


    #implemented region
    sim_1trait_ran <- sim_1trait_ran %>%
      mutate(
        avg = slide(.x = .$avg, .f = ~diff(.x), .before = eval_step) %>%
          ifelse((sapply(., length) == 0), NA, .) %>%
          lapply(X = ., FUN = function(x) mean(x, na.rm = T)) %>%
          rapply(., f = function(x) ifelse(is.nan(x), 0, x), how = "replace") %>%
          unlist(.),
        var = slide(.x = .$var, .f = ~diff(.x), .before = eval_step) %>%
          ifelse((sapply(., length) == 0), NA, .) %>%
          lapply(X = ., FUN = function(x) mean(x, na.rm = T)) %>%
          rapply(., f = function(x) ifelse(is.nan(x), 0, x), how = "replace") %>%
          unlist(.)
      )
    sim_1trait_dis <- sim_1trait_dis %>%
      mutate(
        avg = slide(.x = .$avg, .f = ~diff(.x), .before = eval_step) %>%
          ifelse((sapply(., length) == 0), NA, .) %>%
          lapply(X = ., FUN = function(x) mean(x, na.rm = T)) %>%
          rapply(., f = function(x) ifelse(is.nan(x), 0, x), how = "replace") %>%
          unlist(.),
        var = slide(.x = .$var, .f = ~diff(.x), .before = eval_step) %>%
          ifelse((sapply(., length) == 0), NA, .) %>%
          lapply(X = ., FUN = function(x) mean(x, na.rm = T)) %>%
          rapply(., f = function(x) ifelse(is.nan(x), 0, x), how = "replace") %>%
          unlist(.)
      )
    sim_1trait_sta <- sim_1trait_sta %>%
      mutate(
        avg = slide(.x = .$avg, .f = ~diff(.x), .before = eval_step) %>%
          ifelse((sapply(., length) == 0), NA, .) %>%
          lapply(X = ., FUN = function(x) mean(x, na.rm = T)) %>%
          rapply(., f = function(x) ifelse(is.nan(x), 0, x), how = "replace") %>%
          unlist(.),
        var = slide(.x = .$var, .f = ~diff(.x), .before = eval_step) %>%
          ifelse((sapply(., length) == 0), NA, .) %>%
          lapply(X = ., FUN = function(x) mean(x, na.rm = T)) %>%
          rapply(., f = function(x) ifelse(is.nan(x), 0, x), how = "replace") %>%
          unlist(.)
      )
    if(sides == "LEFT" || sides == "RIGHT"){
      sim_1trait_dir <- sim_1trait_dir %>%
        mutate(
          avg = slide(.x = .$avg, .f = ~diff(.x), .before = eval_step) %>%
            ifelse((sapply(., length) == 0), NA, .) %>%
            lapply(X = ., FUN = function(x) mean(x, na.rm = T)) %>%
            rapply(., f = function(x) ifelse(is.nan(x), 0, x), how = "replace") %>%
            unlist(.),
          var = slide(.x = .$var, .f = ~diff(.x), .before = eval_step) %>%
            ifelse((sapply(., length) == 0), NA, .) %>%
            lapply(X = ., FUN = function(x) mean(x, na.rm = T)) %>%
            rapply(., f = function(x) ifelse(is.nan(x), 0, x), how = "replace") %>%
            unlist(.)
        )

    } else if(sides == "BOTH"){
      sim_1trait_dirl <- sim_1trait_dirl %>%
        mutate(
          avg = slide(.x = .$avg, .f = ~diff(.x), .before = eval_step) %>%
            ifelse((sapply(., length) == 0), NA, .) %>%
            lapply(X = ., FUN = function(x) mean(x, na.rm = T)) %>%
            rapply(., f = function(x) ifelse(is.nan(x), 0, x), how = "replace") %>%
            unlist(.),
          var = slide(.x = .$var, .f = ~diff(.x), .before = eval_step) %>%
            ifelse((sapply(., length) == 0), NA, .) %>%
            lapply(X = ., FUN = function(x) mean(x, na.rm = T)) %>%
            rapply(., f = function(x) ifelse(is.nan(x), 0, x), how = "replace") %>%
            unlist(.)
        )
      sim_1trait_dirr <- sim_1trait_dirr %>%
        mutate(
          avg = slide(.x = .$avg, .f = ~diff(.x), .before = eval_step) %>%
            ifelse((sapply(., length) == 0), NA, .) %>%
            lapply(X = ., FUN = function(x) mean(x, na.rm = T)) %>%
            rapply(., f = function(x) ifelse(is.nan(x), 0, x), how = "replace") %>%
            unlist(.),
          var = slide(.x = .$var, .f = ~diff(.x), .before = eval_step) %>%
            ifelse((sapply(., length) == 0), NA, .) %>%
            lapply(X = ., FUN = function(x) mean(x, na.rm = T)) %>%
            rapply(., f = function(x) ifelse(is.nan(x), 0, x), how = "replace") %>%
            unlist(.)
        )

    }

  }
  else if(sliding_eval == TRUE && rescale == FALSE){

    #test region
    #print(sim_1trait_sta %>%
    #          slide(.x = .$avg, .f = ~mean(.x, na.rm = T), .before = eval_step) %>%
    #            unlist(.)
    #      )
    #print(sim_1trait_sta %>%
    #        slide(.x = .$avg, .f = ~mean(.x, na.rm = T), .before = 0) %>%
    #        unlist(.)
    #)
    #print(sim_1trait_sta %>%
    #        slide(.x = .$avg, .f = ~mean(.x, na.rm = T), .before = 1) %>%
    #        unlist(.)
    #)

    #implemented region
    sim_1trait_ran <- sim_1trait_ran %>%
      mutate(
        avg = slide(.x = .$avg, .f = ~mean(.x, na.rm = T), .before = eval_step) %>%
          unlist(.),
        var = slide(.x = .$var, .f = ~mean(.x, na.rm = T), .before = eval_step) %>%
          unlist(.)
      )
    sim_1trait_dis <- sim_1trait_dis %>%
      mutate(
        avg = slide(.x = .$avg, .f = ~mean(.x, na.rm = T), .before = eval_step) %>%
          unlist(.),
        var = slide(.x = .$var, .f = ~mean(.x, na.rm = T), .before = eval_step) %>%
          unlist(.)
      )
    sim_1trait_sta <- sim_1trait_sta %>%
      mutate(
        avg = slide(.x = .$avg, .f = ~mean(.x, na.rm = T), .before = eval_step) %>%
          unlist(.),
        var = slide(.x = .$var, .f = ~mean(.x, na.rm = T), .before = eval_step) %>%
          unlist(.)
      )

    if(sides == "LEFT" || sides == "RIGHT"){
      sim_1trait_dirl <- sim_1trait_dirl %>%
        mutate(
          avg = slide(.x = .$avg, .f = ~mean(.x, na.rm = T), .before = eval_step) %>%
            unlist(.),
          var = slide(.x = .$var, .f = ~mean(.x, na.rm = T), .before = eval_step) %>%
            unlist(.)
        )

    } else if(sides == "BOTH"){
      sim_1trait_dirl <- sim_1trait_dirl %>%
        mutate(
          avg = slide(.x = .$avg, .f = ~mean(.x, na.rm = T), .before = eval_step) %>%
            unlist(.),
          var = slide(.x = .$var, .f = ~mean(.x, na.rm = T), .before = eval_step) %>%
            unlist(.)
        )
      sim_1trait_dirr <- sim_1trait_dirr %>%
        mutate(
          avg = slide(.x = .$avg, .f = ~mean(.x, na.rm = T), .before = eval_step) %>%
            unlist(.),
          var = slide(.x = .$var, .f = ~mean(.x, na.rm = T), .before = eval_step) %>%
            unlist(.)
        )

    }

  }

  #bind them together
  if(sides == "LEFT" || sides == "RIGHT"){

    #sim_1trait_final <- rbind(sim_1trait_ran, sim_1trait_dis, sim_1trait_sta, sim_1trait_bid, sim_1trait_dir)
    sim_1trait_final <- rbind(sim_1trait_ran, sim_1trait_dis, sim_1trait_sta, sim_1trait_dir)

  } else if(sides == "BOTH"){

    #sim_1trait_final <- rbind(sim_1trait_ran, sim_1trait_dis, sim_1trait_sta, sim_1trait_bidlr, sim_1trait_bidrl, sim_1trait_dirl, sim_1trait_dirr)
    sim_1trait_final <- rbind(sim_1trait_ran, sim_1trait_dis, sim_1trait_sta, sim_1trait_dirl, sim_1trait_dirr)

  }

  #output
  return(sim_1trait_final)

}
