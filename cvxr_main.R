main_problem <- function(.x, type = c("minrisk","optimal"),
                         risk = c("ev", "mad", "minimax", "lpm", "es", "cdar"),
                         expected_reward = rep(0, ncol(.x)),
                         lower = rep(0, ncol(.x)),
                         upper = rep(1, ncol(.x)),
                         target = 0, budget = NULL, leverage = NULL,
                         solver = NULL, verbose = FALSE, options = list(alpha = 0.05, threshold = 999, moment = 1), ...)
{
  if (risk == "ev"){
    cvxr_solution <- ev_problem(.x, type = type, 
                                expected_reward = expected_reward,
                                lower = lower, 
                                upper = upper, 
                                target = target,
                                budget = budget, leverage = leverage, 
                                solver = solver, verbose = verbose)     
  }
  
  if (risk == "mad"){
    cvxr_solution <- mad_problem(.x, type = type, 
                                 expected_reward = expected_reward,
                                 lower = lower, 
                                 upper = upper,
                                 target = target,
                                 budget = budget, leverage = leverage, 
                                 solver = solver, verbose = verbose)
  } 
  
  if (risk == "minimax"){
    cvxr_solution <- minimax_problem(.x, type = type, 
                                     expected_reward = expected_reward,
                                     lower = lower, 
                                     upper = upper,
                                     target = target,
                                     budget = budget, leverage = leverage, 
                                     solver = solver, verbose = verbose)
  }
  
  if (risk == "lpm"){
    cvxr_solution <- lpm_problem(.x, type = type, 
                                 expected_reward = expected_reward,
                                 lower = lower, 
                                 upper = upper, 
                                 target = target,
                                 budget = budget, leverage = leverage, 
                                 solver = solver, verbose = verbose, 
                                 moment = options$moment, threshold = options$threshold)
  }

  if (risk == "es"){
    cvxr_solution <- es_problem(.x, type = type,
                           expected_reward = expected_reward,
                           lower = lower,
                           upper = upper,
                           target = target, 
                           budget = budget, leverage = leverage,
                           solver = solver, verbose = verbose, alpha = options$alpha)
  }  
  
  if (risk == "cdar"){
    cvxr_solution <- cdar_problem(.x, type = type,
                             expected_reward = expected_reward,
                             lower = lower,
                             upper = upper,
                             target = target, budget = budget, leverage = leverage,
                             solver = solver, verbose = verbose, alpha = options$alpha)
  }
  
  return(round(cvxr_solution$weights, 4))
}


