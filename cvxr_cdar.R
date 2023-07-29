library(CVXR)
library(parma)

data(dji30ret, package = "rugarch", overwrite = T)
.x <- as.matrix(dji30ret[, 1:4])


cdar_problem <- function(.x, type = c("minrisk","optimal"),
                       expected_reward = rep(0, ncol(.x)),
                       lower = rep(0, ncol(.x)),
                       upper = rep(1, ncol(.x)),
                       target = 0, budget = NULL, leverage = NULL,
                       solver = NULL, verbose = FALSE, alpha = 0.05, ...)
{
  m <- NCOL(.x)
  n <- NROW(.x)
  
  A <- diag(n)
  for (i in 1:n-1){
    A[i+1, i] = -1
  }
  
  if (type == "minrisk") {
    w <- Variable(m)
    z1 <- Variable(n)
    v1 <- Variable(1)
    u1 <- Variable(n)
    if (is.null(budget) & !is.null(leverage)) {
      # exact 1-norm problem
      w_positive <- Variable(m)
      w_negative <- Variable(m)
      z <- rep(1, m)
#      z <- Variable(m, boolean = TRUE)
      ones <- rep(1, m)
      constraints <- list(w_positive >= 0,
                          w_negative >= 0,
                          w == w_positive - w_negative,
                          w_positive <= leverage * z,
                          w_negative <= leverage * (ones - z),
                          sum(w_positive) + sum(w_negative) == leverage,
                          sum(expected_reward * w) >= target,
                          w >= lower,
                          w <= upper,
                          z1 - u1 + v1 >= 0,
                          z1 >= 0,
                          u1 >= 0,
                          .x %*% w + A %*% u1 >= 0)
    } else if (!is.null(budget) & !is.null(leverage)) {
      constraints <- list(sum(w) == budget,
                          norm1(w) <= leverage,
                          sum(expected_reward * w) >= target,
                          w >= lower,
                          w <= upper,
                          z1 - u1 + v1 >= 0,
                          z1 >= 0,
                          u1 >= 0,
                          .x %*% w + A %*% u1 >= 0)
    } else if (!is.null(budget) & is.null(leverage)) {
      constraints <- list(sum(w) == budget,
                          sum(expected_reward * w) >= target,
                          w >= lower,
                          w <= upper,
                          z1 - u1 + v1 >= 0,
                          z1 >= 0,
                          u1 >= 0,
                          .x %*% w + A %*% u1 >= 0)
    } else {
      stop("\neither budget or leverage or both must be non NULL!")
    }
  } else {
    w <- Variable(m)
    v <- Variable(1, pos = TRUE)
    z1 <- Variable(n)
    v1 <- Variable(1)
    u1 <- Variable(n)
    if (is.null(budget) & !is.null(leverage)) {
      # exact 1-norm problem
      w_positive <- Variable(m)
      w_negative <- Variable(m)
      z <- Variable(m, boolean = TRUE)
      ones <- as.integer(rep(1, m))
      constraints <- list(w == w_positive - w_negative,
                          w_positive >= 0,
                          w_negative >= 0,
                          w_positive <= leverage * z,
                          w_negative <= leverage * (ones - z),
                          sum(w_positive) + sum(w_negative) == v * leverage,
                          sum(expected_reward * w) == 1,
                          w >= v * lower,
                          w <= v * upper,
                          z1 - u1 + v1 >= 0,
                          z1 >= 0,
                          u1 >= 0,
                          .x %*% w + A %*% u1 >= 0)
    } else if (!is.null(budget) & !is.null(leverage)) {
      constraints <- list(sum(expected_reward * w) == 1,
                          sum(w) == v * budget,
                          norm1(w) <= v * leverage,
                          w >= v * lower,
                          w <= v * upper,
                          z1 - u1 + v1 >= 0,
                          z1 >= 0,
                          u1 >= 0,
                          .x %*% w + A %*% u1 >= 0)
    } else if (!is.null(budget) & is.null(leverage)) {
      constraints <- list(sum(w) == v * budget,
                          sum(expected_reward * w) == 1,
                          w >= v * lower,
                          w <= v * upper,
                          z1 - u1 + v1 >= 0,
                          z1 >= 0,
                          u1 >= 0,
                          .x %*% w + A %*% u1 >= 0)
    } else {
      stop("\neither budget or leverage or both must be non NULL!")
    }
  }
  admissible_solver <- c("CBC", "OSQP", "GLPK_MI", "SCS", "ECOS_BB", "GUROBI")
  risk <- v1 + 1/(n * alpha) * sum(z1)
  problem <- Problem(Minimize(risk), constraints = constraints)
  if (is.null(solver)) solver <- admissible_solver[1]
  solution <- psolve(problem, solver = solver, verbose = verbose, ...)
  if (type == "optimal") {
    w <- solution$getValue(w)/solution$getValue(v)
  } else {
    w <- solution$getValue(w)
  }
  w <- as.numeric(w)
  port_risk <- riskfun(as.matrix(w), .x, risk = "cdar")
  port_reward <- sum(expected_reward * w)
  names(w) <- colnames(.x)
  return(list(weights = w, risk = port_risk,
              reward = port_reward,
              risk_type = "cdar",
              problem_type = type))
}



# optimal solution
cvxr_solution <- cdar_problem(.x, type = c("minrisk","optimal")[2],
                                 expected_reward = colMeans(.x) * 10000,
                                 lower = rep(-0.2, ncol(.x)),
                                 upper = rep(0.4, ncol(.x)),
                                 budget = 1, solver = "CBC")


spec <- parmaspec(scenario = .x, forecast = colMeans(.x),
                  risk = "cdar", riskType = "optimal",
                  LB = rep(-0.2, 4), UB = rep(0.4, 4), options = list(alpha = 0.05))
parma_solution <- parmasolve(spec)

round(cbind(cvxr_solution$weights, weights(parma_solution)), 4)



# minrisk solution
cvxr_solution <- cdar_problem(.x, type = c("minrisk","optimal")[1], 
                                 expected_reward = colMeans(.x) * 10000,
                                 lower = rep(-0.8, ncol(.x)), 
                                 upper = rep(0.8, ncol(.x)), 
                                 target = 2,
                                 budget = 1, solver = "CBC")


spec <- parmaspec(scenario = .x, forecast = colMeans(.x) * 10000, 
                  risk = "cdar", riskType = "minrisk",
                  target = 2, targetType = "inequality",
                  LB = rep(-0.8, 4), UB = rep(0.8, 4), options = list(alpha = 0.05))
parma_solution <- parmasolve(spec)

round(cbind(cvxr_solution$weights, weights(parma_solution)), 4)                                



# minrisk solution with no leverage
cvxr_solution <- cdar_problem(.x, type = c("minrisk","optimal")[1], 
                                 expected_reward = colMeans(.x) * 10000,
                                 lower = rep(0, ncol(.x)), 
                                 upper = rep(0.8, ncol(.x)), 
                                 target = 2,
                                 budget = 1, leverage = NULL, solver = "CBC")


spec <- parmaspec(scenario = .x, forecast = colMeans(.x) * 10000, 
                  risk = "cdar", riskType = "minrisk",
                  target = 2, targetType = "inequality",
                  LB = rep(0, 4), UB = rep(0.8, 4), budget = 1)
parma_solution <- parmasolve(spec)

round(cbind(cvxr_solution$weights, weights(parma_solution)), 4)



# optimal solution with no leverage
cvxr_solution <- cdar_problem(.x, type = c("minrisk","optimal")[2], 
                                 expected_reward = colMeans(.x) * 10000,
                                 lower = rep(0, ncol(.x)), 
                                 upper = rep(0.8, ncol(.x)), 
                                 budget = 1, leverage = NULL, solver = "CBC")


spec <- parmaspec(scenario = .x, forecast = colMeans(.x) * 10000, 
                  risk = "cdar", riskType = "optimal",
                  LB = rep(0, 4), UB = rep(0.8, 4), budget = 1, options = list(alpha = 0.05))
parma_solution <- parmasolve(spec)

round(cbind(cvxr_solution$weights, weights(parma_solution)), 4)
