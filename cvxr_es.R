# expected shortfall
# for the optimal problem need to
# scale the expected reward as expected_reward/min(abs(expected_reward))


data(dji30ret, package = "rugarch", overwrite = T)
.x <- as.matrix(dji30ret[, 1:4])

es_problem <- function(.x, type = c("minrisk","optimal"),
                       expected_reward = rep(0, ncol(.x)),
                       lower = rep(0, ncol(.x)),
                       upper = rep(1, ncol(.x)),
                       target = 0, budget = NULL, leverage = NULL,
                       solver = NULL, verbose = FALSE, alpha = 0.05, ...)
{
    m <- NCOL(.x)
    n <- NROW(.x)
    if (type == "minrisk") {
        w <- Variable(m)
        q <- Variable(1)
        d <- Variable(n, pos = TRUE)
        if (is.null(budget) & !is.null(leverage)) {
            # exact 1-norm problem
            w_positive <- Variable(m)
            w_negative <- Variable(m)
            z <- rep(1, ncol(.x))
#            z <- Variable(m, boolean = TRUE)
            ones <- rep(1, m)
            constraints <- list(w_positive >= 0,
                                w_negative >= 0,
                                w == w_positive - w_negative,
                                w_positive <= leverage * z,
                                w_negative <= leverage * (ones - z),
                                sum(w_positive) + sum(w_negative) == leverage,
                                sum(expected_reward * w) >= target,
                                (.x %*% w + q) >= -d,
                                w >= lower,
                                w <= upper)
        } else if (!is.null(budget) & !is.null(leverage)) {
            constraints <- list(sum(w) == budget,
                                norm1(w) <= leverage,
                                sum(expected_reward * w) >= target,
                                (.x %*% w + q) >= -d,
                                w >= lower,
                                w <= upper)
        } else if (!is.null(budget) & is.null(leverage)) {
            constraints <- list(sum(w) == budget,
                                sum(expected_reward * w) >= target,
                                (.x %*% w + q) >= -d,
                                w >= lower,
                                w <= upper)
        } else {
            stop("\neither budget or leverage or both must be non NULL!")
        }
    } else {
        w <- Variable(m)
        v <- Variable(1, pos = TRUE)
        q <- Variable(1)
        d <- Variable(n, pos = TRUE)
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
                                (.x %*% w + q) >= -d,
                                w >= v * lower,
                                w <= v * upper)
        } else if (!is.null(budget) & !is.null(leverage)) {
            constraints <- list(sum(expected_reward * w) == 1,
                                sum(w) == v * budget,
                                norm1(w) <= v * leverage,
                                (.x %*% w + q) >= -d,
                                w >= v * lower,
                                w <= v * upper)
        } else if (!is.null(budget) & is.null(leverage)) {
            constraints <- list(sum(w) == v * budget,
                                sum(expected_reward * w) == 1,
                                (.x %*% w + q) >= -d,
                                w >= v * lower,
                                w <= v * upper)
        } else {
            stop("\neither budget or leverage or both must be non NULL!")
        }
    }
    admissible_solver <- c("CBC","GLPK_MI")
    risk <- (1/(n * alpha) * sum(d)) + q
    problem <- Problem(Minimize(risk), constraints = constraints)
    if (is.null(solver)) solver <- admissible_solver[1]
    solution <- psolve(problem, solver = solver, verbose = verbose, ...)
    if (type == "optimal") {
        w <- solution$getValue(w)/solution$getValue(v)
        value_at_risk <- solution$getValue(q)/solution$getValue(v)
    } else {
        w <- solution$getValue(w)
        value_at_risk <- solution$getValue(q)
    }
    w <- as.numeric(w)
    port <- .x %*% w
    port_risk <- riskfun(w, .x, risk = "es", alpha = alpha)
    port_reward <- sum(expected_reward * w)
    names(w) <- colnames(.x)
    return(list(weights = w, risk = port_risk,
                value_at_risk = value_at_risk,
                reward = port_reward,
                risk_type = "es",
                problem_type = type))
}


# optimal solution with leverage
cvxr_solution <- es_problem(.x, type = c("minrisk","optimal")[2], 
                            expected_reward = colMeans(.x) * 10000,
                            lower = rep(-0.2, ncol(.x)), 
                            upper = rep(0.4, ncol(.x)), 
                            budget = NULL, leverage = 0.5, solver = "CBC", alpha = 0.05)


spec <- parmaspec(scenario = .x, forecast = colMeans(.x), 
                  risk = "CVaR", riskType = "optimal",
                  LB = rep(-0.2,4), UB = rep(0.4, 4), leverage = 0.5, alpha = 0.05)
parma_solution <- parmasolve(spec)

round(cbind(cvxr_solution$weights, weights(parma_solution)), 4)



# minrisk solution with leverage
cvxr_solution <- es_problem(.x, type = c("minrisk","optimal")[1], 
                            expected_reward = colMeans(.x) * 10000,
                            lower = rep(-0.8, ncol(.x)), 
                            upper = rep(0.8, ncol(.x)), 
                            target = 2,
                            budget = NULL, leverage = 1, solver = "GUROBI", alpha = 0.05)


spec <- parmaspec(scenario = .x, forecast = colMeans(.x) * 10000, 
                  risk = "CVaR", riskType = "minrisk",
                  target = 2, targetType = "inequality",
                  LB = rep(-0.8,4), UB = rep(0.8, 4), leverage = 1, alpha = 0.05)
parma_solution <- parmasolve(spec)

round(cbind(cvxr_solution$weights, weights(parma_solution)), 4)



# minrisk solution with no leverage
cvxr_solution <- es_problem(.x, type = c("minrisk","optimal")[1], 
                            expected_reward = colMeans(.x) * 10000,
                            lower = rep(0, ncol(.x)), 
                            upper = rep(0.8, ncol(.x)), 
                            target = 2,
                            budget = 1, leverage = NULL, solver = "CBC", alpha = 0.05)


spec <- parmaspec(scenario = .x, forecast = colMeans(.x) * 10000, 
                  risk = "CVaR", riskType = "minrisk",
                  target = 2, targetType = "inequality",
                  LB = rep(0, 4), UB = rep(0.8, 4), budget = 1, alpha = 0.05)
parma_solution <- parmasolve(spec)

round(cbind(cvxr_solution$weights, weights(parma_solution)), 4)



# optimal solution with no leverage
cvxr_solution <- es_problem(.x, type = c("minrisk","optimal")[2], 
                            expected_reward = colMeans(.x) * 10000,
                            lower = rep(0, ncol(.x)), 
                            upper = rep(0.8, ncol(.x)), 
                            budget = 1, leverage = NULL, solver = "CBC", alpha = 0.05)


spec <- parmaspec(scenario = .x, forecast = colMeans(.x) * 10000, 
                  risk = "CVaR", riskType = "optimal",
                  LB = rep(0, 4), UB = rep(0.8, 4), budget = 1, alpha = 0.05)
parma_solution <- parmasolve(spec)

round(cbind(cvxr_solution$weights, weights(parma_solution)), 4)

